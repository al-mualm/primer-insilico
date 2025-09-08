# -*- coding: utf-8 -*-
# Streamlit: Fast in-silico PCR (UCSC + Local BLAST), Arabic UI

import io, os, re, time, math, sqlite3, shutil, hashlib, subprocess
from pathlib import Path
from typing import List, Dict
from collections import defaultdict

import requests
import streamlit as st
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup  # for UCSC HTML parsing

# -------------------- Page & styles --------------------
st.set_page_config(page_title="محاكاة PCR", layout="wide")
st.markdown("""
<style>
html, body, [class*="css"] { direction: rtl; font-family: "Noto Naskh Arabic","Tajawal","Cairo",sans-serif; }
.header { background: linear-gradient(135deg,#ecf3ff 0%,#f0fff6 100%); border-radius:18px; padding:16px 20px; border:1px solid #e6eef7; margin-bottom:14px;}
.header h1 { margin:0 0 6px 0; font-size:1.6rem; }
.badges span { display:inline-block; margin:4px 6px 0 0; padding:6px 10px; border-radius:999px; border:1px solid #e3e8ef; background:#fff; font-size:.9rem; }
.card { border:1px solid #e9eef4; background:#fff; border-radius:14px; padding:14px 16px; margin:10px 0;}
.small { font-size:.9rem; color:#475569;}
pre.mono { background:#f8fafc; border:1px solid #e9eef4; padding:8px 10px; border-radius:10px; overflow-x:auto; direction:ltr; text-align:left;}
hr.soft { border:none; border-top:1px dashed #e3e8ef; margin:8px 0 12px 0;}
.review { border-bottom:1px dashed #e3e8ef; padding:8px 0; }
.warn { background:#fff7ed; border:1px solid #fed7aa; padding:8px 12px; border-radius:10px; }
</style>
""", unsafe_allow_html=True)

# -------------------- Tiny SQLite DB (visitors + reviews) --------------------
DB_PATH = Path("app_data.db").resolve()

def _init_db():
    with sqlite3.connect(DB_PATH) as c:
        c.execute("""CREATE TABLE IF NOT EXISTS visitors(
                        id INTEGER PRIMARY KEY,
                        ts DATETIME DEFAULT CURRENT_TIMESTAMP)""")
        c.execute("""CREATE TABLE IF NOT EXISTS reviews(
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        name TEXT, rating INTEGER, comment TEXT,
                        ts DATETIME DEFAULT CURRENT_TIMESTAMP)""")
        c.commit()

def _inc_visit_once():
    if "vis_once" not in st.session_state:
        with sqlite3.connect(DB_PATH) as c:
            c.execute("INSERT INTO visitors DEFAULT VALUES")
            c.commit()
        st.session_state["vis_once"] = True

def _visitors():
    with sqlite3.connect(DB_PATH) as c:
        return c.execute("SELECT COUNT(*) FROM visitors").fetchone()[0]

def add_review(name, rating, comment):
    with sqlite3.connect(DB_PATH) as c:
        c.execute("INSERT INTO reviews(name,rating,comment) VALUES (?,?,?)",
                  (name, rating, comment))
        c.commit()

def list_reviews(n=50):
    with sqlite3.connect(DB_PATH) as c:
        return c.execute("SELECT name,rating,comment,ts FROM reviews ORDER BY id DESC LIMIT ?",
                         (n,)).fetchall()

_init_db(); _inc_visit_once()

# -------------------- Header --------------------
st.markdown(f"""
<div class="header">
  <h1>محاكاة PCR السريعة</h1>
  <div class="small">
    صُمِّم هذا الموقع بواسطة <b>Mahmood Al-Mualm</b> — <b>محمود أحمد محي المعلّم</b>.<br/>
    البريد: <a href="mailto:mahmoodalmoalm@gmail.com">mahmoodalmoalm@gmail.com</a> — الهاتف/واتساب: <a href="tel:+9647730585329">+964 7730585329</a>
  </div>
  <div class="badges">
    <span>UCSC In-silico PCR</span><span>BLAST+ محلي (اختياري)</span><span>In-silico Gel</span><span>واجهة عربية</span><span>الزوار: {_visitors()}</span>
  </div>
</div>
""", unsafe_allow_html=True)

# -------------------- Genomes --------------------
# UCSC engines (label -> (org, db))
UCSC_GENOMES = {
    "Human (hg38)": ("Human", "hg38"),
    "Human (hg19)": ("Human", "hg19"),
    "Mouse (mm39)": ("Mouse", "mm39"),
    "Rat (rn7)": ("Rat", "rn7"),
    "Zebrafish (danRer11)": ("Zebrafish", "danRer11"),
    "Fruit fly (dm6)": ("D. melanogaster", "dm6"),
    "Worm (ce11)": ("C. elegans", "ce11"),
    "Yeast (sacCer3)": ("S. cerevisiae", "sacCer3"),
}

# Local BLAST DB paths (edit these when you set up the VM)
# You can also add "Bacteria panel": "/data/blastdb/bacteria_panel"
LOCAL_DBS = {
    "Human (hg38)": "/data/blastdb/hg38",
    "Mouse (mm39)": "/data/blastdb/mm39",
    "Rat (rn7)": "/data/blastdb/rn7",
}

# -------------------- Hardened UCSC client --------------------
UCSC_ENDPOINTS = [
    "https://genome.ucsc.edu/cgi-bin/hgPcr",
    "https://genome-euro.ucsc.edu/cgi-bin/hgPcr",
    "https://genome-asia.ucsc.edu/cgi-bin/hgPcr",
]

def _clean_primer(seq: str) -> str:
    s = re.sub(r"[^ACGTNacgtn]", "", seq or "")
    return s.upper()

@st.cache_data(ttl=24*60*60, show_spinner=False)
def _cached_ucsc(cache_key: str, endpoint: str, params: dict, timeout: int) -> str:
    r = requests.get(
        endpoint, params=params, timeout=timeout,
        headers={"User-Agent": "Primer-Insilico/1.0 (+contact: mahmoodalmoalm@gmail.com)"}
    )
    r.raise_for_status()
    return r.text

def _request_ucsc_any(endpoints: List[str], params: dict, timeout: int = 12) -> str:
    last_err = None
    for ep in endpoints:
        for attempt in range(3):
            try:
                ck = hashlib.sha1((ep + "|" + str(sorted(params.items()))).encode()).hexdigest()
                return _cached_ucsc(ck, ep, params, timeout)
            except Exception as e:
                last_err = e
                time.sleep(1.2 * (attempt + 1))
    raise RuntimeError(f"تعذّر الاتصال بخادم UCSC بعد عدة محاولات. آخر خطأ: {last_err}")

def ucsc_in_silico_pcr(fwd: str, rev: str, org: str, db: str,
                       max_bp: int = 4000, top_n: int = 5, timeout: int = 12) -> List[Dict]:
    fwd = _clean_primer(fwd); rev = _clean_primer(rev)
    if not fwd or not rev:
        raise ValueError("البادئات فارغة بعد التنظيف. استخدم A/C/G/T/N فقط.")
    max_bp = max(50, min(int(max_bp or 4000), 4000))

    params = {"org": org, "db": db, "wp_target": "genome",
              "wp_f": fwd, "wp_r": rev, "wp_size": max_bp, "Submit": "submit"}

    html = _request_ucsc_any(UCSC_ENDPOINTS, params, timeout=timeout)
    soup = BeautifulSoup(html, "lxml")
    pre_txt = "\n".join(pre.get_text("\n") for pre in soup.find_all("pre")) or soup.get_text("\n")

    hits = []
    blocks = re.split(r"(?m)^>", pre_txt)
    for blk in blocks:
        blk = blk.strip()
        if not blk:
            continue
        header, *seq_lines = blk.splitlines()
        m = re.search(r"(chr[\w\.]+):(\d+)-(\d+).*?([+-])?", header)
        if not m:
            continue
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        strand = m.group(4) if m.group(4) in ["+", "-"] else "+"
        seq = "".join(s.strip() for s in seq_lines if set(s.strip().upper()) <= set("ACGTN")).upper()
        size = abs(end - start) + 1
        hits.append({"chrom": chrom, "start": start, "end": end, "size": size, "strand": strand, "sequence": seq})

    hits.sort(key=lambda h: h["size"])
    return hits[:max(1, int(top_n or 5))]

# -------------------- Local BLAST+ engine --------------------
def _check_tool(name: str) -> bool:
    return shutil.which(name) is not None

def blast_hits(primer_seq: str, db_path: str, max_targets: int = 2000, timeout: int = 20):
    cmd = [
        "blastn", "-task", "blastn-short",
        "-db", db_path,
        "-query", "-",  # stdin
        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send sstrand evalue bitscore",
        "-max_target_seqs", str(max_targets),
        "-evalue", "1000",
        "-dust", "no", "-soft_masking", "false",
        "-word_size", "7"
    ]
    q = f">q\n{primer_seq.strip().upper()}\n".encode()
    res = subprocess.run(cmd, input=q, capture_output=True, timeout=timeout, check=True)
    hits = []
    for line in res.stdout.decode().strip().splitlines():
        cols = line.split("\t")
        d = {
            "sseqid": cols[1],
            "pident": float(cols[2]),
            "length": int(cols[3]),
            "sstart": int(cols[6]),
            "send":   int(cols[7]),
            "sstrand": cols[8],
        }
        hits.append(d)
    return hits

def find_amplicons_local(fwd: str, rev: str, db_path: str,
                         max_bp: int = 4000, min_match: int = 16, top_n: int = 5) -> List[Dict]:
    fwd = _clean_primer(fwd); rev = _clean_primer(rev)
    f_hits = [h for h in blast_hits(fwd, db_path) if h["length"] >= min_match]
    r_hits = [h for h in blast_hits(rev, db_path) if h["length"] >= min_match]

    amplicons = []
    for fh in f_hits:
        for rh in r_hits:
            if fh["sseqid"] != rh["sseqid"]:
                continue
            # fwd on plus, rev on minus
            if fh["sstrand"] == "plus" and rh["sstrand"] == "minus":
                start = min(fh["sstart"], fh["send"])
                end   = max(rh["sstart"], rh["send"])
                size  = end - start + 1
                if 0 < size <= max_bp and end >= start:
                    amplicons.append({"chrom": fh["sseqid"], "start": start, "end": end, "size": size})
    # unique & sort
    seen = set(); uniq = []
    for a in sorted(amplicons, key=lambda x: x["size"]):
        key = (a["chrom"], a["start"], a["end"])
        if key not in seen:
            seen.add(key); uniq.append(a)
    return uniq[:max(1, top_n)]

def fetch_seq_local(contig: str, start: int, end: int, db_path: str) -> str:
    cmd = ["blastdbcmd", "-db", db_path, "-entry", contig, "-range", f"{start}-{end}"]
    return subprocess.check_output(cmd).decode()

# -------------------- Gel plotting --------------------
def _log(bp, a=100.0, b=50.0):
    return a - b * math.log10(max(bp, 1))

def gel(lanes: Dict[str, List[int]], ladder="auto"):
    # ladder defs
    lad100 = list(range(100, 1600, 100))
    lad1k = [250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000]
    mx = max([max(v) for v in lanes.values() if v] or [0])
    ladder = "100bp" if (ladder == "auto" and mx <= 800) else ("1kb" if ladder == "auto" else ladder)
    lad = lad100 if ladder == "100bp" else lad1k

    def y(vals): return [_log(x) for x in vals]
    all_y = y(lad + [s for arr in lanes.values() for s in arr])
    ymin, ymax = min(all_y) - 5, max(all_y) + 5

    fig, ax = plt.subplots(figsize=(max(5.2, 1.5*len(lanes)), 6), dpi=170)
    lane_w, gap = 1.0, 1.5
    xs = []
    x = 0.0
    # ladder lane
    xs.append(("Ladder", x))
    for yy, sz in zip(y(lad), lad): ax.hlines(yy, x, x+lane_w, linewidth=3); ax.text(x-0.25, yy, f"{sz}", va="center", ha="right", fontsize=8)
    x += lane_w + gap

    # sample lanes
    for name, sizes in lanes.items():
        xs.append((name, x))
        for yy, sz in zip(y(sizes), sizes):
            ax.hlines(yy, x, x+lane_w, linewidth=4)
            ax.text(x + lane_w + 0.2, yy, f"{sz} bp", va="center", fontsize=8)
        x += lane_w + gap

    for name, xv in xs:
        ax.text(xv + lane_w/2, ymin+2, name, ha="center", fontsize=9)

    ax.set_xlim(-1.2, x + 0.5)
    ax.set_ylim(ymin, ymax)
    ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="png", bbox_inches="tight"); plt.close(fig); buf.seek(0)
    return buf

# -------------------- UI --------------------
page = st.sidebar.radio("الصفحة", ["المحاكاة", "التقييمات"], index=0)
engine = st.sidebar.radio("المحرك", ["UCSC (سريع ومجاني)", "BLAST+ المحلي"], index=0)
ucsc_timeout = st.sidebar.slider("مهلة الاتصال بـ UCSC (ثوانٍ)", 6, 20, 12)

if page == "المحاكاة":
    st.sidebar.markdown("### إعدادات المرجع")
    genome_label = st.selectbox("المرجع (Genome)", list(UCSC_GENOMES.keys()))
    org, db = UCSC_GENOMES[genome_label]

    pmin = st.number_input("أصغر حجم (bp) (للعرض فقط)", 20, 50000, 80, 10)
    pmax = st.number_input("أكبر حجم (bp) — حد البحث", 50, 50000, 4000, 50)
    st.caption("ملاحظة: في UCSC يُفضَّل ≤ 4000 bp لسرعة أعلى.")

    st.markdown("### إدخال البادئات (Primers)")
    c1, c2 = st.columns(2)
    with c1: fwd = st.text_input("Forward Primer", "")
    with c2: rev = st.text_input("Reverse Primer", "")

    if st.button("تشغيل المحاكاة", type="primary"):
        if not fwd or not rev:
            st.error("يرجى إدخال كلا البادئين."); st.stop()
        if len(fwd) < 16 or len(rev) < 16:
            st.markdown('<div class="warn">تحذير: البادئات أقصر من 16 نيوكليوتيد — قد تكون النتائج غير محددة.</div>', unsafe_allow_html=True)

        if engine == "UCSC (سريع ومجاني)":
            with st.status("جارٍ الاتصال بخوادم UCSC…", expanded=True) as status:
                try:
                    hits = ucsc_in_silico_pcr(fwd, rev, org, db, max_bp=pmax, top_n=5, timeout=ucsc_timeout)
                    status.update(label="تم الاستلام من UCSC", state="complete")
                except Exception as e:
                    status.update(label="فشل الاتصال بكل مرايا UCSC.", state="error")
                    st.error("تعذّر الوصول إلى UCSC الآن. أعد المحاولة لاحقًا أو استخدم محرك BLAST+ المحلي إن كان مُتاحًا.")
                    st.exception(e); st.stop()

            if not hits:
                st.info("لم يتم العثور على نواتج ضمن الحد الأقصى للحجم."); st.stop()

            # Report
            st.markdown("## التقارير المفصّلة (أفضل ٥)")
            for i, h in enumerate(hits, 1):
                st.markdown(f"""
                <div class="card">
                  <b>Product {i}</b> — {h['chrom']}:{h['start']:,}-{h['end']:,}
                  (الحجم: <b>{h['size']} bp</b>, الاتجاه: {h.get('strand','?')})
                  <hr class="soft"/>
                  <pre class="mono">{h['sequence'][:80]}{"..." if len(h['sequence'])>80 else ""}</pre>
                </div>
                """, unsafe_allow_html=True)

            # Gel: single lane with multiple bands
            sizes = sorted([h["size"] for h in hits])
            st.markdown("### الجل الافتراضي")
            st.image(gel({"Sample": sizes}), use_column_width=True)

            # FASTA
            fasta = []
            for i, h in enumerate(hits, 1):
                fasta.append(f">{db}|{h['chrom']}:{h['start']}-{h['end']}|size={h['size']}bp|prod{i}\n")
                s = h["sequence"]
                for j in range(0, len(s), 70): fasta.append(s[j:j+70] + "\n")
            st.download_button("تنزيل نواتج الأمبليكونات (FASTA)",
                               data=("".join(fasta)).encode(),
                               file_name=f"pcr_{db}.fasta", mime="text/plain")

        else:
            # Local BLAST+ engine
            db_path = LOCAL_DBS.get(genome_label, "")
            if not db_path:
                st.error("لا يوجد مسار قاعدة بيانات محلية لهذا المرجع. حدّث LOCAL_DBS في الكود."); st.stop()
            if not Path(db_path + ".nin").exists() and not Path(db_path + ".00.nin").exists():
                st.error(f"ملفات قاعدة BLAST غير موجودة: {db_path}. شغّل makeblastdb أولاً."); st.stop()
            if not (_check_tool("blastn") and _check_tool("blastdbcmd")):
                st.error("أدوات BLAST+ غير مُثبّتة على الخادم (blastn/blastdbcmd). استخدم UCSC أو ثبّت BLAST+."); st.stop()

            with st.status("جارٍ تشغيل BLAST+ محليًّا…", expanded=True) as status:
                try:
                    amps = find_amplicons_local(fwd, rev, db_path, max_bp=pmax, min_match=16, top_n=5)
                    status.update(label="تم إيجاد مواضع البادئات", state="complete")
                except Exception as e:
                    status.update(label="حدث خطأ أثناء BLAST", state="error"); st.exception(e); st.stop()

            if not amps:
                st.info("لا توجد أمبليكونات ضمن الحد الأقصى للحجم."); st.stop()

            # fetch sequences & prepare lanes
            hits = []
            for i, a in enumerate(amps, 1):
                try:
                    fa = fetch_seq_local(a["chrom"], a["start"], a["end"], db_path)
                    seq = "".join(x.strip() for x in fa.splitlines() if not x.startswith(">"))
                except Exception:
                    seq = ""
                hits.append({"chrom": a["chrom"], "start": a["start"], "end": a["end"], "size": a["size"], "sequence": seq})

            st.markdown("## التقارير المفصّلة (أفضل ٥)")
            for i, h in enumerate(hits, 1):
                st.markdown(f"""
                <div class="card">
                  <b>Product {i}</b> — {h['chrom']}:{h['start']:,}-{h['end']:,}
                  (الحجم: <b>{h['size']} bp</b>)
                  <hr class="soft"/>
                  <pre class="mono">{(h['sequence'] or '')[:80]}{ '...' if (h['sequence'] and len(h['sequence'])>80) else ''}</pre>
                </div>
                """, unsafe_allow_html=True)

            # Gel:
            # If you later let users pick multiple organisms (multiple DBs), pass a dict of lanes.
            sizes = sorted([h["size"] for h in hits])
            st.markdown("### الجل الافتراضي")
            st.image(gel({"Sample": sizes}), use_column_width=True)

            # FASTA
            fasta = []
            for i, h in enumerate(hits, 1):
                fasta.append(f">{db_path}|{h['chrom']}:{h['start']}-{h['end']}|size={h['size']}bp|prod{i}\n")
                s = h["sequence"]
                for j in range(0, len(s), 70): fasta.append(s[j:j+70] + "\n")
            st.download_button("تنزيل نواتج الأمبليكونات (FASTA)",
                               data=("".join(fasta)).encode(),
                               file_name=f"pcr_local.fasta", mime="text/plain")

elif page == "التقييمات":
    st.markdown("## آراء المستخدمين")
    with st.form("rev"):
        name = st.text_input("الاسم (اختياري)", "")
        rating = st.slider("التقييم", 1, 5, 5)
        comment = st.text_area("تعليقك", "", height=120)
        ok = st.form_submit_button("إرسال")
    if ok:
        if not comment.strip():
            st.warning("يرجى كتابة تعليق.")
        else:
            add_review(name.strip() or "مستخدم", rating, comment.strip())
            st.success("شكرًا!"); st.experimental_rerun()

    st.markdown("### أحدث المراجعات")
    for nm, rt, cm, ts in list_reviews(100):
        st.markdown(
            f"<div class='review'><b>{nm}</b> — ⭐️ {rt}/5 "
            f"<span class='small'>({ts})</span><br/>{cm}</div>",
            unsafe_allow_html=True
        )

