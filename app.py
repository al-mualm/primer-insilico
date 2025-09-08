# -*- coding: utf-8 -*-
# In-silico PCR with UCSC first (JSON API → mirrors), Arabic UI, gel, FASTA download
# If UCSC is temporarily protected by Cloudflare, we auto-rotate mirrors and detect blocks.

import io, re, time, math, hashlib
from typing import List, Dict
import requests
import streamlit as st
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup

# ---------------- Page & styles ----------------
st.set_page_config(page_title="محاكاة PCR", layout="wide")
st.markdown("""
<style>
html, body, [class*="css"] { direction: rtl; font-family: "Noto Naskh Arabic","Tajawal","Cairo",sans-serif; }
.header { background: linear-gradient(135deg,#ecf3ff 0%,#f0fff6 100%); border-radius:18px; padding:16px 20px; border:1px solid #e6eef7; margin-bottom:14px;}
.header h1 { margin:0 0 6px 0; font-size:1.6rem; }
.card { border:1px solid #e9eef4; background:#fff; border-radius:14px; padding:14px 16px; margin:10px 0;}
.small { font-size:.9rem; color:#475569;}
pre.mono { background:#f8fafc; border:1px solid #e9eef4; padding:8px 10px; border-radius:10px; overflow-x:auto; direction:ltr; text-align:left;}
.badges span { display:inline-block; margin:4px 6px 0 0; padding:6px 10px; border-radius:999px; border:1px solid #e3e8ef; background:#fff; font-size:.9rem; }
hr.soft { border:none; border-top:1px dashed #e3e8ef; margin:8px 0 12px 0;}
.warn { background:#fff7ed; border:1px solid #fed7aa; padding:8px 12px; border-radius:10px; }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="header">
  <h1>محاكاة PCR — UCSC سريعة + جيل افتراضي + تنزيل FASTA</h1>
  <div class="small">
    صُمِّم بواسطة <b>Mahmood Al-Mualm</b> — <b>محمود أحمد محي المعلّم</b> ·
    البريد: <a href="mailto:mahmoodalmoalm@gmail.com">mahmoodalmoalm@gmail.com</a> ·
    الهاتف/واتساب: <a href="tel:+9647730585329">+964 7730585329</a>
  </div>
  <div class="badges"><span>UCSC In-silico PCR</span><span>In-silico Gel</span><span>واجهة عربية</span></div>
</div>
""", unsafe_allow_html=True)

# ---------------- UCSC config ----------------
UCSC_JSON = "https://api.genome.ucsc.edu/hgPcr"   # JSON/FASTA gateway used by UCSC’s API front
UCSC_HTML = [
    "https://genome.ucsc.edu/cgi-bin/hgPcr",
    "https://genome-euro.ucsc.edu/cgi-bin/hgPcr",
    "https://genome-asia.ucsc.edu/cgi-bin/hgPcr",
]

# label → (org, db)
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

UA = {"User-Agent": "Primer-Insilico/1.2 (+mahmoodalmoalm@gmail.com)"}

def _clean(seq: str) -> str:
    return re.sub(r"[^ACGTNacgtn]", "", (seq or "")).upper()

def _looks_like_cloudflare(text: str) -> bool:
    # Cloudflare bot/turnstile page fingerprints
    t = text.lower()
    return ("turnstile" in t) or ("cloudflare" in t) or ("challenge" in t and "cf" in t)

@st.cache_data(ttl=24*60*60, show_spinner=False)
def _http_get_cached(url: str, params: Dict, timeout: int) -> requests.Response:
    r = requests.get(url, params=params, timeout=timeout, headers=UA)
    r.raise_for_status()
    return r

# ------------- UCSC JSON first -------------
def ucsc_via_json(fwd: str, rev: str, org: str, db: str, max_bp: int, timeout: int = 12) -> List[Dict]:
    """
    Try UCSC API endpoint (returns JSON or plain FASTA).
    If JSON: we expect a structure with PCR products; otherwise we may get FASTA in text.
    """
    params = dict(org=org, db=db, wp_f=fwd, wp_r=rev, wp_size=int(max_bp))
    r = _http_get_cached(UCSC_JSON, params, timeout)
    ct = r.headers.get("Content-Type","").lower()
    text = r.text

    # If API returns FASTA directly, parse it.
    if "text/plain" in ct or text.startswith(">"):
        return parse_fasta_products(text)

    # If JSON, try to parse common shapes
    try:
        data = r.json()
    except Exception:
        # not JSON; maybe HTML (blocked)
        if _looks_like_cloudflare(text):
            raise RuntimeError("واجهنا حماية Cloudflare على UCSC API.")
        # try to parse from <pre> if present
        return parse_html_products(text)

    # A couple of known shapes (UCSC sometimes returns nested fields)
    # We normalize to a list of {chrom, start, end, size, strand, sequence}
    products = []

    # Guess 1: 'results' or 'pcr' list with dict items
    for key in ("results","pcr","items","products"):
        if key in data and isinstance(data[key], list):
            for item in data[key]:
                prod = _coerce_item(item)
                if prod: products.append(prod)

    # Guess 2: maybe it's a dict keyed by product names
    if not products and isinstance(data, dict):
        for k,v in data.items():
            if isinstance(v, dict) or isinstance(v, list):
                prod = _coerce_item(v)
                if prod:
                    products.append(prod)

    return products

def _coerce_item(item) -> Dict:
    """
    Try to convert any UCSC-like product into our normalized dict.
    Accepts dicts with fields: chrom/chr, start, end, size, strand, seq/sequence
    """
    if isinstance(item, list) and item and isinstance(item[0], dict):
        item = item[0]

    if not isinstance(item, dict):
        return {}

    chrom = item.get("chrom") or item.get("chr") or item.get("target") or ""
    start = item.get("start") or item.get("txStart") or item.get("s") or None
    end   = item.get("end")   or item.get("txEnd")   or item.get("e") or None
    strand= item.get("strand") or item.get("dir") or "+"
    seq   = item.get("sequence") or item.get("seq") or ""

    # If only seq exists and header is embedded, try to pull coords out of header.
    if seq and not (chrom and start and end):
        m = re.search(r">(chr[\w\.\-]+):(\d+)-(\d+)", seq.splitlines()[0])
        if m:
            chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))

    if chrom and start and end:
        try:
            start, end = int(start), int(end)
            size = abs(end - start) + 1
            return {"chrom": chrom, "start": start, "end": end, "size": size, "strand": strand, "sequence": (seq or "").replace("\r","")}
        except Exception:
            return {}
    return {}

# ------------- UCSC mirrors (HTML/FASTA) -------------
def ucsc_via_html(fwd: str, rev: str, org: str, db: str, max_bp: int, timeout: int = 12) -> List[Dict]:
    params = {
        "org": org,
        "db": db,
        "wp_target": "genome",
        "wp_f": fwd,
        "wp_r": rev,
        "wp_size": int(max_bp),
        "Submit": "submit"
    }
    last_err = None
    for url in UCSC_HTML:
        try:
            r = _http_get_cached(url, params, timeout)
            txt = r.text
            if _looks_like_cloudflare(txt):
                # try next mirror
                last_err = RuntimeError("تم حظر الطلب بواسطة Cloudflare على هذا المرآة.")
                continue

            # Prefer parsing FASTA blocks inside <pre>, fallback to plain pre text
            soup = BeautifulSoup(txt, "lxml")
            pre_txt = "\n".join(p.get_text("\n") for p in soup.find_all("pre")) or txt
            products = parse_fasta_products(pre_txt)
            if not products:
                products = parse_html_products(pre_txt)
            if products:
                return products
            last_err = RuntimeError("لم أتعرف على نتائج صالحة من هذا المرآة.")
        except Exception as e:
            last_err = e
            time.sleep(0.8)
    if last_err:
        raise last_err
    return []

# ------------- Parsers -------------
def parse_fasta_products(text: str) -> List[Dict]:
    """
    Parse UCSC PCR output when it includes FASTA-like products:
    >chr1:100-200(+)
    ACTG...
    """
    products = []
    blocks = re.split(r"(?m)^>", text)
    for blk in blocks:
        blk = blk.strip()
        if not blk:
            continue
        header, *seq_lines = blk.splitlines()
        m = re.search(r"(chr[\w\.\-]+):(\d+)-(\d+)(?:\(([-+])\))?", header)
        if not m:
            # Accept headers like "chr1:100-200" without strand
            m = re.search(r"(chr[\w\.\-]+):(\d+)-(\d+)", header)
        if not m:
            continue
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        strand = m.group(4) if m.lastindex and m.group(m.lastindex) in ["+","-"] else "+"
        # keep only ACGTN lines
        seq = "".join(s.strip() for s in seq_lines if re.fullmatch(r"[ACGTNacgtn]+", s.strip()))
        size = abs(end - start) + 1
        products.append({"chrom": chrom, "start": start, "end": end, "size": size, "strand": strand, "sequence": seq.upper()})
    return products

def parse_html_products(text: str) -> List[Dict]:
    """
    As a last resort, parse coordinates from preformatted text without FASTA.
    """
    products = []
    # lines like: chr1:12345-12456 (+)
    for line in text.splitlines():
        m = re.search(r"(chr[\w\.\-]+):(\d+)-(\d+).*?([+-])?", line)
        if not m:
            continue
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        strand = m.group(4) if m.group(4) in ["+","-"] else "+"
        size = abs(end - start) + 1
        products.append({"chrom": chrom, "start": start, "end": end, "size": size, "strand": strand, "sequence": ""})
    return products

# ------------- Gel drawing -------------
def _gel_y(bp, a=100.0, b=50.0):  # log spacing
    return a - b * math.log10(max(bp, 1))

def render_gel(sizes: List[int], ladder="auto") -> io.BytesIO:
    lad100 = list(range(100, 1600, 100))
    lad1k  = [250,500,750,1000,1500,2000,3000,4000,5000,6000,8000,10000]
    mx = max(sizes) if sizes else 0
    if ladder == "auto":
        ladder = "100bp" if mx <= 800 else "1kb"
    lad = lad100 if ladder == "100bp" else lad1k

    def y(v): return [_gel_y(x) for x in v]
    all_y = y(lad + sizes)
    ymin, ymax = min(all_y)-5, max(all_y)+5

    fig, ax = plt.subplots(figsize=(5.2, 6), dpi=170)
    lane_w, gap = 1.0, 2.0
    x_lad = 0.0
    x_s   = x_lad + lane_w + gap

    for yy, sz in zip(y(lad), lad):
        ax.hlines(yy, x_lad, x_lad+lane_w, linewidth=3)
        ax.text(x_lad-0.25, yy, f"{sz}", va="center", ha="right", fontsize=8)

    for yy, sz in zip(y(sizes), sizes):
        ax.hlines(yy, x_s, x_s+lane_w, linewidth=4)
        ax.text(x_s+lane_w+0.2, yy, f"{sz} bp", va="center", fontsize=8)

    ax.text(x_lad+lane_w/2, ymin+2, "Ladder", ha="center", fontsize=9)
    ax.text(x_s+lane_w/2,    ymin+2, "Sample", ha="center", fontsize=9)
    ax.set_xlim(-1.2, x_s+lane_w+1.2)
    ax.set_ylim(ymin, ymax)
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf

# ------------- UI -------------
with st.sidebar:
    st.subheader("إعدادات المرجع")
    genome_label = st.selectbox("المرجع (UCSC)", list(UCSC_GENOMES.keys()), index=0)
    org, db = UCSC_GENOMES[genome_label]
    min_bp  = st.number_input("أصغر حجم (bp) (للعرض فقط)", 50, 2000, 80)
    max_bp  = st.number_input("أكبر حجم (bp)", 50, 10000, 4000)
    st.caption("ملاحظة: أفضل سرعة مع USC ≥ 4000 bp.")

st.markdown("## (Primers) إدخال البادئات")
c1, c2 = st.columns(2)
with c1:
    rev = st.text_input("Reverse primer", "")
with c2:
    fwd = st.text_input("Forward primer", "")

run = st.button("تشغيل المحاكاة (UCSC)")

if run:
    fwd = _clean(fwd)
    rev = _clean(rev)
    if not fwd or not rev:
        st.error("يرجى إدخال كلا البادئين.")
        st.stop()
    if len(fwd) < 16 or len(rev) < 16:
        st.markdown('<div class="warn">تحذير: يفضّل أن تكون البادئات ≥ 16 nt للحصول على نتائج جيدة.</div>', unsafe_allow_html=True)

    all_hits: List[Dict] = []
    with st.status("جارٍ الاستعلام عبر UCSC (API أولًا ثم المرايا) …", expanded=False) as s:
        # 1) JSON API
        try:
            s.update(label="الاتصال بـ UCSC API (JSON)…")
            hits = ucsc_via_json(fwd, rev, org, db, max_bp, timeout=12)
            all_hits.extend(hits or [])
        except Exception as e:
            st.write("تعذّر استخدام واجهة JSON:", str(e))

        # 2) Mirrors if needed
        if not all_hits:
            s.update(label="تجربة مرايا UCSC (نص/فاستا)…")
            try:
                hits = ucsc_via_html(fwd, rev, org, db, max_bp, timeout=12)
                all_hits.extend(hits or [])
            except Exception as e:
                st.write("المرايا HTML/FASTA لم تنجح:", str(e))

        if not all_hits:
            st.error("لم يتم العثور على نواتج ضمن الحد الأقصى للحجم.\n\nجرب واحدًا مما يلي:\n• تقليل الحد الأقصى للحجم\n• استخدام مرجع آخر (hg19 بدلًا من hg38 …)\n• إعادة المحاولة بعد دقائق (إن كان هناك حماية مؤقتة)")
            st.stop()

    # keep top 5 by size (smallest first)
    all_hits = sorted(all_hits, key=lambda h: h.get("size", 10**9))[:5]

    st.markdown("## النتائج (أفضل ٥)")
    for i, h in enumerate(all_hits, 1):
        seq = (h.get("sequence") or "").strip()
        seq_preview = seq[:120] + ("…" if len(seq) > 120 else "")
        st.markdown(
            f"""<div class="card">
<b>Product {i}</b> — {h.get('chrom','?')}:{h.get('start','?'):,}-{h.get('end','?'):,}
(الحجم: <b>{h.get('size','?')} bp</b>, الاتجاه: {h.get('strand','+')})
<hr class="soft"/>
<pre class="mono">{seq_preview or '(لا توجد تسلسلات في هذه المخرجات)'} </pre>
</div>""",
            unsafe_allow_html=True
        )

    sizes = [h["size"] for h in all_hits if "size" in h]
    if sizes:
        st.markdown("### الجل الافتراضي")
        st.image(render_gel(sorted(sizes)), use_column_width=True)

    # FASTA download
    fasta_lines = []
    for i, h in enumerate(all_hits, 1):
        hdr = f">{db}|{h.get('chrom','?')}:{h.get('start','?')}-{h.get('end','?')}|size={h.get('size','?')}bp|prod{i}\n"
        fasta_lines.append(hdr)
        seq = (h.get("sequence") or "").strip().replace("\r","")
        if seq:
            for j in range(0, len(seq), 70):
                fasta_lines.append(seq[j:j+70] + "\n")
    st.download_button(
        "تنزيل (FASTA)",
        data=("".join(fasta_lines)).encode(),
        file_name=f"ucsc_{db}_amplicons.fasta",
        mime="text/plain"
    )

