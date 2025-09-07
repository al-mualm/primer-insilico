# -*- coding: utf-8 -*-
import io
import time
import math
import textwrap
import requests
import xml.etree.ElementTree as ET
from collections import defaultdict

import streamlit as st
import matplotlib.pyplot as plt

# ============ إعداد واجهة التطبيق ============

st.set_page_config(
    page_title="محاكاة تفاعل PCR (عربي)",
    page_icon="🧬",
    layout="wide"
)

# بسيط: تحسين المظهر والخط واتجاه RTL
st.markdown("""
<style>
/* الخط والاتجاه */
html, body, [class*="css"]  {
  direction: rtl;
  font-family: "Noto Naskh Arabic", "Tajawal", "Cairo", "Helvetica", sans-serif;
}

/* صندوق العنوان */
.hero {
  background: linear-gradient(135deg, #f5f7ff 0%, #eef9ff 100%);
  border-radius: 18px;
  padding: 18px 22px;
  border: 1px solid #e7eef7;
  margin-bottom: 18px;
}

/* بطاقات النتائج */
.card {
  border: 1px solid #e9eef4;
  background: #ffffff;
  border-radius: 14px;
  padding: 14px 16px;
  margin: 10px 0;
}

/* تذييل */
.footer {
  margin-top: 20px;
  padding: 12px 16px;
  border-top: 1px dashed #e3e8ef;
  color: #334155;
  font-size: 0.94rem;
}
.small {
  font-size: 0.86rem;
  color: #475569;
}
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hero">
  <h2>🧬 محاكاة <span style="white-space:nowrap">PCR</span> بالعربية</h2>
  <div class="small">
    أدخل بادئًا أماميًا (Forward) وبادئًا عكسيًا (Reverse)، سنستخدم BLAST في NCBI للعثور على مواقع التطابق، 
    ربط الضربات لاستخراج الأمبليكونات، ثم سنعرض <b>جل افتراضي</b>.
    <br/>
    إذا ظهرت النتائج من <b>نفس الـaccession</b> (نفس الكائن/العزلة) ⇒ <b>حارة واحدة</b> بعدة أشرطة. <br/>
    إذا كانت النتائج من <b>accession</b>ات مختلفة (أنواع/عزلات مختلفة) ⇒ <b>عدة حارات</b> (حارة لكل accession).
  </div>
</div>
""", unsafe_allow_html=True)

# ============ ثوابت وخصائص NCBI ============

NCBI_TOOL  = "primer-insilico-ar"
NCBI_EMAIL = "mahmoodalmoalm@gmail.com"   # اجعلها بريدك للاحترام لدى NCBI

# ============ دوال مساعدة ============

def to_fasta(seq, title="query"):
    seq = "".join([c for c in seq.strip().upper() if c in "ACGTURYKMSWBDHVN"])
    wrapped = "\n".join(textwrap.wrap(seq, width=70))
    return f">{title}\n{wrapped}\n"

def blast_put_and_wait(fasta_query, database="refseq_genomic", task="blastn-short",
                       program="blastn", sleep_s=3.0, max_wait_s=120):
    """إرسال إلى BLAST (خدمة NCBI) وانتظار النتيجة (XML)."""
    put_params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": fasta_query,
        "TASK": task,
        "FORMAT_TYPE": "XML",
        "MEGABLAST": "on",
        "FILTER": "L",
        "UNGAPPED_ALIGNMENT": "on",
        "TOOL": NCBI_TOOL,
        "EMAIL": NCBI_EMAIL,
    }
    r = requests.post("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data=put_params, timeout=60)
    r.raise_for_status()
    rid = None
    for line in r.text.splitlines():
        if "RID =" in line:
            rid = line.split("=", 1)[1].strip()
            break
    if not rid:
        raise RuntimeError("تعذّر الحصول على RID من BLAST.")

    start_t = time.time()
    while True:
        g = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi",
                         params={"CMD":"Get","RID":rid,"FORMAT_TYPE":"XML"}, timeout=60)
        g.raise_for_status()
        txt = g.text
        if "Status=WAITING" in txt:
            if time.time() - start_t > max_wait_s:
                raise TimeoutError("انتهى وقت الانتظار لنتيجة BLAST.")
            time.sleep(sleep_s)
            continue
        if "Status=FAILED" in txt:
            raise RuntimeError("فشل تشغيل BLAST على خادم NCBI.")
        if "Status=UNKNOWN" in txt:
            raise RuntimeError("RID غير معروف في BLAST.")
        # في حال لا توجد ضربات:
        if ("QBlastInfoBegin" in txt and "QBlastInfoEnd" in txt and
           "ThereAreHits=yes" not in txt and "ThereAreHits=Yes" not in txt and "<Hit>" not in txt):
            return None
        return txt

def parse_blast_hits(blast_xml_text):
    """
    يعيد قائمة ضربات: [(accession, defline, hit_from, hit_to, strand, evalue, pct_id)]
    defline تحتوي غالباً على وصف يذكر النوع/العزلة.
    """
    if not blast_xml_text:
        return []
    try:
        root = ET.fromstring(blast_xml_text)
    except ET.ParseError:
        idx = blast_xml_text.find("<")
        if idx >= 0:
            root = ET.fromstring(blast_xml_text[idx:])
        else:
            return []

    hits = []
    for hit in root.iterfind(".//Hit"):
        acc = hit.findtext("Hit_accession") or ""
        dfl = hit.findtext("Hit_def") or ""     # مهم لوسم الكائن
        for hsp in hit.iterfind("Hit_hsps/Hsp"):
            hit_from = int(hsp.findtext("Hsp_hit-from"))
            hit_to   = int(hsp.findtext("Hsp_hit-to"))
            strand   = "+" if hit_from < hit_to else "-"
            evalue   = float(hsp.findtext("Hsp_evalue"))
            identity = int(hsp.findtext("Hsp_identity"))
            alignlen = int(hsp.findtext("Hsp_align-len"))
            pct_id   = 100.0 * identity / max(1, alignlen)
            hits.append((acc, dfl, hit_from, hit_to, strand, evalue, pct_id))
    return hits

def pair_hits(forward_hits, reverse_hits, product_min=80, product_max=2000):
    """
    أزواج (fwd, rev) على نفس الaccession باتجاهين متعاكسين ضمن نطاق الحجم.
    تُرجع قائمة أمبليكونات: dict فيه accession/defline/الحجم وحدود الضربتين.
    """
    by_acc_rev = defaultdict(list)
    for acc, dfl, hfrom, hto, strand, evalue, pct in reverse_hits:
        by_acc_rev[acc].append((dfl, hfrom, hto, strand, evalue, pct))

    amplicons = []
    for acc, dfl, ffrom, fto, fstrand, fev, fpct in forward_hits:
        if acc not in by_acc_rev:
            continue
        for rdfl, rfrom, rto, rstrand, rev, rpct in by_acc_rev[acc]:
            if fstrand == rstrand:
                continue
            f_end = max(ffrom, fto)   # نهاية ضربة + عادة أكبر
            r_end = min(rfrom, rto)   # الطرف الأقرب في الضربة -
            size = abs(r_end - f_end) + 1
            if product_min <= size <= product_max:
                amplicons.append({
                    "accession": acc,
                    "defline": dfl or rdfl,
                    "product_size": int(size),
                    "fwd": {"start": min(ffrom, fto), "end": max(ffrom, fto), "strand": fstrand,
                            "evalue": fev, "pct_identity": fpct},
                    "rev": {"start": min(rfrom, rto), "end": max(rfrom, rto), "strand": rstrand,
                            "evalue": rev, "pct_identity": rpct},
                })
    amplicons.sort(key=lambda x: (x["accession"], x["product_size"]))
    return amplicons

def efetch_sequence(accession, start, end, rettype="fasta"):
    r = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params={
            "db": "nuccore",
            "id": accession,
            "rettype": rettype,
            "retmode": "text",
            "seq_start": start,
            "seq_stop": end,
            "strand": 1,
            "tool": NCBI_TOOL,
            "email": NCBI_EMAIL,
        },
        timeout=60
    )
    r.raise_for_status()
    return r.text

# ---- رسم الجل (بدون ألوان مخصصة) ----

def log_migration(bp, a=100.0, b=50.0):
    return a - b * math.log10(bp)

def draw_single_lane_gel(product_sizes, ladder="auto"):
    """حارة واحدة (مجتمع لكل الأمبليكونات)."""
    ladder_100 = list(range(100, 1600, 100))
    ladder_1k  = [250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000]
    if ladder == "auto":
        max_size = max(product_sizes) if product_sizes else 0
        ladder = "100bp" if max_size <= 800 else "1kb"
    ladder_sizes = ladder_100 if ladder == "100bp" else ladder_1k

    def y_of(s): return [log_migration(x) for x in s]
    y_all = y_of(ladder_sizes + product_sizes) if product_sizes else y_of(ladder_sizes)
    y_min, y_max = min(y_all) - 5, max(y_all) + 5

    fig, ax = plt.subplots(figsize=(5.2, 6), dpi=170)
    lane_w, gap = 1.0, 2.0
    lane0_x = 0.0
    lane1_x = lane0_x + lane_w + gap

    # Ladder
    for yy, sz in zip(y_of(ladder_sizes), ladder_sizes):
        ax.hlines(yy, lane0_x, lane0_x + lane_w, linewidth=3)
        ax.text(lane0_x - 0.25, yy, f"{sz}", va="center", ha="right", fontsize=8)
    # Sample
    for yy, sz in zip(y_of(product_sizes), product_sizes):
        ax.hlines(yy, lane1_x, lane1_x + lane_w, linewidth=4)
        ax.text(lane1_x + lane_w + 0.2, yy, f"{sz} bp", va="center", fontsize=8)

    ax.set_xlim(-1.2, lane1_x + lane_w + 1.2)
    ax.set_ylim(y_min, y_max); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(lane0_x + lane_w/2, y_min + 2, "Ladder", ha="center", fontsize=9)
    ax.text(lane1_x + lane_w/2, y_min + 2, "Sample", ha="center", fontsize=9)
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="png", bbox_inches="tight"); plt.close(fig)
    buf.seek(0); return buf

def draw_multi_lane_gel(lanes, lane_labels, ladder="auto"):
    """
    lanes: list[list[int]]  (أحجام الأمبليكونات لكل حارة)
    lane_labels: list[str]  (اسم accession/الوسم لكل حارة)
    """
    ladder_100 = list(range(100, 1600, 100))
    ladder_1k  = [250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000]

    all_sizes = [s for lane in lanes for s in lane]
    max_size = max(all_sizes) if all_sizes else 0
    if ladder == "auto":
        ladder = "100bp" if max_size <= 800 else "1kb"
    ladder_sizes = ladder_100 if ladder == "100bp" else ladder_1k

    def y_of(s): return [log_migration(x) for x in s]
    y_all = y_of(ladder_sizes + all_sizes) if all_sizes else y_of(ladder_sizes)
    y_min, y_max = min(y_all) - 5, max(y_all) + 5

    fig, ax = plt.subplots(figsize=(6.0 + 1.2*len(lanes), 6), dpi=170)
    lane_w, gap = 1.0, 1.6
    # Ladder lane at 0
    lane0_x = 0.0
    for yy, sz in zip(y_of(ladder_sizes), ladder_sizes):
        ax.hlines(yy, lane0_x, lane0_x + lane_w, linewidth=3)
        ax.text(lane0_x - 0.25, yy, f"{sz}", va="center", ha="right", fontsize=8)

    # Sample lanes
    for i, lst in enumerate(lanes, start=1):
        x = lane0_x + i*(lane_w + gap)
        for yy, sz in zip(y_of(lst), lst):
            ax.hlines(yy, x, x + lane_w, linewidth=4)
            ax.text(x + lane_w + 0.15, yy, f"{sz} bp", va="center", fontsize=8)
        label = lane_labels[i-1] if i-1 < len(lane_labels) else f"Sample {i}"
        ax.text(x + lane_w/2, y_min + 2, label, ha="center", fontsize=9)

    ax.text(lane0_x + lane_w/2, y_min + 2, "Ladder", ha="center", fontsize=9)
    right_edge = lane0_x + (len(lanes)+1)*(lane_w + gap)
    ax.set_xlim(-1.2, right_edge + 0.8)
    ax.set_ylim(y_min, y_max); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="png", bbox_inches="tight"); plt.close(fig)
    buf.seek(0); return buf

# ============ واجهة الإدخال ============

with st.sidebar:
    st.markdown("### إعدادات البحث")
    database = st.selectbox("قاعدة بيانات BLAST", ["refseq_genomic", "refseq_representative_genomes", "nt"])
    product_min = st.number_input("أصغر حجم (bp)", min_value=20,  max_value=50000, value=80,   step=10)
    product_max = st.number_input("أكبر حجم (bp)", min_value=20,  max_value=50000, value=2000, step=10)
    st.caption("للاستخدام الكثيف يفضّل BLAST محلي.")

st.markdown("### إدخال البادئات (Primers)")
colA, colB = st.columns(2)
with colA:
    fwd = st.text_input("Forward Primer", placeholder="مثال: ACGTACGTACGTACGTACGT")
with colB:
    rev = st.text_input("Reverse Primer", placeholder="مثال: TGCATGCATGCATGCATGCA")

if st.button("تشغيل المحاكاة", type="primary"):
    if not fwd or not rev:
        st.error("يرجى إدخال كلا البادئين.")
        st.stop()

    with st.status("جارٍ الاتصال بـ NCBI BLAST وتشغيل الاستعلامات…", expanded=True) as status:
        st.write("🔹 إرسال الاستعلام للبريمير الأمامي…")
        fwd_xml = blast_put_and_wait(to_fasta(fwd, "forward"), database=database)
        time.sleep(1.2)
        st.write("🔹 إرسال الاستعلام للبريمير العكسي…")
        rev_xml = blast_put_and_wait(to_fasta(rev, "reverse"), database=database)
        status.update(label="تم استلام النتائج", state="complete")

    fwd_hits = parse_blast_hits(fwd_xml)
    rev_hits = parse_blast_hits(rev_xml)

    if not fwd_hits or not rev_hits:
        st.warning("لم يتم العثور على ضربات كافية لربط أمبليكون. جرّب قاعدة بيانات أخرى أو بادئات مختلفة.")
        st.stop()

    amplicons = pair_hits(fwd_hits, rev_hits, product_min=product_min, product_max=product_max)
    if not amplicons:
        st.info("لا توجد أزواج ضمن نطاق الحجم المحدد.")
        st.stop()

    # عرض النتائج في بطاقات
    st.markdown("### النتائج المحتملة")
    for i, amp in enumerate(amplicons, 1):
        with st.container():
            st.markdown(f"""
            <div class="card">
            <b>{i}. {amp['accession']}</b> — <b>{amp['product_size']} bp</b><br/>
            <span class="small">{amp['defline']}</span><br/>
            <span class="small">Fwd: {amp['fwd']['start']}-{amp['fwd']['end']} ({amp['fwd']['strand']}) | Rev: {amp['rev']['start']}-{amp['rev']['end']} ({amp['rev']['strand']})</span>
            </div>
            """, unsafe_allow_html=True)

    # 1) هل كل الأمبليكونات من نفس accession ؟ (حارة واحدة متعددة الأشرطة)
    unique_accs = sorted({a["accession"] for a in amplicons})
    if len(unique_accs) == 1:
        sizes = [a["product_size"] for a in amplicons]
        st.markdown("### الجل الافتراضي: حارة واحدة (نفس الكائن/العزلة)")
        img = draw_single_lane_gel(sorted(sizes), ladder="auto")
        st.image(img, use_column_width=True)
    else:
        # 2) عدة accessions ⇒ عدة حارات
        by_acc = defaultdict(list)
        for a in amplicons:
            by_acc[a["accession"]].append(a["product_size"])
        lane_labels = []
        lanes = []
        # نحد العدد الظاهر لسهولة القراءة (يمكن تعديل/إزالة القيد)
        for acc, lst in list(by_acc.items())[:8]:
            lane_labels.append(acc)
            lanes.append(sorted(lst))
        st.markdown("### الجل الافتراضي: عدة حارات (حارة لكل accession)")
        img = draw_multi_lane_gel(lanes, lane_labels, ladder="auto")
        st.image(img, use_column_width=True)

    # تنزيل الـ FASTA لكل الأمبليكونات
    st.markdown("### تنزيل التسلسلات (FASTA)")
    fasta_all = []
    for a in amplicons:
        start = min(a["fwd"]["start"], a["rev"]["start"])
        end   = max(a["fwd"]["end"],   a["rev"]["end"])
        try:
            fasta = efetch_sequence(a["accession"], start, end, rettype="fasta")
            fasta_all.append(fasta)
        except Exception as e:
            st.warning(f"لم يمكن جلب التسلسل من EFetch: {e}")
    if fasta_all:
        payload = "".join(fasta_all).encode()
        st.download_button("تنزيل جميع الأمبليكونات (FASTA)", data=payload,
                           file_name="amplicons.fasta", mime="text/plain")

# ============ تذييل ============
st.markdown("""
<div class="footer">
  تم إنشاء هذا الموقع بواسطة <b>Mahmood Al-Mualm</b> ( <b>محمود أحمد محي المعلّم</b> ).<br/>
  البريد الإلكتروني: <a href="mailto:mahmoodalmoalm@gmail.com">mahmoodalmoalm@gmail.com</a> — 
  الهاتف/واتساب: <a href="tel:+9647730585329">+964 7730585329</a>
</div>
""", unsafe_allow_html=True)

