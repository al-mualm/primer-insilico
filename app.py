import io
import time
import math
import textwrap
import requests
import xml.etree.ElementTree as ET

import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

# ------------------------ إعدادات عامة ------------------------
st.set_page_config(page_title="محاكاة PCR بالعربية", page_icon="🧬", layout="centered")

st.markdown("## 🧬 محاكاة PCR بالعربية (NCBI)")
st.markdown("أدخل البرايمرين (Forward/Reverse)، وسيقوم التطبيق بالبحث في قواعد بيانات NCBI (BLAST) وربط الضربات لاستخراج الأمبليكون، ثم عرض جل إلكترُوفوريز افتراضي.")

# ------------------------ وظائف مساعدة ------------------------
NCBI_TOOL = "primer-insilico-ar"
NCBI_EMAIL = "mahmoodalmoalm@gmail.com"  # غيّره إلى بريدك الإلكتروني للتعرّف لدى NCBI

def to_fasta(seq, title="query"):
    seq = "".join([c for c in seq.strip().upper() if c in "ACGTURYKMSWBDHVN"])
    wrapped = "\n".join(textwrap.wrap(seq, width=70))
    return f">{title}\n{wrapped}\n"

def blast_put_and_wait(fasta_query, database="refseq_genomic", task="blastn-short", program="blastn", sleep_s=3.0, max_wait_s=120):
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
        raise RuntimeError("تعذّر الحصول على RID من استجابة BLAST.")
    start_t = time.time()
    while True:
        get_params = {"CMD":"Get","RID":rid,"FORMAT_TYPE":"XML"}
        g = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi", params=get_params, timeout=60)
        g.raise_for_status()
        txt = g.text
        if "Status=WAITING" in txt:
            if time.time()-start_t > max_wait_s:
                raise TimeoutError("انتهى وقت الانتظار لنتيجة BLAST.")
            time.sleep(sleep_s)
            continue
        if "Status=FAILED" in txt:
            raise RuntimeError("فشل تشغيل BLAST على خادم NCBI.")
        if "Status=UNKNOWN" in txt:
            raise RuntimeError("غير معروف RID لدى BLAST.")
        if "QBlastInfoBegin" in txt and "QBlastInfoEnd" in txt and "ThereAreHits=yes" not in txt and "ThereAreHits=Yes" not in txt and "<Hit>" not in txt:
            return None
        return txt

def parse_blast_hits(blast_xml_text):
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
        acc = hit.findtext("Hit_accession")
        for hsp in hit.iterfind("Hit_hsps/Hsp"):
            hit_from = int(hsp.findtext("Hsp_hit-from"))
            hit_to = int(hsp.findtext("Hsp_hit-to"))
            strand = "+" if hit_from < hit_to else "-"
            evalue = float(hsp.findtext("Hsp_evalue"))
            identity = int(hsp.findtext("Hsp_identity"))
            align_len = int(hsp.findtext("Hsp_align-len"))
            pct_id = 100.0 * identity / max(1, align_len)
            hits.append((acc, hit_from, hit_to, strand, evalue, pct_id))
    return hits

def pair_hits(forward_hits, reverse_hits, product_min=80, product_max=2000):
    amplicons = []
    by_acc_rev = {}
    for acc, hfrom, hto, strand, evalue, pct in reverse_hits:
        by_acc_rev.setdefault(acc, []).append((hfrom, hto, strand, evalue, pct))
    for acc, ffrom, fto, fstrand, fev, fpct in forward_hits:
        if acc not in by_acc_rev:
            continue
        for rfrom, rto, rstrand, rev, rpct in by_acc_rev[acc]:
            if fstrand == rstrand:
                continue
            f_end = max(ffrom, fto)
            r_end = min(rfrom, rto)
            size = abs(r_end - f_end) + 1
            if product_min <= size <= product_max:
                amplicons.append({
                    "accession": acc,
                    "product_size": int(size),
                    "fwd": {"start": min(ffrom, fto), "end": max(ffrom, fto), "strand": fstrand, "evalue": fev, "pct_identity": fpct},
                    "rev": {"start": min(rfrom, rto), "end": max(rfrom, rto), "strand": rstrand, "evalue": rev, "pct_identity": rpct},
                })
    amplicons.sort(key=lambda x: x["product_size"])
    return amplicons

def efetch_sequence(accession, start, end, rettype="fasta"):
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": rettype,
        "retmode": "text",
        "seq_start": start,
        "seq_stop": end,
        "strand": 1,
        "tool": NCBI_TOOL,
        "email": NCBI_EMAIL,
    }
    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params, timeout=60)
    r.raise_for_status()
    return r.text

def log_migration(bp, a=100.0, b=50.0):
    return a - b * math.log10(bp)

def draw_in_silico_gel(product_sizes, ladder="auto"):
    ladder_100 = list(range(100, 1600, 100))
    ladder_1k = [250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000]
    if ladder == "auto":
        max_size = max(product_sizes) if product_sizes else 0
        ladder = "100bp" if max_size <= 800 else "1kb"
    ladder_sizes = ladder_100 if ladder == "100bp" else ladder_1k

    def y_of(sizes):
        return [log_migration(s) for s in sizes]

    y_all = y_of(ladder_sizes + product_sizes)
    y_min, y_max = min(y_all) - 5, max(y_all) + 5

    fig, ax = plt.subplots(figsize=(4.5, 6), dpi=160)
    lane_width = 1.0
    lane_gap = 2.0
    lane0_x = 0.0
    lane1_x = lane0_x + lane_width + lane_gap

    for yy, sz in zip(y_of(ladder_sizes), ladder_sizes):
        ax.hlines(yy, lane0_x, lane0_x + lane_width, linewidth=3)
        ax.text(lane0_x - 0.25, yy, f"{sz}", va="center", ha="right", fontsize=8)

    for yy, sz in zip(y_of(product_sizes), product_sizes):
        ax.hlines(yy, lane1_x, lane1_x + lane_width, linewidth=4)
        ax.text(lane1_x + lane_

