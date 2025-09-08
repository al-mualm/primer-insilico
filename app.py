# -*- coding: utf-8 -*-
# FINAL VERSION v3.1 - Corrected Syntax

import io
import re
import time
import math
from typing import List, Dict

import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
import requests
from bs4 import BeautifulSoup
import firebase_admin
from firebase_admin import credentials, firestore

# ---------------- Page & styles ----------------
st.set_page_config(page_title="محاكاة PCR", layout="wide")
st.markdown("""
<style>
html, body, [class*="css"] { direction: rtl; font-family: "Noto Naskh Arabic","Tajawal","Cairo",sans-serif; }
.header { background: linear-gradient(135deg,#ecf3ff 0%,#f0fff6 100%); border-radius:18px; padding:16px 20px; border:1px solid #e6eef7; margin-bottom:14px;}
.header h1 { margin:0 0 6px 0; font-size:1.6rem; }
.intro-box { background-color: #f8f9fa; border-radius: 10px; padding: 15px; margin-bottom: 15px; border: 1px solid #dee2e6; }
.small { font-size:.9rem; color:#475569;}
.badges span { display:inline-block; margin:4px 6px 0 0; padding:6px 10px; border-radius:999px; border:1px solid #e3e8ef; background:#fff; font-size:.9rem; }
code { direction: ltr; text-align: left; display: block; white-space: pre; }
</style>
""", unsafe_allow_html=True)

# --- Firebase Firestore Initialization for Visitor Counter ---
def init_firebase():
    """Initializes Firebase app if not already initialized."""
    try:
        firebase_creds = dict(st.secrets["firebase"])
        if not all(k in firebase_creds for k in ["type", "project_id", "private_key_id", "private_key", "client_email", "client_id"]):
             st.error("Firebase credentials are not correctly configured in secrets.")
             return None
        cred = credentials.Certificate(firebase_creds)
        if not firebase_admin._apps:
            firebase_admin.initialize_app(cred)
        return firestore.client()
    except Exception as e:
        st.error(f"Failed to initialize Firebase: {e}. Ensure your secrets are set correctly.")
        return None

def get_and_increment_visitor_count(db):
    """Gets the current count, increments it, and returns the new count."""
    if db is None: return "N/A"
    doc_ref = db.collection('app_stats').document('visitors')
    doc = doc_ref.get()
    if doc.exists:
        new_count = doc.to_dict().get('count', 0) + 1
        doc_ref.set({'count': new_count})
        return new_count
    else:
        doc_ref.set({'count': 1})
        return 1

if 'visitor_count' not in st.session_state:
    db = init_firebase()
    st.session_state.visitor_count = get_and_increment_visitor_count(db) if db else "Error"

# --- Main Header ---
st.markdown("""
<div class="header">
  <h1>محاكاة PCR — بحث سريع عبر الويب + جيل افتراضي</h1>
  <div class="small">
    صُمِّم بواسطة <b>Mahmood Al-Mualm</b> — <b>محمود أحمد محي المعلّم</b> ·
    البريد: <a href="mailto:mahmoodalmoalm@gmail.com">mahmoodalmoalm@gmail.com</a>
  </div>
  <div class="badges"><span>Web In-silico PCR</span><span>In-silico Gel</span><span>واجهة عربية</span></div>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="intro-box">
<p>عزيزي/عزيزتي الباحث(ة)، هذه الأداة توفر طريقة بسيطة باللغة العربية لمحاكاة تضخيم البادئات على مراجع جينومية معروفة عبر الويب.</p>
<p><b>يرجى الانتباه:</b> هذه محاكاة حاسوبية فقط وهي خطوة تمهيدية مساعدة، وليست بديلاً عن التحقق المخبري.</p>
</div>
""", unsafe_allow_html=True)

# ---------------- UCSC config & Session ----------------
UCSC_JSON = "https://api.genome.ucsc.edu/hgPcr"
UCSC_HTML = ["https://genome.ucsc.edu/cgi-bin/hgPcr", "https://genome-euro.ucsc.edu/cgi-bin/hgPcr", "https://genome-asia.ucsc.edu/cgi-bin/hgPcr"]
UCSC_GENOMES = {"Human (hg38)": ("Human", "hg38"), "Human (hg19)": ("Human", "hg19"), "Mouse (mm39)": ("Mouse", "mm39"), "Rat (rn7)": ("Rat", "rn7")}
API_KEY = st.secrets.get("scrapingbee", {}).get("api_key")
http_session = requests.Session()

# ---------------- Helper Functions ----------------
def _clean(seq: str) -> str: return re.sub(r"[^ACGTNacgtn]", "", (seq or "")).upper()
def _looks_like_cloudflare(text: str) -> bool: return "cloudflare" in text.lower() or "challenge" in text.lower()

@st.cache_data(ttl=24*60*60, show_spinner=False)
def _http_get_cached(url: str, params: Dict, timeout: int) -> requests.Response:
    if not API_KEY: st.error("ScrapingBee API Key not found in secrets."); st.stop()
    proxy_url = 'https://app.scrapingbee.com/api/v1/'
    target_url = requests.Request('GET', url, params=params).prepare().url
    proxy_params = {'api_key': API_KEY, 'url': target_url}
    r = http_session.get(proxy_url, params=proxy_params, timeout=timeout)
    r.raise_for_status()
    return r

def calculate_primer_properties(primer_sequence: str) -> Dict:
    if not primer_sequence: return {"Length": 0, "GC%": 0.0, "Tm (°C)": 0.0}
    return {"Length": len(primer_sequence), "GC%": round(gc_fraction(primer_sequence) * 100, 2), "Tm (°C)": round(mt.Tm_NN(primer_sequence), 2)}

def run_pcr_search(fwd: str, rev: str, org: str, db: str, max_bp: int) -> List[Dict]:
    all_hits = []
    try: all_hits.extend(ucsc_via_json(fwd, rev, org, db, max_bp) or [])
    except Exception: pass
    if not all_hits:
        try: all_hits.extend(ucsc_via_html(fwd, rev, org, db, max_bp) or [])
        except Exception as e: st.error(f"Failed to fetch data from all sources: {e}")
    return all_hits

def ucsc_via_json(fwd, rev, org, db, max_bp):
    params = dict(org=org, db=db, wp_f=fwd, wp_r=rev, wp_size=int(max_bp))
    r = _http_get_cached(UCSC_JSON, params, 45)
    data = r.json()
    products = []
    for key in ("results","pcr","items","products"):
        if key in data and isinstance(data[key], list):
            for item in data[key]:
                if prod := _coerce_item(item): products.append(prod)
    return products

def ucsc_via_html(fwd, rev, org, db, max_bp):
    last_err = None
    for url in UCSC_HTML:
        try:
            params = {"org": org, "db": db, "wp_f": fwd, "wp_r": rev, "wp_size": int(max_bp), "Submit": "submit"}
            r = _http_get_cached(url, params, 45)
            txt = r.text
            if _looks_like_cloudflare(txt):
                last_err = RuntimeError("Request blocked by Cloudflare."); time.sleep(1); continue
            soup = BeautifulSoup(txt, "lxml")
            pre_txt = "\n".join(p.get_text("\n") for p in soup.find_all("pre")) or txt
            if products := (parse_fasta_products(pre_txt) or parse_html_products(pre_txt)): return products
        except Exception as e: last_err = e; time.sleep(1)
    if last_err: raise last_err
    return []

def _coerce_item(item):
    if isinstance(item, list) and item and isinstance(item[0], dict): item = item[0]
    if not isinstance(item, dict): return {}
    chrom, start, end = item.get("chrom") or item.get("chr"), item.get("start"), item.get("end")
    if seq := item.get("sequence", "") and not (chrom and start and end):
        if m := re.search(r">(chr[\w\.\-]+):(\d+)-(\d+)", seq.splitlines()[0]):
            chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if chrom and start and end:
        size = abs(int(end) - int(start)) + 1
        return {"chrom": chrom, "start": int(start), "end": int(end), "size": size, "strand": item.get("strand", "+"), "sequence": item.get("sequence", "").replace("\r","")}
    return {}

def parse_fasta_products(text: str) -> List[Dict]:
    products = []
    for blk in re.split(r"(?m)^>", text):
        if not (blk := blk.strip()): continue
        header, *seq_lines = blk.splitlines()
        if not (m := re.search(r"(chr[\w\.\-]+):(\d+)-(\d+)(?:\(([-+])\))?", header)): continue
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        strand = m.group(4) if m.lastindex and m.group(m.lastindex) in ["+","-"] else "+"
        seq = "".join(s.strip() for s in seq_lines if re.fullmatch(r"[ACGTNacgtn]+", s.strip()))
        products.append({"chrom": chrom, "start": start, "end": end, "size": abs(end-start)+1, "strand": strand, "sequence": seq.upper()})
    return products

def parse_html_products(text: str) -> List[Dict]:
    products = []
    for line in text.splitlines():
        if not (m := re.search(r"(chr[\w\.\-]+):(\d+)-(\d+).*?([+-])?", line)): continue
        chrom, start, end, strand = m.group(1), int(m.group(2)), int(m.group(3)), (m.group(4) if m.group(4) in ["+","-"] else "+")
        products.append({"chrom": chrom, "start": start, "end": end, "size": abs(end-start)+1, "strand": strand, "sequence": ""})
    return products

def _gel_y(bp, a=100.0, b=50.0): return a - b * math.log10(max(bp, 1))
def render_gel(sizes: List[int]) -> io.BytesIO:
    lad100, lad1k = list(range(100, 1600, 100)), [250,500,750,1000,1500,2000,3000,4000,5000,6000,8000,10000]
    ladder = lad100 if (max(sizes) if sizes else 0) <= 1500 else lad1k
    all_y = [_gel_y(x) for x in ladder + sizes]
    ymin, ymax = min(all_y)-5, max(all_y)+5
    fig, ax = plt.subplots(figsize=(5.2, 6), dpi=170)
    x_lad, x_s, lane_w = 0.0, 3.0, 1.0
    for yy, sz in zip([_gel_y(x) for x in ladder], ladder):
        ax.hlines(yy, x_lad, x_lad+lane_w, linewidth=3); ax.text(x_lad-0.25, yy, f"{sz}", va="center", ha="right", fontsize=8)
    for yy, sz in zip([_gel_y(x) for x in sizes], sizes):
        ax.hlines(yy, x_s, x_s+lane_w, linewidth=4); ax.text(x_s+lane_w+0.2, yy, f"{sz} bp", va="center", fontsize=8)
    ax.text(x_lad+lane_w/2, ymin+2, "Ladder", ha="center", fontsize=9); ax.text(x_s+lane_w/2, ymin+2, "Sample", ha="center", fontsize=9)
    ax.set_xlim(-1.2, x_s+lane_w+1.2); ax.set_ylim(ymin, ymax); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); buf = io.BytesIO(); fig.savefig(buf, format="png", bbox_inches="tight"); plt.close(fig); buf.seek(0)
    return buf

# --- UI Starts Here ---
with st.sidebar:
    st.subheader("إعدادات المرجع")
    genome_label = st.selectbox("المرجع", list(UCSC_GENOMES.keys()), index=0)
    org, db = UCSC_GENOMES[genome_label]
    max_bp  = st.number_input("أكبر حجم (bp)", 50, 10000, 4000)
    st.markdown("---")
    st.write(f"**عدد الزوار:** {st.session_state.get('visitor_count', 'Loading...')}")

st.markdown("## (Primers) إدخال البادئات")
c1, c2 = st.columns(2)
fwd_in, rev_in = c1.text_input("Forward primer", "CATACCACAATTGCAT"), c2.text_input("Reverse primer", "AAGAAGAAGAGAGGGGG")
fwd, rev = _clean(fwd_in), _clean(rev_in)

if fwd and rev:
    st.markdown("#### خصائص البادئات")
    df = pd.DataFrame([calculate_primer_properties(fwd), calculate_primer_properties(rev)], index=["Forward", "Reverse"])
    st.table(df)

if st.button("تشغيل المحاكاة"):
    if not fwd or not rev: st.error("يرجى إدخال كلا البادئين."); st.stop()
    with st.spinner("...جارٍ البحث عبر الويب"):
        all_hits = run_pcr_search(fwd, rev, org, db, max_bp)
    if not all_hits: st.error("لم يتم العثور على نواتج."); st.stop()
    
    all_hits = sorted(all_hits, key=lambda h: h.get("size", 10**9))[:5]
    
    st.markdown("---"); st.markdown("## النتائج (أفضل ٥)")
    for i, h in enumerate(all_hits, 1):
        product_header = f"{h.get('chrom','?')}:{h.get('start','?'):,}-{h.get('end','?'):,}"
        st.markdown(f"#### 🎯 المنتج {i}: `{product_header}`")
        st.write(f"**طول المنتج:** {h.get('size', '?')} bp")
    
    if sizes := [h["size"] for h in all_hits if "size" in h]:
        st.markdown("### الجل الافتراضي"); st.image(render_gel(sorted(sizes)), use_column_width=True)
