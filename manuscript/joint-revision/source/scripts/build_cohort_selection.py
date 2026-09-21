"""Build the aggregate UCMM selection chart directly from the cleaning code, without intermediate files."""
from pathlib import Path
import csv, re, html, io, subprocess
from reportlab.graphics.shapes import Drawing, Rect, String, Line, Polygon
from reportlab.graphics import renderPDF, renderSVG
from reportlab.lib.colors import HexColor
root=Path(__file__).resolve().parents[3]
app=root/'lrcbart-case-study-mm'; src=root/'manuscript/source'
# Evaluate only the in-memory preparation; skip its save section.
r_code = r"""
e <- new.env()
lines <- readLines('lrcbart-case-study-mm/private_data/data_cleaning_ucmm.R')
stop_at <- grep('^# 7[.] VALIDATE', lines)[1] - 2L
invisible(capture.output(eval(parse(text=lines[seq_len(stop_at)]), envir=e)))
write.csv(e$regimen_audit, stdout(), row.names=FALSE)
cat('---FLOW---\n')
write.csv(e$flow, stdout(), row.names=FALSE)
"""
raw = subprocess.run(['Rscript', '-e', r_code], cwd=root, check=True,
                     capture_output=True, text=True).stdout
regimen_csv, flow_csv = raw.split('---FLOW---\n')
rows = list(csv.DictReader(io.StringIO(regimen_csv)))
flow = {r['stage']: int(r['n']) for r in csv.DictReader(io.StringIO(flow_csv))}
counts={c:sum(int(r['n_c2']) for r in rows if r['classification']==c) for c in {r['classification'] for r in rows}}
assert sum(int(r['n_c3']) for r in rows)==flow['Analysis']
# Fail rather than publish stale annotated counts if the source cohort changes.
assert flow == dict(Source=797, Nonnegative_times=790, C1=303, C2=302, Analysis=200)
assert counts == dict(triplet=200, doublet=66, four_or_more=19, single_agent=4, other_or_unclear=13)
W,H=720,590; d=Drawing(W,H)
ink=HexColor('#173447'); edge=HexColor('#72858f'); blue=HexColor('#edf3f7'); green=HexColor('#e1f1ea'); cream=HexColor('#faf1e8')
def text(x,y,s,size=13,bold=False,anchor='middle'):
 d.add(String(x,H-y,s,fontName='Helvetica-Bold' if bold else 'Helvetica',fontSize=size,fillColor=ink,textAnchor=anchor))
def box(x,y,w,h,lines,fill=blue):
 d.add(Rect(x,H-y-h,w,h,rx=6,ry=6,fillColor=fill,strokeColor=edge,strokeWidth=1))
 for i,line in enumerate(lines): text(x+w/2,y+23+i*19,line,13 if i else 14,i==0)
def arrow(x1,y1,x2,y2):
 d.add(Line(x1,H-y1,x2,H-y2,strokeColor=edge,strokeWidth=1.3))
 pts=[x2,H-y2,x2-4,H-y2+7,x2+4,H-y2+7] if x1==x2 else [x2,H-y2,x2-7,H-y2-4,x2-7,H-y2+4]
 d.add(Polygon(pts,fillColor=edge,strokeColor=edge))
ys=[15,120,225,330,500]
nodes=[('UCMM source extract',flow['Source']),('After survival-time check',flow['Nonnegative_times']),('C1: induction date recorded',flow['C1']),('C2: all-regimen cohort',flow['C2']),('C3: triplet analysis cohort',flow['Analysis'])]
for i,(label,n) in enumerate(nodes):
 box(20,ys[i],320,64,[label,f'n = {n}'],green if i in (3,4) else blue)
 if i<4: arrow(180,ys[i]+64,180,ys[i+1])
ex=[['Excluded: 7','Negative OS or PFS time'],['Excluded: 487','Missing induction-start date'],['Excluded: 1','E-Rd / Elo-Rd'],['Excluded: 102','Doublets 66; single agent 4','Four or more agents 19','Other / unclear 13']]
for i,lines in enumerate(ex):
 mid=(ys[i]+64+ys[i+1])/2
 arrow(180,mid,395,mid)
 box(400,mid-(32 if i<3 else 53),300,64 if i<3 else 106,lines,cream)
text(360,580,'Triplets count treatment agents, including steroids; supportive zoledronic acid is not counted.',10)
fig=src/'figures/ucmm_cohort_selection';renderPDF.drawToFile(d,str(fig.with_suffix('.pdf')))
# Keep the HTML chart inline so the codebook stays portable.
f=app/'private_data/ucmm_elokrd_codebook.html';s=f.read_text()
svg=renderSVG.drawToString(d);svg=svg[svg.index('<svg'):];svg=svg.replace('<svg ', '<svg style="width:100%;max-width:900px;height:auto;" ',1)
chart='<div id="ucmm-cohort-selection"><div class="section-title">UCMM cohort selection</div>'+svg+'<p class="study-text">The analysis cohort includes 200 triplet-treated controls; merged with 30 EloKRd patients, n=230.</p></div>'
if 'id="ucmm-cohort-selection"' in s:
 s=re.sub(r'<div id="ucmm-cohort-selection">.*?</p></div>',chart,s,count=1,flags=re.S)
else:
 marker='<div class="section-title" style="margin-top:0;">Drug Concordance'
 s=s.replace(marker,chart+'\n'+marker,1)
# Correct cohort headers and cells in cross-dataset concordance.
start=s.index('<thead><tr><th style="width:6%;">Abbrev.');end=s.index('</table>',start)
part=s[start:end]
header='<thead><tr><th style="width:6%;">Abbrev.</th><th style="width:18%;">Full Regimen (Generic Names)</th><th style="width:11%;">EloKRD (n=30)</th><th style="width:10%;">KRd (n=25)</th><th style="width:16%;">UCMM C1: Has Tx (n=303)</th><th style="width:12%;">UCMM C2: all regimens (n=302)</th><th style="width:9%;">UCMM C3: triplets (n=200)</th><th style="width:18%;">BCCRC (n=386)</th></tr></thead>'
part=re.sub(r'<thead><tr>.*?</tr></thead>',header,part,count=1,flags=re.S)
values_c2={'E-Rd':'&mdash;','KRd':'38','VRd':'107','VCd':'22','Vd':'41','VMP':'1','Rd':'20','MP':'1','Thal':'TD 4, VTD 4, VTD-PACE 4, VTD-Cy 1','Dara-VRd':'1','Dara-KRd':'2'}
values_c3={'E-Rd':'&mdash;','KRd':'38','VRd':'107','VCd':'22','Vd':'&mdash;','VMP':'1','Rd':'&mdash;','MP':'&mdash;','Thal':'VTD 4','Dara-VRd':'&mdash;','Dara-KRd':'&mdash;'}
def fixrow(m):
 row=m.group(0); label=re.search(r'<b>(.*?)</b>',row).group(1)
 if label not in values_c2:return row
 cells=list(re.finditer(r'<td\b[^>]*>.*?</td>',row))
 replacements={5:values_c2[label],6:values_c3[label]}
 for pos in sorted(replacements,reverse=True):
  c=cells[pos]; value=replacements[pos]
  cell='<td>'+value+'</td>' if value=='&mdash;' else '<td class="avail">'+value+'</td>'
  row=row[:c.start()]+cell+row[c.end():]
 return row
part=re.sub(r'<tr><td><b>.*?</tr>',fixrow,part);s=s[:start]+part+s[end:]
note='<p class="study-text" style="margin-top:8px;font-size:0.82em;color:#666;"><b>Key overlaps:</b> VRd, VCd, Vd, and Rd appear in both UCMM and BCCRC. KRd is shared between the KRd trial and UCMM. The UCMM columns are nested: <b>C1</b> has a recorded induction date (n=303), <b>C2</b> excludes E-Rd/Elo-Rd (n=302), and <b>C3</b> retains verified triplets (n=200). Steroids count as treatment agents; supportive zoledronic acid does not. Classification uses the recorded induction components rather than the source dtq flag.<br>BCCRC (Dx 2009&ndash;2013) reflects an older treatment era: Dex-only was common (n=181). UCMM (2010+) and the two trials (2013+, 2017+) reflect more modern regimens with proteasome inhibitors and immunomodulatory drugs.</p>'
s=re.sub(r'<p class="study-text" style="margin-top:8px;font-size:0.82em;color:#666;">\s*<b>Key overlaps:</b>.*?</p>',note,s,count=1,flags=re.S)
# Replace the UCMM overview with the current cohort flow and regimen audit.
table='<div class="tab-content" id="tab-realworld"><div class="overview-wrap"><div class="section-title" style="margin-top:0;">Synthetic Control &mdash; UCMM</div><p class="study-text">The cleaning funnel retains records with nonnegative survival times, requires a recorded induction date, excludes E-Rd/Elo-Rd, and then retains verified triplets: <b>C1 n=303</b>, <b>C2 n=302</b>, and <b>C3 n=200</b>.</p><div class="section-title">Cohort 3 &mdash; verified triplets (n=200)</div><p class="study-text">Components were checked against induction_abstracted. Steroids count as treatment agents; supportive zoledronic acid does not.</p><table class="cmp-table"><thead><tr><th>Recorded regimen</th><th>C3 patients</th><th>Classification note</th></tr></thead><tbody>'
for r in rows:
 if int(r['n_c3']):table+='<tr><td>'+html.escape(r['regimen'])+'</td><td>'+r['n_c3']+'</td><td>'+html.escape(r['note'])+'</td></tr>'
table+='<tr><td><b>Total</b></td><td><b>200</b></td><td>Verified three-agent treatment regimens.</td></tr></tbody></table></div></div>\n'
a=s.index('<div class="tab-content" id="tab-realworld">')
b=s.index('<div class="tab-content" id="tab-timeline">',a)
s=s[:a]+table+s[b:]
f.write_text(s)
print('Generated chart and updated codebook; triplet analysis n=',flow['Analysis'])
