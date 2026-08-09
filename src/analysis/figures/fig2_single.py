#!/usr/bin/env python3
"""
Figure 2 - one detailed worked example, as NATIVE PowerPoint shapes.

CD4+ T cell branch, before (v2024-08-13) vs after (v2026-06-23). Every node is a
real rounded rectangle with editable text and every is_a edge is a real connector,
so the figure can be edited directly in PowerPoint.

Visual grammar
    box fill/border = what changed about the TERM
        blue fill     new term
        teal border   relabelled
        grey          unchanged
        dashed        context (parent shown for orientation)
    arrow colour = what changed about the EDGE
        amber         this is_a edge was re-parented

Structure verified against the real releases:
    OLD  Abnormal T cell morphology > Abnormal T cell count
         > Abnormal T cell subset distribution
         > Abnormal proportion of CD4-positive T cells
         > Decreased / Increased proportion of CD4-positive T cells
         (no absolute-count term existed)
    NEW  Abnormal T cell morphology > Abnormal T cell subset number (HP:0025540)
         > Abnormal total CD4+ T cell number (HP:0025183)
         > Decreased (HP:5210415) / Increased (HP:5210416) total CD4+ T cell number
         > count (HP:5210418 / HP:5210417) + proportion (HP:0032218 / HP:0032219)

Deps:  pip install python-pptx lxml
Run:   python fig2_single.py
"""

import pathlib
HERE = pathlib.Path(__file__).resolve().parent

from lxml import etree
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

TPL = str(HERE / "esid4hpo_figures.pptx")   # deck whose theme/size to reuse
OUT = str(HERE / "_work" / "figures" / "ESID4HPO_Figure2_single.pptx")

NAVY=RGBColor(0x00,0x37,0x54); BLUE=RGBColor(0x4E,0x7E,0x96)
GREY=RGBColor(0xE7,0xEA,0xEC); AMBER=RGBColor(0xC8,0xA8,0x70)
TEAL=RGBColor(0x6C,0xBC,0xC5); INK=RGBColor(0x2B,0x2F,0x33)
MUT=RGBColor(0x5F,0x70,0x78); BORD=RGBColor(0xC4,0xC9,0xCE)
WHITE=RGBColor(0xFF,0xFF,0xFF); ARROW=RGBColor(0x9A,0xA3,0xAD)
TINT=RGBColor(0xF2,0xF5,0xF7)

STYLE={ "ctx":(WHITE,BORD,1.0,MUT,False,True),
        "exist":(GREY,BORD,1.0,INK,False,False),
        "new":(BLUE,BLUE,1.0,WHITE,True,False),
        "rep":(GREY,AMBER,2.25,INK,False,False),
        "relab":(GREY,TEAL,2.25,INK,False,False)}


def drop_all(prs):
    lst=prs.slides._sldIdLst
    for s in list(lst):
        rId=s.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
        prs.part.drop_rel(rId); lst.remove(s)


def white_bg(sl):
    """BIH master maps bg1 -> dk1 (navy); figure panels need white."""
    cSld=sl._element.find(qn("p:cSld"))
    for o in cSld.findall(qn("p:bg")): cSld.remove(o)
    bg=etree.SubElement(cSld,qn("p:bg")); bgPr=etree.SubElement(bg,qn("p:bgPr"))
    sf=etree.SubElement(bgPr,qn("a:solidFill")); c=etree.SubElement(sf,qn("a:srgbClr"))
    c.set("val","FFFFFF"); etree.SubElement(bgPr,qn("a:effectLst")); cSld.insert(0,bg)


def tb(s,x,y,w,h,txt,size=9,bold=False,color=INK,align=PP_ALIGN.LEFT,italic=False):
    t=s.shapes.add_textbox(Inches(x),Inches(y),Inches(w),Inches(h))
    f=t.text_frame; f.word_wrap=True
    f.margin_left=f.margin_right=f.margin_top=f.margin_bottom=0
    p=f.paragraphs[0]; p.text=txt; p.alignment=align
    for r in p.runs:
        r.font.size=Pt(size); r.font.bold=bold; r.font.color.rgb=color; r.font.italic=italic
    return t


def box(s,x,y,w,h,label,kind,size):
    fill,line,lw,tc,bold,dash=STYLE[kind]
    sh=s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,Inches(x),Inches(y),Inches(w),Inches(h))
    sh.fill.solid(); sh.fill.fore_color.rgb=fill
    sh.line.color.rgb=line; sh.line.width=Pt(lw); sh.shadow.inherit=False
    if dash:
        ln=sh.line._get_or_add_ln()
        d=etree.SubElement(ln,qn("a:prstDash")); d.set("val","dash")
    try: sh.adjustments[0]=0.12
    except Exception: pass
    f=sh.text_frame; f.word_wrap=True
    f.margin_left=f.margin_right=Inches(0.03); f.margin_top=f.margin_bottom=Inches(0.01)
    f.vertical_anchor=MSO_ANCHOR.MIDDLE
    p=f.paragraphs[0]; p.text=label; p.alignment=PP_ALIGN.CENTER
    for r in p.runs:
        r.font.size=Pt(size); r.font.bold=bold; r.font.color.rgb=tc
    return sh


def arrow(s,x1,y1,x2,y2,color=ARROW,width=1.0):
    cn=s.shapes.add_connector(MSO_CONNECTOR.STRAIGHT,Inches(x1),Inches(y1),Inches(x2),Inches(y2))
    cn.line.color.rgb=color; cn.line.width=Pt(width)
    ln=cn.line._get_or_add_ln()
    te=etree.SubElement(ln,qn("a:tailEnd")); te.set("type","triangle"); te.set("w","med"); te.set("len","med")
    return cn


prs=Presentation(TPL); drop_all(prs)
layout=next((l for l in prs.slide_layouts if l.name=="Bilder zweispaltig"), prs.slide_layouts[0])
s=prs.slides.add_slide(layout); white_bg(s)
for ph in list(s.placeholders): ph._element.getparent().remove(ph._element)

frame=s.shapes.add_shape(MSO_SHAPE.RECTANGLE,Inches(0.16),Inches(0.10),Inches(9.68),Inches(7.28))
frame.fill.background(); frame.line.color.rgb=BORD; frame.line.width=Pt(0.75); frame.shadow.inherit=False

tb(s,0.30,0.20,9.4,0.28,
   "Restructuring the CD4+ T cell branch: grouping by biology, not by unit of measurement",
   size=12.5,bold=True,color=NAVY,align=PP_ALIGN.CENTER)

# ---------------- BEFORE (left) ----------------
BX, BW = 0.34, 3.30
tb(s,BX,0.60,BW,0.24,"Before   (v2024-08-13)",size=10.5,bold=True,color=MUT,align=PP_ALIGN.CENTER)
H=0.40; y=0.96; VG=0.30
rows_b=[
 [("Abnormal T cell morphology","ctx")],
 [("Abnormal T cell count","ctx")],
 [("Abnormal T cell subset distribution","exist")],
 [("Abnormal proportion of CD4-positive T cells","exist")],
 [("Decreased proportion of CD4-positive T cells","exist"),
  ("Increased proportion of CD4-positive T cells","exist")],
]
posb={}
for i,row in enumerate(rows_b):
    n=len(row); gap=0.10
    bw=(BW-gap*(n-1))/n
    size=7.5 if n==1 else 6.5
    yy=y+i*(H+VG)
    for j,(lab,kind) in enumerate(row):
        xx=BX+j*(bw+gap)
        box(s,xx,yy,bw,H,lab,kind,size)
        posb[(i,j)]=(xx+bw/2,yy,yy+H)
for (ci,cj),(pi,pj) in [((1,0),(0,0)),((2,0),(1,0)),((3,0),(2,0)),((4,0),(3,0)),((4,1),(3,0))]:
    cx,ct,_=posb[(ci,cj)]; px,_,pb=posb[(pi,pj)]
    arrow(s,cx,ct,px,pb)

# note: no absolute-count term existed
nb=s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,Inches(BX),Inches(y+5*(H+VG)-0.05),Inches(BW),Inches(0.46))
nb.fill.solid(); nb.fill.fore_color.rgb=TINT
nb.line.color.rgb=BORD; nb.line.width=Pt(1.0); nb.shadow.inherit=False
ln=nb.line._get_or_add_ln(); d=etree.SubElement(ln,qn("a:prstDash")); d.set("val","dash")
f=nb.text_frame; f.word_wrap=True; f.vertical_anchor=MSO_ANCHOR.MIDDLE
p=f.paragraphs[0]; p.text="No term for an absolute CD4+ count:\nonly proportions could be recorded"
p.alignment=PP_ALIGN.CENTER
for r in p.runs: r.font.size=Pt(7.5); r.font.color.rgb=MUT; r.font.italic=True

# ---------------- AFTER (right) ----------------
AX, AW = 4.10, 5.60
tb(s,AX,0.60,AW,0.24,"After   (v2026-06-23)",size=10.5,bold=True,color=MUT,align=PP_ALIGN.CENTER)
rows_a=[
 [("Abnormal T cell morphology","ctx")],
 [("Abnormal T cell subset number","relab")],
 [("Abnormal total CD4+ T cell number","new")],
]
posa={}
for i,row in enumerate(rows_a):
    yy=y+i*(H+VG)
    box(s,AX,yy,AW,H,row[0][0],row[0][1],8.0)
    posa[(i,0)]=(AX+AW/2,yy,yy+H)

# level 3: decreased | increased  (wide separation, per reviewer comment)
GRP_GAP=0.62
gw=(AW-GRP_GAP)/2
y3=y+3*(H+VG)
box(s,AX,y3,gw,H,"Decreased total CD4+ T cell number","new",7.5)
posa[(3,0)]=(AX+gw/2,y3,y3+H)
box(s,AX+gw+GRP_GAP,y3,gw,H,"Increased total CD4+ T cell number","new",7.5)
posa[(3,1)]=(AX+gw+GRP_GAP+gw/2,y3,y3+H)

# level 4: count | proportion within each group, on a tinted container that makes
# the "same biological finding" grouping visible
y4=y+4*(H+VG); H4=0.50
inner=0.08
lw_=(gw-inner)/2
for gi,gx in [(0,AX),(1,AX+gw+GRP_GAP)]:
    cont=s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,Inches(gx-0.05),Inches(y4-0.07),
                            Inches(gw+0.10),Inches(H4+0.14))
    cont.fill.solid(); cont.fill.fore_color.rgb=TINT
    cont.line.color.rgb=BORD; cont.line.width=Pt(0.75); cont.shadow.inherit=False
    s.shapes._spTree.remove(cont._element); s.shapes._spTree.insert(2,cont._element)  # send to back
    lab_c = "Decreased total CD4+ T cell count" if gi==0 else "Increased total CD4+ T cell count"
    lab_p = "Decreased total CD4+ T cell proportion" if gi==0 else "Increased total CD4+ T cell proportion"
    box(s,gx,y4,lw_,H4,lab_c,"new",6.5)
    posa[(4,gi*2)]=(gx+lw_/2,y4,y4+H4)
    box(s,gx+lw_+inner,y4,lw_,H4,lab_p,"relab",6.5)
    posa[(4,gi*2+1)]=(gx+lw_+inner+lw_/2,y4,y4+H4)

for (ci,cj),(pi,pj),st in [((1,0),(0,0),"amber"),((2,0),(1,0),None),
                           ((3,0),(2,0),None),((3,1),(2,0),None),
                           ((4,0),(3,0),None),((4,1),(3,0),"amber"),
                           ((4,2),(3,1),None),((4,3),(3,1),"amber")]:
    cx,ct,_=posa[(ci,cj)]; px,_,pb=posa[(pi,pj)]
    arrow(s,cx,ct,px,pb, AMBER if st=="amber" else ARROW, 1.75 if st else 1.0)

tb(s,AX,y4+H4+0.10,gw,0.34,
   "same biological finding\n(low CD4+ T cells)",size=7.5,color=MUT,align=PP_ALIGN.CENTER,italic=True)
tb(s,AX+gw+GRP_GAP,y4+H4+0.10,gw,0.34,
   "same biological finding\n(high CD4+ T cells)",size=7.5,color=MUT,align=PP_ALIGN.CENTER,italic=True)

# ---------------- legend ----------------
ly=5.28; lx=0.34
for fill,line,txt in [(BLUE,BLUE,"new term"),(GREY,AMBER,"re-parented"),
                      (GREY,TEAL,"relabelled"),(GREY,BORD,"unchanged"),(WHITE,BORD,"context")]:
    sq=s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,Inches(lx),Inches(ly),Inches(0.18),Inches(0.18))
    sq.fill.solid(); sq.fill.fore_color.rgb=fill
    sq.line.color.rgb=line; sq.line.width=Pt(1.5); sq.shadow.inherit=False
    tb(s,lx+0.23,ly+0.005,1.5,0.2,txt,size=8.5,color=INK)
    lx+=0.23+0.068*len(txt)+0.20

# ---------------- caption ----------------
cap=("Before curation, CD4+ T cell abnormalities could only be expressed as proportions; no term existed for an absolute "
     "count, and the hierarchy grouped terms by unit of measurement (all proportions together), splitting by direction only "
     "at the leaves. After curation, an \u201cabnormal number\u201d parent was introduced for each direction of change, so that the "
     "absolute count and the proportion of the same biological finding \u2014 e.g. a low CD4+ T cell count and a low CD4+ T cell "
     "proportion \u2014 are now sibling terms under one parent. This groups biologically equivalent observations rather than "
     "measurement types, so a query for \u201cdecreased CD4+ T cells\u201d retrieves both, while the count/proportion distinction is "
     "preserved at the leaves. The example shows all three types of curation: new terms, re-parenting of existing terms, and "
     "relabelling from laboratory nomenclature to clinical wording.")
tb(s,0.34,5.62,9.10,1.30,cap,size=8.5,color=INK)

pathlib.Path(OUT).parent.mkdir(parents=True, exist_ok=True)
prs.save(OUT)
print("saved", OUT)
