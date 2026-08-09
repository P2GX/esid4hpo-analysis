#!/usr/bin/env python3
"""
Build the Figure 2 PowerPoint (one worked example per slide) from the BIH template.

Each slide: title, before/after hierarchy diagrams, colour legend, caption.
Diagrams are produced first by fig2_hierarchy.py (Graphviz, dpi=200).

NOTE ON BACKGROUNDS
    The BIH master maps bg1 -> dk1, i.e. every layout is navy with white text.
    Figure panels need a white canvas, so each slide gets an explicit solid
    white background (see white_bg()). If you add slides by hand in PowerPoint
    they will be navy unless you set the background yourself.

Deps:  pip install python-pptx pillow
Run:   python fig2_build_pptx.py
"""

import pathlib

from PIL import Image
from lxml import etree
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

HERE = pathlib.Path(__file__).resolve().parent
TEMPLATE = HERE / "BIH_PPT_EN_stripped.pptx"      # put the BIH template here
FIGDIR = HERE / "_work" / "figures" / "fig2"      # output of fig2_hierarchy.py
OUT = HERE / "_work" / "figures" / "ESID4HPO_Figure2_panels.pptx"

DPI = 200          # must match the dpi used when rendering the Graphviz PNGs

NAVY = RGBColor(0x00, 0x37, 0x54)
BLUE = RGBColor(0x4E, 0x7E, 0x96)
GREY = RGBColor(0xE7, 0xEA, 0xEC)
AMBER = RGBColor(0xC8, 0xA8, 0x70)
TEAL = RGBColor(0x6C, 0xBC, 0xC5)
INK = RGBColor(0x2B, 0x2F, 0x33)
MUT = RGBColor(0x5F, 0x70, 0x78)
BORD = RGBColor(0xC4, 0xC9, 0xCE)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

LEGEND = [(BLUE, BLUE, "new term"),
          (GREY, AMBER, "re-parented"),
          (GREY, TEAL, "relabelled"),
          (GREY, BORD, "unchanged"),
          (WHITE, BORD, "context")]


def drop_all_slides(prs):
    """Remove the template's demo slides, including their parts."""
    lst = prs.slides._sldIdLst
    for sldId in list(lst):
        rId = sldId.get(
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
        prs.part.drop_rel(rId)
        lst.remove(sldId)


def white_bg(slide):
    """Force a solid white slide background (BIH master is navy)."""
    cSld = slide._element.find(qn("p:cSld"))
    for old in cSld.findall(qn("p:bg")):
        cSld.remove(old)
    bg = etree.SubElement(cSld, qn("p:bg"))
    bgPr = etree.SubElement(bg, qn("p:bgPr"))
    solid = etree.SubElement(bgPr, qn("a:solidFill"))
    clr = etree.SubElement(solid, qn("a:srgbClr"))
    clr.set("val", "FFFFFF")
    etree.SubElement(bgPr, qn("a:effectLst"))
    cSld.insert(0, bg)


def add_slide(prs, layout, title, before_png, after_png, caption):
    s = prs.slides.add_slide(layout)
    white_bg(s)

    # title from the layout placeholder; drop the rest so we can position freely
    for ph in list(s.placeholders):
        if ph.placeholder_format.idx == 0:
            ph.text = title
            for p in ph.text_frame.paragraphs:
                for r in p.runs:
                    r.font.size = Pt(24)
                    r.font.bold = True
                    r.font.color.rgb = NAVY
        else:
            ph._element.getparent().remove(ph._element)

    # ---- images: all arithmetic in INCHES as floats, then -> EMU ----
    TOP_IN, H_IN = 1.52, 3.90
    halves = [(0.55, 5.85), (6.95, 5.85)]        # (left_in, width_in)
    labels = ["Before   (v2024-08-13)", "After   (v2026-06-23)"]
    for (x0_in, wmax_in), png, lab in zip(halves, [before_png, after_png], labels):
        tb = s.shapes.add_textbox(Inches(x0_in), Inches(1.10),
                                  Inches(wmax_in), Inches(0.35))
        p = tb.text_frame.paragraphs[0]
        p.text = lab
        p.alignment = PP_ALIGN.CENTER
        r = p.runs[0]
        r.font.size = Pt(13)
        r.font.bold = True
        r.font.color.rgb = MUT

        iw, ih = Image.open(png).size
        img_w_in, img_h_in = iw / DPI, ih / DPI
        sc = min(wmax_in / img_w_in, H_IN / img_h_in)     # preserve aspect ratio
        w_in, h_in = img_w_in * sc, img_h_in * sc
        s.shapes.add_picture(str(png),
                             Inches(x0_in + (wmax_in - w_in) / 2.0),
                             Inches(TOP_IN + (H_IN - h_in) / 2.0),
                             width=Inches(w_in), height=Inches(h_in))

    # ---- legend ----
    ly, lx = Inches(5.62), 0.55
    for fill, line, txt in LEGEND:
        box = s.shapes.add_shape(1, Inches(lx), ly, Inches(0.22), Inches(0.22))
        box.fill.solid()
        box.fill.fore_color.rgb = fill
        box.line.color.rgb = line
        box.line.width = Pt(1.75)
        box.shadow.inherit = False
        tb = s.shapes.add_textbox(Inches(lx + 0.30), ly - Inches(0.045),
                                  Inches(2.6), Inches(0.3))
        p = tb.text_frame.paragraphs[0]
        p.text = txt
        r = p.runs[0]
        r.font.size = Pt(11)
        r.font.color.rgb = INK
        lx += 0.30 + 0.10 + 0.075 * len(txt)

    # ---- caption ----
    cb = s.shapes.add_textbox(Inches(0.55), Inches(6.05), Inches(12.2), Inches(1.1))
    cb.text_frame.word_wrap = True
    p = cb.text_frame.paragraphs[0]
    p.text = caption
    r = p.runs[0]
    r.font.size = Pt(12)
    r.font.color.rgb = INK
    return s


SLIDES = [
    ("Flow cytometry: separating absolute counts from proportions",
     "flow_before.png", "flow_after.png",
     "Before curation, CD4+ T cell subsets could only be recorded as proportions; no term existed for an absolute "
     "count, and \u201cAbnormal T cell subset distribution\u201d sat under \u201cAbnormal T cell count\u201d. After curation, an "
     "etiology-free \u201cnumber\u201d parent was introduced for each direction of change, with absolute count and proportion "
     "as distinct children (e.g. HP:5210418 Decreased total CD4+ T cell count vs. HP:0032218 Decreased total CD4+ "
     "T cell proportion). Existing proportion terms were relabelled from laboratory nomenclature to clinical wording "
     "and re-parented accordingly."),
    ("Antibodies: specific antibody response separated from immunoglobulin concentration",
     "ab_before.png", "ab_after.png",
     "Before curation, \u201cDecreased circulating level of specific antibody\u201d was classified as a subtype of decreased "
     "antibody concentration. After curation it was relabelled \u201cImpaired specific antibody response\u201d (HP:0012475) and "
     "re-parented directly under humoral immunity, since an impaired response to an antigen is not a subtype of a "
     "reduced concentration. Pneumococcal antibody deficiency was re-parented under the anti-polysaccharide branch, "
     "and new terms were added for Haemophilus influenzae type b and for graded responses to unconjugated "
     "Salmonella typhi vaccine."),
    ("Infections: cross-classification by organism, anatomical site and opportunism",
     "inf_before.png", "inf_after.png",
     "Before curation, the anatomical-site axis contained only three classes (CNS, gastrointestinal, skin) and "
     "\u201cPneumonia\u201d was not reachable from the unusual-infection branch. After curation, eight further site classes "
     "were added, and organ-specific infections are cross-classified on three axes simultaneously: HP:5210283 "
     "Cytomegalovirus pneumonitis is a child of Opportunistic viral infection, Unusual cytomegalovirus infection and "
     "Unusual lower respiratory tract infection, so the case is retrieved by a query on any axis."),
    ("Granulomatosis: an etiological axis for granulomatous disease",
     "granuloma_before.png", "granuloma_after.png",
     "Before curation, granulomatosis terms hung off unrelated morphological parents with no etiological "
     "organisation. After curation, an infectious / non-infectious \u00d7 caseating / non-caseating scaffold was created "
     "under Granulomatosis, existing terms were re-parented into it, and clinically important entities were added, "
     "including granulomatous lymphocytic interstitial lung disease (GLILD)."),
]


def main():
    if not TEMPLATE.exists():
        raise SystemExit(f"{TEMPLATE} not found - copy the BIH template next to this script.")
    prs = Presentation(str(TEMPLATE))
    drop_all_slides(prs)
    layout = prs.slide_layouts[18]      # "Bilder zweispaltig"

    for title, before, after, caption in SLIDES:
        add_slide(prs, layout, title, FIGDIR / before, FIGDIR / after, caption)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    prs.save(str(OUT))
    print(f"saved {OUT}  ({len(prs.slides._sldIdLst)} slides)")


if __name__ == "__main__":
    main()
