# Appendix visual QA

Date: 2026-09-20

## Scope

I inspected the full-size rendered PNGs `appendix-01.png` through
`appendix-27.png` and `appendix-63.png` through `appendix-68.png` with
`view_image`. The review covered equation rendering, tree-prior and Schur
displays, the fixed-region theorem pages, table labels, page flow, clipping,
and overlap. No source, PDF, or rendered image was modified.

## Result

**Pass.** I found no clipped equations, missing symbols, overlapping objects,
or table labels outside their tables on the inspected pages. The displayed
math remains legible at full size, including the conditional-update equations
on pages 14–18 and 63–65, the tree-prior formula on page 63, and the ESS/Schur
material on pages 08–09 and 14–15.

Pages 01–07 and 10–18 have clean text and equation flow. Pages 08–09 are
intentional rotated landscape table pages; both tables fit their page width and
remain fully visible. Pages 19–27 continue the theorem, benchmark, survival,
simulation-design, and Monte Carlo sections without clipping or overlap. The
page break at the continuation beginning page 20 and the continuation of the
results paragraph at page 67 are ordinary paragraph breaks and do not obscure
content.

Pages 63–65 present the joint-treatment model, restricted tree prior,
conditional treatment update, calibration description, and reproducibility
details cleanly. Tables S14 and S15 on page 66, Table S16 on page 67, and Table
S17 on page 68 are within the text block, with visible captions, headers, rules,
and bottom boundaries. The compact tables are readable at the rendered size;
none is clipped or overlapped by surrounding text.

## Boundary of review

This was a visual inspection of the requested appendix pages only. It did not
alter source files, rebuild the manuscript, or assess pages outside the listed
PNG set.
