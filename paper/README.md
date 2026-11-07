# Manuscript

`main.tex` is the manuscript; `main.pdf` is the built version.

Every figure, table and quoted result is generated from the released benchmark files:

```bash
python paper/scripts/make_assets.py      # writes paper/generated/ (numbers.tex, tab_*.tex, inputs.sha256) and paper/figures/
cd paper && latexmk -pdf main.tex        # pdflatex + bibtex; needs siunitx, booktabs, natbib, cleveref
```

`paper/generated/inputs.sha256` lists the SHA-256 of every input file the script read. The references in
`references.bib` record, in a non-printing `checked` field, where each entry was verified.

The AI-use disclosure approved by the author (2026-11-07) is the section "Use of generative AI"; a journal version
places it where the journal asks (the Methods section for the Journal of Computational Electronics, the Acknowledgment
for IEEE journals). "Data and code availability" points to release v1.0.0 and its Zenodo archive, doi:10.5281/zenodo.24040430.
