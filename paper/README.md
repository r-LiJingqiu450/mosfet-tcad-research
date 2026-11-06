# Manuscript

`main.tex` is the manuscript; `main.pdf` is the built version.

Every figure, table and quoted result is generated from the released benchmark files:

```bash
python paper/scripts/make_assets.py      # writes paper/generated/ (numbers.tex, tab_*.tex, inputs.sha256) and paper/figures/
cd paper && latexmk -pdf main.tex        # pdflatex + bibtex; needs siunitx, booktabs, natbib, cleveref
```

`paper/generated/inputs.sha256` lists the SHA-256 of every input file the script read. The references in
`references.bib` record, in a non-printing `checked` field, where each entry was verified.

Before submission, the author fills in: the author line, the "Use of AI" confirmation, the acknowledgments, and the
repository URL and DOI in "Data and code availability".
