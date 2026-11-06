# WP7 report: public release package and manuscript

**Date:** 2026-11-05 (MST) · **Performed by:** Claude, directly · **Branch:** `benchmark/wp7-release-2026-11-05`, from the WP6
report commit `d41e826`.

**Decisions by the user (2026-11-05):**
- the manuscript in LaTeX in the repository;
- MIT for code and CC BY 4.0 for data and text;
- a new public repository with Zenodo: I build a clean, benchmark-only snapshot, and the user creates the repository and
  turns on Zenodo.

**Untouched:** `phase3/`, `semiconductor_sim/`, every case file, every result of WP1–WP6, the frozen WP6 instruments.
Two files in benchmark tooling and tests changed; they are described under Deviations.

## Deliverables

| Deliverable | Where | State |
|---|---|---|
| Manuscript (13 pages; LaTeX, 4 figures, 5 tables, 25 references) | `paper/main.tex`, `paper/main.pdf` | Complete draft for the author's review. The author line, the "Use of AI" confirmation, acknowledgments, repository URL and DOI are placeholders. |
| Figure, table and number generator | `paper/scripts/make_assets.py` → `paper/generated/`, `paper/figures/` | Every quoted result (96 macros), every table and every figure is written from release files; the SHA-256 of its 132 inputs are in `paper/generated/inputs.sha256`. |
| References | `paper/references.bib` | Every entry checked on 2026-11-05 (see below); the check is recorded in a non-printing `checked` field. |
| Release files | `release/` | `README.md`, `LICENSE` (MIT), `LICENSE-DATA.md` (CC BY 4.0), `CITATION.cff`, `.zenodo.json`, `RELEASE_NOTES.md`. Author names, URL and DOI are placeholders. |
| Release builder | `benchmark/tools/build_release.py` | Builds the benchmark-only snapshot from one commit. |
| Release candidate | branch `release/v1.0.0-candidate` (commit `8cd1bed`, no parent; author `Claude <noreply@anthropic.com>`) | 1072 files, 71 MB. Built from `d41d142`. |
| Release tests | `benchmark/tests/test_release_integrity.py` | Case sources present with their SHA-256; imports resolve; inside a snapshot, `SHA256SUMS` matches every file. |

## What the snapshot holds

The snapshot keeps every file at its repository path, so that every SHA-256 a case cites and every import resolves inside it.

| Content | Files |
|---|---:|
| `benchmark/`: every tracked file, including the WP2/WP4 production records and the WP6 run records and transcripts | 931 + WP7 additions |
| Earlier project results that WP1 case files cite as sources (`results/phase3/…`) | 89 |
| Project modules the benchmark code imports or that case files name as producing scripts (`phase3/…`, `semiconductor_sim/common/…`) | 26 |
| `paper/` and the top-level release files | the rest |

- **How the module list was found.** Every module of `benchmark/wp2`, `benchmark/wp4` and `benchmark/tools` was imported, and the repository files it loaded were recorded. A static import pass then covered the three producing scripts.
- **Personal data and secrets.** I searched every benchmark file, including the gzip transcripts, for:
  - the user's e-mail address and GitHub name;
  - token patterns (`sk-ant-`, `ghp_`, `github_pat_`, `xox…`, AWS keys, `Bearer …`);
  - session URLs and OAuth or API-key variables.

  None was found. The only path-like match is `/Users/user/...` inside DEVSIM's own manual text, which some agents
  extracted into their workspaces. The agent transcripts' init events carry no account or organization identifiers; their session ids belong to ephemeral headless sessions.

## Checks

| Check | Result |
|---|---|
| Full test suite in the repository | 71 passed, 1 skipped (the snapshot-only check) |
| The tarball unpacked in an unrelated directory, no git, tests run there | 70 passed, 2 skipped (the two checks that need the research history) |
| Registered evaluations rerun inside that snapshot (`evaluate_wp2.py`, `evaluate_wp4.py`, `evaluate_wp6.py --stage 1`) and the asset script | all outputs byte-identical to the released files (SHA-256) |
| Manuscript built inside that snapshot | 13 pages, no undefined references or citations |
| Exporter (`export_cases.py`) in the repository | reproduces every case file and the manifest byte-identically |
| Figures | rendered and inspected; palette validated for colour-vision deficiency (3 categorical slots, all pairs), with marker shape as a second encoding |

**References.** Each entry was checked against one of three sources:
- a publisher or index record: IEEE/ISC, Crossref via the project literature matrix, CaltechAUTHORS, university repositories, Scilit, BibBase;
- the arXiv abstract page, for the six 2025–2026 agent and TCAD-agent papers;
- the PDF the user supplied.

Three works are deliberately not cited as evidence:
- Polsky & Rimshans (1981): mentioned only as cited by Palm & Van de Wiele;
- Laux & Fischetti (1985): unread;
- Ward, Odeh & Cohen (1990): checked, but dropped because the manuscript no longer makes a claim about it.

## Manuscript claims checked against the evidence

While writing I corrected these statements before they reached the PDF:

- **The femtoampere figure.** The low end is 1.7e-15 A/cm, i.e. sub-attoampere per µm. The text now says "as small as 1e-15 A per centimetre of width".
- **"Weak inversion".**
  - ψ̂ = 12 and 16 are depletion for dev1, and ψ̂ = 12 is depletion for dev2.
  - The text now says "below threshold, from depletion to weak inversion".
  - Tier C is described by its gate voltages, not as strong inversion.
- **G_E and G_M.** I had first called them two forms of the output. They are in fact one coupled solve done with two linear solvers: equilibrated SuperLU (G_E) and UMFPACK through DEVSIM (G_M).
  - κ_G is now given with its definition.
  - The "forward error ≈ κη" sentence was replaced by a general statement. κ_G was computed with the computed, not the true, G, so it does not explain the magnitude quantitatively.
- **The drain-side flux-sum mechanism.** I had offered a floating-point reading of it. It was unverified and its magnitude did not fit, so it was removed. The text keeps only what was observed: the error does not shrink under refinement, and the source and energy forms agree to ≤ 1e-9.
- **P4's independence.** It is independent in implementation, not in authorship; the same author wrote both paths.
- **The agent study.**
  - Verification: 31 of 32 statements describe a second current measure, and one a sheet-density estimate. This is recorded in the new post hoc file `benchmark/wp6/results/verification_coding_stage1.json`.
  - "All failing paths in every run" was corrected to: the contact current in 33/33 runs and the AC analysis in 30/33.

## Deviations

1. **`benchmark/tools/export_cases.py` now refuses to run without the research history** (`d7b39f9`).
   - Testing the snapshot showed that the exporter deleted the case files before it checked its sources against git history, so outside the research repository it deleted the cases and then failed.
   - It now checks both base commits before changing anything. Its output is unchanged, and the determinism test passes in the repository.
   - The two history-dependent tests in `test_benchmark.py` are skipped where the history is absent.
2. **One new post hoc coding file** for the agent study, `verification_coding_stage1.json`. It is labelled post hoc and decides nothing.
3. **The release snapshot includes 26 project modules and 89 earlier result files at their original paths**, beyond `benchmark/`. Without them the case provenance and the reference code would not resolve. They are released under the same licences.

## What the user needs to do (the remaining release steps)

1. Give the author line: names, affiliations and, optionally, ORCID iDs. I will then:
   - fill it into `paper/main.tex`, `release/CITATION.cff`, `release/.zenodo.json`, `release/LICENSE` and `release/LICENSE-DATA.md`;
   - rebuild the PDF and the snapshot;
   - create the final branch `release/v1.0.0`.
2. Review the manuscript, in particular the "Use of AI" statement, which needs the author's confirmation.
3. Create the public repository, for example `mosfet-lowcurrent-dd-benchmark`, and turn on the Zenodo GitHub integration for it.
4. Push the release branch as its `main` branch, either yourself or by adding the new repository to this session so I can do it. Then publish a GitHub release `v1.0.0`; Zenodo archives it and mints the DOI.
5. Send me the DOI and URL. I will add them to `CITATION.cff`, the READMEs and the manuscript as a follow-up release.

## Limits

The manuscript is a complete draft, not a submission:
- no journal has been chosen;
- it has had no independent review;
- the "first" claims are deliberately absent, and its novelty statements are limited to "not found in the literature we reviewed".

Every limit of WP1–WP6 carries over:
- agreement between computational paths only;
- a narrow physics envelope;
- DEVSIM is the only external code scored;
- the agent study tested one configuration.
