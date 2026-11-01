# Independent discovery evidence corrections

This note records small corrections required before the next research phase.
Historical reports remain unchanged; this note is the authoritative correction
layer.

## Candidate S1 checksum

The file
`verification/candidates/phase2b_v1_0/stage1.json` was hashed directly in the
restored checkout on 2026-11-01:

`93244365320f005eb61ae242f5416cbb4cf9a9af92b1e63917370e336adc0bd0`

The full 64-character value above is the value to use in handoffs and archive
manifests. The corresponding S2 candidate is
`7883dfdd80fe5fc3b768ed3e21088f6b212682dd74048c033a935dd1162afeac`.
All tracked handoff and candidate-manifest references in this restored tree now
match the direct hash. Any shortened or malformed value in an external handoff
must be treated as superseded by this note.

## Gummel and Scharfetter–Gummel records

The inherited records were checked character-for-character against the archived
Phase 2B citation table, the machine-readable source list, and the DOI URL
forms stored with the project. The bibliographic records are:

| Work | DOI retained | Verification status |
|---|---|---|
| A. A. Gummel, “A self-consistent iterative scheme for one-dimensional steady state transistor calculations,” *IEEE Transactions on Electron Devices* 11(10), 455–465 (1964) | `10.1109/T-ED.1964.1536099` | DOI string and citation agree across the archived records; no replacement string is justified |
| D. L. Scharfetter and H. K. Gummel, “Large-signal analysis of a silicon Read diode oscillator,” *IEEE Transactions on Electron Devices* 16(1), 64–77 (1969) | `10.1109/T-ED.1969.16754` | DOI string and citation agree across the archived records; no replacement string is justified |

Direct DOI/Crossref resolution was attempted from this environment and rejected
by the configured proxy with HTTP 403. Accordingly, this phase does not claim a
fresh publisher-side metadata retrieval. It corrects the provenance statement
and prevents an unverified alternate DOI from being introduced.

## Phase 3 regression label

The combined-suite logs in the independent restoration directory are labelled
as Phase 3 regression results. They are separate from the historical startup
validation logs and preserve the four skips caused by absent official S1/S2
golden files.

## Network-energy assumptions

The zero-bias network identity is conditional on a fixed equilibrium state,
constant mobility, the stated SG linearization, symmetric edge conductances,
and the unit-voltage normalization (`source = 0`, `drain = 1`). It is an
algebraic consistency check on one discrete network. Agreement with the finite
volume path is a separate model-verification result and is not implied by the
identity itself.
