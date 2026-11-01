# Archive checkpoint

The branch checkpoint is preserved in two non-destructive artifacts outside the
checkout:

- `/workspace/phase3_codex_next_question_artifacts/mosfet-tcad-phase3-next-scientific-question.zip`
- `/workspace/phase3_codex_next_question_artifacts/mosfet-tcad-phase3-next-scientific-question.gitbundle`

The companion `SHA256SUMS.txt` records the hashes of the ZIP, bundle, bundle
verification log, and fresh-restoration logs.  The ZIP contains the tracked
source snapshot plus `ARCHIVE_INFO.md`, `ENVIRONMENT.txt`, `GIT_INFO.txt`,
`REPRODUCE.md`, the source manifest and a transport copy of the Git bundle.

Verification performed:

1. `unzip -t` returned exit code 0.
2. A fresh extraction passed `MANIFEST_SHA256.txt` and reproduced the seven-row
   finite-bias analysis with the same extrapolated errors.
3. `git bundle verify` passed; cloning the bundle into an empty directory,
   checking out this branch, and `git fsck --full --strict` all passed.

These are workspace download paths.  If the interface does not expose binary
attachments directly, the GitHub branch and the restore commands in
`REPRODUCE.md` remain the recovery path; no inaccessible path is treated as a
successful attachment transfer.
