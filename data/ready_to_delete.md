# Ready to delete — retired reference-file versions

Integer versioning: the largest integer is the newest. Superseded versions are
listed here rather than edited in place, so a pipeline run always pins an exact
version and an old result stays reproducible against the file it actually used.

**Do not delete anything on this list until no committed manifest references
it.** A manifest records the SHA-256 of the reference files a run consumed;
deleting a file that a manifest names makes that run unverifiable.

| File | Superseded by | Date retired | Why | Safe to delete |
|---|---|---|---|---|
| `reference/pin1_covalent_cys113_anchors_1.csv` | `_2` | 2026-07-27 | The Reddi 2023 sulfamate-acetamide anchor was `UNVERIFIED`; Figure 5 of the paper resolved it into two verified anchors (4d, 4g). Anchor count 6 → 7, verified 4 → 6. | not yet — no runs consumed it, but keep until M0 is git-tagged |
| `reference/warhead_classes_1.csv` | `_2` | 2026-07-27 | `sulfamate_acetamide` moved `UNVERIFIED` → `VERIFIED`; `sulfonate_acetamide` (compound 4a) added. Enumerable classes 1 → 3. | not yet — same reason |

| `reference/pin1_reference_binders_1.csv` | `_2` | 2026-07-27 | Two rows were attributed to ChEMBL **CHEMBL3391**, which is *Threonine—tRNA ligase 1*, **not Pin1** (Pin1 is CHEMBL2288) — verified directly against the ChEMBL API. Both removed. Eight new actives added from CHEMBL2288 + literature; the resolved Reddi 4d/4g structures folded in. 11 rows → 18. | not yet — no runs consumed it |

| `decoys/decoys_covalent_1.csv` | `_2` | 2026-07-27 | Only 32 of 302 decoys (10.6%) carried any electrophile, so ~90% could not be covalently docked at all — gnina needs a reactive atom to bond to. The surviving comparison would have been "electrophiles vs inert molecules", which docking wins trivially. `_2` requires every decoy to carry a warhead motif. | not yet — no runs consumed it |

| `append_only/.../D4_1.parquet` | `D4_2` | 2026-07-27 | Enumerated with a broken `snar_chloroazine` fragment that omitted the chlorine, so the warhead-validity gate correctly killed all 198 of that class. `D4_2` has the corrected fragment. | not yet |

## Note on the governed data root

Large derived outputs live under `/data/lab_vm/append_only/inhibition/` and are
never overwritten. When an approach is re-run and supersedes an earlier
integer-versioned output there, record the retired version here too — the
append-only tree grows, so this file is how it gets pruned safely.

Nothing under `/data/lab_vm/immutable/inhibition/` is ever listed here. That
tree is read-only source (6VAJ, model weights, fragment DBs, decoys) and is not
retired by this process.
