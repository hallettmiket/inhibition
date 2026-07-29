# Status

!!! info "Generated page"
    Rendered at build time from the repo's source of truth. Edit the underlying file, not this page.

## Pinned external sources

| Source | Pin |
|---|---|
| `crem_fragment_db_chembl33` | `e863c5de768df66d…` |
| `crem_fragment_db_enamine` | `ce0b380903b82553…` |
| `diffsbdd_checkpoint` | `07f86764bf569aaf…` |
| `diffsbdd_repo` | `5d0d38d16c8932a0…` |
| `gnina` | `3340c1f49cd3c7c8…` |
| `receptor_6vaj` | `820fd5969131bef8…` |
| `reinvent4_repo` | `04de385d33f95e97…` |
| `reinvent_priors` | `03e6cbe8a53e59a4…` |

## Open questions

The choreography's honest limits, kept next to its results rather than in a file nobody opens.

*(warhead library unavailable in this environment: warhead_classes_2.csv missing column(s): ['adduct_attachment_smarts', 'has_leaving_group'])*

### Not yet acquired

- **CReM fragment DB** — the radius variant must be chosen and pinned before T_2 can enumerate; a different radius is a different neighbourhood definition, not a tuning knob.
- **Decoy set** — built by the enrichment gate, which has not run.

### Not yet verified

- **Byun 2023 BDHI fragment** — SI-only. The BDHI *class* is verified (PubChem CID 21983498), which is what the reactivity window needs, so this blocks enumeration rather than the window.
