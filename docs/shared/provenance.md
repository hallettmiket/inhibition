# Provenance

Four layers, each answering a different question. Together they mean nobody has
to reconstruct what happened from memory.

| Layer | Answers | Where |
|---|---|---|
| Decision record | *why* we chose this | `decisions/` |
| Runbook | *how* to choose again | `docs/runbooks/` |
| Manifest | *what* a run consumed | `manifest.json`, in the append-only tree |
| Reference provenance | where the data came from | `data/reference/.provenance.md` |

## Declarative acquisition

Every external input is declared in `config/sources.yaml` — URL, destination,
which approaches need it, license, and why it matters. Staging a fresh machine:

```bash
python -m shared.sources stage
```

Observed hashes and resolved commits land in `config/sources.lock.json`. First
acquisition **observes**; every run after **enforces**.

!!! warning "Config and pins never share a file"
    An earlier version wrote pins back into `sources.yaml` via `yaml.safe_dump`
    and erased every explanatory comment in it. Config is hand-authored; pins
    are machine-written. See [D0007](../decisions/index.md).

## Run manifests

Each stage writes a `manifest.json` recording run id, git commit **and whether
the tree was dirty**, config hashes, input and output hashes, and tool versions.

That last point carries more weight than it looks. "The numbers changed" is only
diagnosable if you can tell whether the *code*, the *config*, or the *input*
changed. Without a manifest all three look identical.

!!! danger "Dirty tree means provisional"
    `git.dirty == true` means the recorded commit does not fully describe the
    code that ran. Outputs carrying it are provisional.

## One-shot judgment is fine — write the runbook

Choosing a PDB entry or reading a compound out of a figure needs judgment, and
using the model for that is the point. But the *procedure* has to survive, so
the next case does not depend on whoever was in the room. That is what
[runbooks](../runbooks/index.md) are for.
