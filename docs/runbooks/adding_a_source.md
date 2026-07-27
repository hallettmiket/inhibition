# Runbook — bringing a new external file into the choreography

**Use this when:** the choreography needs any file it does not produce itself —
a structure, a model checkpoint, a fragment database, a decoy set, a reference
CSV.

**Why it matters:** an input acquired by a command someone typed once is an
input nobody can reacquire. The choreography is meant to be re-run and
re-parameterized; every external byte has to be declarable, fetchable, and
pinned, or "reproducible" is a claim rather than a property.

---

## Procedure

### 1. Decide whether it is a *source* or *reference data*

Two different destinations, and the distinction is about size and stability:

| | Source | Reference data |
|---|---|---|
| Examples | PDB entry, model weights, fragment DB | the frozen binder set, warhead classes |
| Size | large | small (KB) |
| Lives in | `immutable/inhibition/` — **not git** | `data/reference/` — **in git** |
| Declared in | `config/sources.yaml` | committed directly, with `.provenance.md` |
| Versioned by | hash pin in `sources.lock.json` | integer suffix (`_1`, `_2`) |

If it is small, hand-curated, and needs review in a PR, it is reference data.
If it is large or machine-fetched, it is a source.

### 2. For a source — declare it before fetching it

Add to `config/sources.yaml`:

```yaml
  my_new_input:
    description: "one line — what it is"
    kind: http              # http | git | generated
    url: "https://..."      # null if not yet known
    dest: /data/lab_vm/immutable/inhibition/<subdir>/<file>
    required_by: [t2]       # which approaches break without it
    license: MIT            # if there is one
    notes: >
      Why it matters, and anything a future reader would otherwise have to
      rediscover — a variant choice, a known gotcha, a stale upstream note.
```

Then:

```bash
python -m shared.sources stage --only my_new_input
```

First acquisition **observes** the hash and records it in
`config/sources.lock.json`. Every run after **enforces** it.

!!! danger "Never write pins into sources.yaml"
    `sources.yaml` is hand-authored and never written by code. An earlier
    version of `shared/sources.py` round-tripped it through `yaml.safe_dump` to
    record hashes and erased all 22 comments in the file — before it had ever
    been committed, so git could not restore it. See
    [D0007](../decisions/index.md).

### 3. If it cannot be fetched yet, declare it as pending anyway

A source with `url: null` or `kind: generated` is reported as **pending** rather
than being absent. This is deliberate: "we know we need this and have not got
it" is information, and it surfaces in the [Status](../overview/status.md) page
and the GUI's Open Questions panel. An input that appears from nowhere later is
worse than one declared missing now.

### 4. For reference data — commit it, and write provenance

- Name it with an integer suffix: `<name>_1.csv`. Never edit in place; a change
  makes `_2`, and `_1` is retired via `data/ready_to_delete.md`.
- Record in `data/reference/.provenance.md`: where each row came from
  (PubChem CID, ChEMBL ID, DOI), what was **excluded** and why, and any
  verification performed.
- **Never invent a value.** Unresolvable entries are marked `UNVERIFIED` and the
  loaders refuse them. See
  [resolving unverified structures](resolving_unverified_structures.md).
- Add a status column with real tiers if "verified" is not binary for your data.
  A status field nothing enforces is a comment — make a loader check it.

### 5. Watch for CSV escaping

Chemical data breaks naive CSV constantly. Both of these have already bitten:

- **Commas inside chemical names** — `1,4-naphthoquinone`,
  `3-bromo-4,5-dihydroisoxazole`. Quote the field, or use a semicolon.
- **Commas inside SMARTS** — `[c,C]([Br])=[N]` is a legal pattern and an illegal
  unquoted CSV cell. Quote every SMARTS and SMILES field defensively.

Round-trip through `pandas.read_csv` before committing.

### 6. Record the decision if it is consequential

Choosing *which* checkpoint, *which* fragment-DB radius, or *which* PDB entry is
a decision, not a mechanical step. Write a record in `decisions/` — the choice
of a variant is exactly the sort of thing that looks arbitrary six months later.

### 7. Confirm it is reachable from a clean state

```bash
python -m shared.sources check
```

Should report your source as pinned. If a colleague on a fresh machine cannot
run `stage` and get the same bytes, the job is not finished.

---

## Worked example — the DiffSBDD checkpoint, 2026-07-27

**What went wrong first.** The checkpoint was fetched with an ad-hoc `curl` in a
terminal, alongside an ad-hoc `git clone` of DiffSBDD and an ad-hoc `mkdir` of
the governed tree. None of it was in the repo. It worked — and was
unreproducible.

**The fix.** Declared in `config/sources.yaml` with URL, destination, license
(MIT), the Zenodo record number (8183747), and a note that it is inference-only.
Re-staged through `shared.sources`, which found the file already present, hashed
it, and pinned `07f86764…` into `sources.lock.json`.

**Verification that the pin has teeth.** Corrupted the recorded hash, re-ran,
and confirmed it raised rather than silently re-downloading:

```
SourceError: receptor_6vaj: ... hash 820fd5969131bef8 != pinned 0000000000000000.
An input changed underneath a pinned run; investigate before overwriting anything.
```

**A detail worth copying.** `_fetch_http` uses `curl --fail`. Without it, curl
writes a 404 HTML page to the destination and exits 0 — producing a
"checkpoint" that is really an error page, which then fails much later,
somewhere far less obvious.

---

## Failure modes

| Failure | Consequence | Guard |
|---|---|---|
| Ad-hoc fetch | input unreacquirable | declare in `sources.yaml` |
| No hash pin | upstream changes silently | lockfile enforcement |
| `curl` without `--fail` | 404 page saved as a checkpoint | `--fail` in `_fetch_http` |
| Pins written into config | comments erased | pins live in the lockfile |
| Unquoted SMARTS/SMILES in CSV | parse error, or worse, silent misparse | quote defensively, round-trip test |
| Editing reference data in place | old results unreproducible | integer versioning + `ready_to_delete.md` |
| Undeclared missing input | appears from nowhere later | declare as `pending` |
