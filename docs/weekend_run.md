# The weekend run — 2026-08-07 23:00 → Monday 07:00

Operational runbook for the sweep → elevation → 100 ns pipeline running
unattended over the weekend of 2026-08-08/09. Written for @tt8804 to read cold
at 06:00 Saturday.

**Goal:** a defensible top-5 for synthesis, presented Monday 09:00.

---

## The 30-second check

```bash
cat /data/lab_vm/modifiable/inhibition/weekend/supervisor_status.txt
```

```
supervisor alive   2026-08-07 23:38:27
sweeps claimed     9   (worklist units left: 123)
sweep results      5
MD runs claimed    1
workers            6 of 6
```

- **`supervisor alive`** older than ~3 minutes → the supervisor is stuck. See
  *Recovery*.
- **`workers`** below 6 → a GPU is idle; the supervisor should refill it within
  two minutes on its own.
- **`worklist units left`** should fall by roughly 4 every 21 minutes.

---

## What is running

| session | GPU | job |
|---|---|---|
| `sweep1` `sweep2` `sweep3` `sweep5` | 1, 2, 3, 5 | 10 ns attack-geometry sweeps |
| `md6` `md7` | 6, 7 | 100 ns runs on survivors |
| `wsuper` | — | restarts dead workers, every 120 s |
| `wkeep` | — | restarts a dead supervisor, every 120 s |
| `gui` | — | Streamlit on 127.0.0.1:8899 |

GPUs **0 and 4 are never used** — standing instruction, and 4 carries another
group's job.

Every worker is parented to the tmux server (PID 3047317), so an operator
disconnect cannot kill them.

### The two stages

**Sweep** — 10 ns of unbiased MD per (molecule, mode), ~21 min measured. Reads
`weekend/worklist.txt`, claims a unit with `mkdir` (atomic — two workers cannot
take the same one), writes to
`00_outputs/blacksmith/attack_sweep/attack_sweep_<N>.csv`.

**Elevation** — a worker polls every 10 minutes, takes the best unclaimed
survivor, and runs 100 ns (~3.2 h measured). When the sweep worklist empties the
sweep workers **convert themselves to MD**, taking elevation from 2 GPUs to 6.

Both refuse to start work they cannot finish before Monday 07:00.

### Timeline

| when | what |
|---|---|
| Sat ~10:00 | sweep list exhausted, all 6 GPUs on 100 ns |
| Mon 03:00 | last 100 ns run may start |
| Mon 07:00 | hard deadline; everything stops |

---

## What counts as a survivor

**At least one sustained episode of attack geometry** — a run of ≥100 ps in
which the warhead is both in the 2.8–4.2 Å window and at a mechanism-appropriate
angle. Ordered by `frac_attack_ready`.

This replaced `frac_attack_ready > 0` on the first night, at @tt8804's insistence
and correctly. The first molecule swept, `t4_e0b03662d460`, scored 0.2% — **one
frame out of 501**, a 20 ps touch — and the old rule called it a survivor. The
debounced visit count in the same row already said `n_visits = 0`; the two
criteria contradicted each other and the weaker one was being used.

**The bar is still open.** @tt8804's position is that it should be far higher.
The counter-argument is that literal 100% occupancy rejects every molecule ever
measured here — the best of the six 100 ns runs was 55.2% attack-ready and the
other five spanned 0–7.3%, references lower still. A bound ligand rattles; it
samples attack geometry in bursts. **Set the bar from the distribution once ~120
sweeps are in**, not from either intuition.

Columns that inform that decision, per sweep:

| column | meaning |
|---|---|
| `frac_attack_ready` | fraction of frames in attack geometry |
| `n_visits` | episodes lasting ≥100 ps — the debounced count |
| `n_visits_raw` | every touch, including single frames |
| `median_episode_ps` | how long a typical visit lasts |
| `start_attack_ready` | did it begin in attack geometry |

A large `n_visits_raw` with `n_visits` near zero means the molecule is skimming
the boundary rather than approaching — visible in the very first results, where
one molecule touched 25 times and sustained none.

---

## Recovery

**A worker died.** Nothing to do; the supervisor refills within 120 s.

**The supervisor is stuck** (status file not updating, process alive):

```bash
pkill -f weekend_supervisor.sh          # the wrapper restarts it in ~10 s
```

**Everything is gone** (tmux server died):

```bash
tmux new-session -d -s wsuper "while :; do bash /data/lab_vm/modifiable/inhibition/weekend_supervisor.sh; sleep 10; done"
tmux new-session -d -s wkeep  "while :; do bash /data/lab_vm/modifiable/inhibition/weekend_keepalive.sh; sleep 120; done"
```

The supervisor rebuilds every worker from there. Claims already made are
preserved, so nothing is swept twice — **but a molecule claimed by a killed
worker is never retried**. To re-run one, remove its claim:

```bash
rmdir /data/lab_vm/modifiable/inhibition/weekend/claims_sweep/<molecule>__r<rank>
```

> **Never edit a shell script while its workers are running.** Bash reads
> scripts by byte offset; an in-place edit corrupts the running process. This
> happened to the supervisor on the first night — it kept running and silently
> stopped updating its status. Write a new file and repoint the supervisor
> instead; `weekend_worker_2.sh` and `_3.sh` exist for exactly this reason.

---

## What this run can and cannot tell you

**Can:** which molecules present their warhead to Cys113 in a sustained way
under explicit-solvent dynamics, and which of those hold position for 100 ns.

**Cannot:**

- **Say anything about reaction rate.** Everything here is the pre-reaction
  complex. Nothing models bond formation; that is FEP, on a handful of
  candidates, later.
- **Validate the ranking that chose the worklist.** The molecules swept were
  picked by `conditional_eb`, which has never passed a convergence test — see
  `docs/prereg_score_selection.md`, where convergence is disqualifying and comes
  first. The sweep is a fair test of the molecules it was given; whether it was
  given the right molecules is a separate, open question (#42, #43).
- **Support ranking on the sweep's own readout.** Agreement between 10 ns and
  100 ns is ρ = +0.60 on attack-readiness — the quantity the sweep ranks on —
  against ρ = +0.83 for mere proximity, which is the number quoted in older
  docs. The pre-registered reading table says +0.60 is good enough to reject the
  bottom and **not** to order the middle (#34).
- **Rest on unreplicated residence numbers.** Every 100 ns run is n = 1, and the
  field's standard for this measurement is replicates (#33 — and note the three
  molecules that came unbound are exactly the three simulated without ions).

---

## Files

| path | what |
|---|---|
| `weekend/supervisor_status.txt` | the health check |
| `weekend/worklist.txt` | `<molecule> <pose_rank>` per line, priority order |
| `weekend/claims_sweep/`, `claims_md/` | one directory per claimed unit |
| `weekend/logs/` | per-worker logs, `supervisor.log` |
| `00_outputs/blacksmith/attack_sweep/` | sweep results |
| `00_outputs/blacksmith/md_residence_3ikd/` | 100 ns results |

---

## Composition of the worklist

125 T4 molecules, 132 units (a molecule contributes one unit per mode worth
sweeping):

| class | molecules |
|---|---|
| bdhi_c4 | 40 |
| acrylamide | 40 |
| bdhi_c5 | 40 |
| others | 5 |

T4 only, by decision on the first night: T3 carries a **single warhead class**
across all 4,065 of its molecules, docks each into essentially one position, and
sits in range 85% of the time against T4's 44%. That is the profile of a small
rigid set scoring well for reasons unrelated to being a good inhibitor, and it is
under review (#40).

Ordered by `conditional_eb` — empirical-Bayes shrinkage — which replaced the
Wilson lower bound because the bound systematically demoted minority binding
modes: it pushes uncertain estimates toward zero, and uncertainty scales with the
mode's population (#42).
