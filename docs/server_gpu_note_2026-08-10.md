# Server GPU usage — what happened, 2026-08-10 ~01:00

Dalwinder Singh flagged that this account was occupying **all eight GPUs on
biodatsci**, asked for a lower number, and said the system reboots in the morning.
This is what was found and done, so the reply to him can be specific.

## What was actually holding the GPUs

Nine `gmx mdrun` jobs under `twu383`, spread across all 8 cards. Only three were
the overnight controls; the other six came from earlier sessions. The reason
usage never dropped was not the jobs themselves but four **relauncher loops**
that refilled a card the moment one freed:

| loop | what it did |
|---|---|
| `maxout_supervisor.sh` in `while :; sleep 10` | kept every GPU saturated |
| `maxout_keepalive.sh` in `while :; sleep 120` | restarted work that stopped |
| `weekend_worker.sh sweep {2,3,5}`, `weekend_worker_2.sh sweep 1` | pulled the next molecule per card |
| `rerun69_worker.sh {6,7}` | same, found on a second pass |

The name is the diagnosis: `maxout` exists to max the box out. Nothing would have
drained on its own.

## What was done

All relauncher loops stopped. **Running trajectories were left alone**; the three
controls kept their GPUs. Footprint went **8 of 8 → 3 of 8**, and cannot grow
again because nothing is left to launch new work.

**Cost, stated plainly.** Stopping the `weekend_worker` loops also killed the four
candidate 100 ns runs those loops owned — I expected the children to survive and
they did not. All four hold `prod.cpt`, so they resume rather than restart:

| candidate | stopped at |
|---|---|
| `t4_be879c67f567` | 69.3 / 100 ns |
| `t4_3d79e382fd35` | 54.5 / 100 ns |
| `t4_5bbb2b4414b8` | 46.8 / 100 ns |
| `t4_4d2a30f4345d` | 37.6 / 100 ns |

Two partial 10 ns sweeps (`t4_ecfea7489e0a`, `t4_c21cef59b29c`) went the same way
and are cheap to re-run.

## The guardrail gap

`shared/compute.py` caps **CPU** workers at 50. There is no equivalent cap on
**GPUs**, and `CLAUDE.md` states the rule as a convention — *"GPUs are shared with
other users: check `nvidia-smi` before taking a card"* — which a supervisor loop
does not read. A worker pool that keeps every card busy satisfies every automated
check in the repo.

Worth fixing as a real guard rather than a line in a doc: a maximum-cards constant
the workers must honour, refusing to launch past it.

## For the reply to Dalwinder

- Cause: supervisor/keepalive loops, now stopped.
- Now at 3 of 8, peaking at 4 when the last control starts.
- Everything finishes tonight; nothing runs for days.
- Nothing was left that can re-occupy cards after the reboot.
