# Running things

## Principle

**Build the code, review it, then run.** The approaches are not run
incrementally as they are written. What *does* run during the build phase is
substrate preparation and validation — and every such step is recorded code
with a manifest, never ad-hoc shell.

## Long jobs go under tmux

CReM enumeration, DiffSBDD inference, REINVENT RL, covalent docking, xtb, and
MM-GBSA/MD all run tmux-backed and checkpointed, so a killed job resumes rather
than restarts.

```bash
tmux new-session -d -s dwi_<stage> "cd ~/repos/inhibition && nice -n 19 <cmd>"
tmux capture-pane -p -t dwi_<stage> | tail
```

!!! warning "Chain stages with `;` not `&&`"
    An `&&` chain silently skips everything after the first failure. This
    already happened once: REINVENT 4 failed and `amber_md` never ran, while
    the session still looked busy.

## Substrate commands

```bash
python -m shared.sources stage        # acquire + pin external inputs
python -m shared.sources check        # verify pins only
python -m shared.receptor_prep        # prepare receptor + boxes  (M1)
python -m shared.decisions check      # validate the decision log
mkdocs serve                          # this site, live, on :8000
```

## Resource conventions

`nice -n 19` on everything. 8× A100-80GB available; use them when it makes
sense. The dashboard and docs bind to localhost — reach them over an SSH tunnel.

## Milestones

| | Done when | Status |
|---|---|---|
| **M0** | repo, envs, reference frozen and tagged | mostly ✓ |
| **M1** | prepared receptor + both boxes | ✓ |
| **M2** | covalent protocol pinned, one dock end-to-end | ✗ |
| **M3** | enrichment gate PASS/FAIL token written | ✗ |
| **M4a–d** | the four shortlists | ✗ |
| **M5** | integration GUI | ✗ |
| **M6** | adversary sign-off + Oracle capture | ✗ |

If M3 **fails**, M4a–d still proceed with docking demoted to a label. The
choreography does not stall on a failed gate.
