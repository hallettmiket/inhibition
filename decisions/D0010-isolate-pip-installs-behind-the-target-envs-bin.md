---
id: D0010
title: Isolate pip installs behind the target env's bin
date: 2026-07-27
status: accepted
approach: shared
decided_by: '@mhallett'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/setup_envs.sh
evidence:
  - "REINVENT's install.py shells out to a bare `pip`, not sys.executable -m pip"
  - 'it resolved to base conda pip and installed reinvent + ~30 packages into base'
  - 'pandas was upgraded 2.x -> 3.0.5, breaking streamlit 1.55.0 and anndata 0.12.8'
  - 'pydantic was upgraded <2 -> 2.13.4, breaking anaconda-cloud-auth 0.1.3'
  - 'repaired by pinning pandas<3 and pydantic<2 back into base'
  - 'plan also stale: REINVENT needs python >=3.11 (built 3.10), and install.py first positional is processor_type not the dependency set'
runbook: null
---

## Context

The five-env design exists so three incompatible torch builds cannot collide.
That protects against *conda* cross-contamination, but not against a third-party
installer shelling out to a bare `pip`.

REINVENT's `install.py` builds `cmd = ["pip", "install", ...]`. Invoked as
`$ENV/bin/python install.py`, the Python is the env's but the `pip` is whatever
PATH resolves first — here base conda's. REINVENT and roughly thirty chemistry
packages went into the base environment, and pandas and pydantic were upgraded
out from under base tooling that pins them.

Nothing failed loudly. The install reported success; the env stayed empty; the
damage was in a different environment entirely.

## Decision

Any third-party installer that may shell out to `pip` is run with the target
env's `bin` FIRST on PATH, and the result is verified by importing from that
env rather than trusting the installer's exit code.

Base conda is treated as off-limits. It belongs to the user's other projects and
is not part of this choreography's environment contract.

## Consequences

`setup_envs.sh` now sets PATH explicitly for the REINVENT build and asserts the
CLI runs afterward. Base was repaired to pandas 2.3.3 / pydantic 1.10.26; the
remaining `pip check` failures there (requests, urllib3, tqdm, numpy) predate
this and are not ours.

The wider lesson matches D0009: verify the artifact you produced. An installer
exit code of 0 says the installer ran, not that the package landed where you
intended.
