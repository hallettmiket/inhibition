---
id: D0079
title: Pose splitting and per-mode ranking are accepted; the selection that consumed them was not per-mode, and four questions remain open
date: 2026-08-11
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/nac_screen_v2.py
  - scripts/attack_sweep.py
  - shared/mode_key.py
  - shared/mode_ranking.py
  - shared/mode_assets.py
evidence:
  - '@tt8804 on the ranking view: "this looks beautiful ... ranking and pose splitting is good"'
  - 'pose splitting: 8,096 modes over 5,747 molecules; 1,362 molecules carry more than one mode'
  - 'the ranking is per mode and pooled; selection for simulation was not -- 233 of 239 modes sent were mode 0, and 242 of 242 molecules contributed their mode 0'
  - '5 modes ranking FIRST in their warhead class were never sent; 50 modes in a class top-10 were never sent'
  - 'the collision that hid it: mode 0 is the bare ident in attack_sweep and _m0 in rank_v2, so merge(on="ident") drops every simulated row without erroring'
  - 'fixed: shared/mode_key.py keys on (parent_ident, mode); attack_sweep reads mode from the pose property and always writes _m<mode>'
  - 'pose_rank - 1 == mode holds for all 1,751 exported poses checked, so the old derivation was latent rather than active'
  - 'unresolved: #47 measures every warhead class with wet-lab anchoring ranking LAST, on a criterion whose difficulty varies by class'
---

# What is accepted

**Pose splitting.** Docked poses are clustered on the reactive atom's position and
the direction the warhead faces — never on energy, never on the NAC geometry the
criterion then measures. 500 poses per molecule become 8,096 modes over 5,747
molecules; 1,362 molecules carry more than one. The representative of a mode is
the medoid of its best-anchored quartile, not its best-scoring member, because
argmax of a noisy score recovered the crystal pose 6.7% of the time against 26.7%
for a typical member.

**Per-mode ranking.** A mode is a candidate row. Modes are scored independently
and ranked in one pooled list, within a warhead class.

**The view.** `modes.html` is the first of two GUIs — the ranked list and the
poses, before anything is simulated — and `combined.html` is the second, the
sweep and 100 ns results for what was. Every mode carries its rank, its pose, and
its siblings' ranks beside it.

# What is fixed

The selection that consumed the ranking was **not** per-mode: it took mode 0,
once per molecule, for 242 of 242 molecules. It stayed invisible because a mode
was written two ways — `t4_x` in `attack_sweep`, `t4_x_m0` in `rank_v2` — so the
obvious join dropped exactly the rows that had been simulated, without erroring.

`shared/mode_key.py` is now the only place a mode key is built and the key is
`(parent_ident, mode)`. `attack_sweep` reads the mode from the pose's own
property rather than deriving it from `pose_rank`, and always writes `_m<mode>`.
`md_residence` records which pose it ran.

# What is NOT resolved, and blocks

These are not bookkeeping. Each changes what the ranking means.

1. **#47 — the criterion ranks the chemistry we most believe in LAST.** Every
   warhead class with wet-lab anchoring — chloroacetamide (crystal + known
   actives), sulfamate_acetamide and sulfonate_acetamide (measured k) — sits at
   the bottom, medians 0.000–0.001. The classes leading the shortlist, bdhi_c5
   and bdhi_c4, have **no measured Pin1 activity at all**, and the best-scoring
   control is the promiscuous naphthoquinone reference. A large part is mechanism
   bias: `SN2_ANGLE_MIN = 150°` is a far narrower target than
   `PERPENDICULAR_MAX_OFF_NORMAL = 30°`, and the isotropic nulls do not equalise
   them. **This is the deepest open problem in the ranking and nothing in this
   record addresses it.**
2. **#53's audit half.** The code is fixed; the science is not. Whether mode 0 was
   chosen deliberately is unestablished, the 5 class-leading modes are still
   unswept, and the shortlist has not been re-derived under per-mode selection —
   so it is not known whether 2.2.0's ordering is an ordering of mode 0s.
3. **#44 — the pose cloud is not persisted.** `nac_screen_v2` docks into a
   `tempfile.mkdtemp` and deletes it, so a mode's individual poses no longer
   exist; only medoids do. Confirmed by a filesystem-wide search. Mode membership
   cannot be re-derived, and the ranking view can never show the cloud behind a
   mode.
4. **#39 — the Sulfopin pose has never been shown to be right.** The one molecule
   whose true binding position is known, and our pose has not been validated
   against it.

# The decision

Pose splitting and per-mode ranking stand. The selection defect is closed. **The
stage is not closed**, because #47 puts the ranking's meaning in question in a way
that no amount of correct plumbing repairs, and because it is not yet known
whether the 2.2.0 shortlist is a shortlist of mode 0s.
