# `pocket_landmarks_1.csv` — provenance

The landmark residues the contact-space pose splitter measures against.

## Where they come from

`exp/14_residue_selection` — a **greedy non-redundant** ranking, not a distance
cutoff and not a literature list. Each residue is added in turn by how much
contact-matrix variance it explains that the residues already chosen do not.

Two properties made this set adoptable:

* **the order is identical across five independent dockings** (Spearman 1.000),
  so it is a property of the pocket rather than of one pose cloud;
* **the top 15 span 90%** of the contact matrix's variance.

## Waters are excluded

`A:40:HOH` ranked 13th on the raw list and is **not** in this file. A landmark
that is modelled inconsistently between structures is not a landmark, and the
exclusion is applied once here rather than being a rule every reader has to
remember.

## How to use it

Read through `shared.reference_set.latest_reference("pocket_landmarks")`, never
by naming a version. `shared/pose_contacts.landmark_residues(n)` does this.

Regenerate with `exp/14_residue_selection/run_all.py`; write a new
`pocket_landmarks_2.csv` rather than editing this one.
