# Framing

## The problem

Find a molecule that **inhibits** human Pin1 (peptidyl-prolyl cis-trans
isomerase NIMA-interacting 1) — critically an inhibitor, not an activator.

A good candidate should ideally:

- **(i)** pass structural-alert and drug-likeness filtering;
- **(ii)** have intrinsic electrophilic reactivity in a safe, precedented range;
- **(iii)** be synthetically accessible;
- **(iv)** have a structurally well-resolved binding mode, not merely a high
  docking score.

Some approaches add **(v)**: engage Cys113 via a chemically valid covalent
mechanism.

Conditions are **approach-specific** and are tracked as parameters, not applied
as universal gates. Not every approach can resolve every condition, and that is
fine — what matters is that the coverage gaps are visible.

## Inhibition versus activation is not resolved computationally

No approach here separates inhibition from activation. Occupancy of the
catalytic site is the working proxy — and **the proxy is not equally strong
across approaches**:

- **T_3 / T_4** engage Cys113 covalently, as sulfopin does. That is a defensible
  inhibition mechanism.
- **T_1 / T_2** rely on reversible pocket occupancy, which is weaker: a
  non-productive binder need not inhibit catalysis.

This asymmetry is surfaced per approach and displayed in the GUI, not buried.

## Conditions x approaches

| Condition | T_1 | T_2 | T_3 | T_4 |
|---|---|---|---|---|
| (i) alerts + drug-likeness | ✓ | ✓ | ✓ | ✓ |
| (ii) reactivity window | — | ~ removes electrophiles | partial | ✓ |
| (iii) synthetic accessibility | ✓ | ✓ | ✓ | ✓ |
| (iv) resolved binding mode | ✓ | ✓ | ✓ | ✓ |
| (v) covalent at Cys113 | — | — | ✓ | ✓ |
| inhibition-proxy strength | weak | weak | strong | strong |

T_4 is the only approach touching all five.

## The data model

Each approach emits a candidate frame `D^i` — rows are candidates, columns are
attributes — keyed on **canonical SMILES**, the join key across approaches.
Each delivers its **top 10** to the integration phase.

**Stamp, do not delete.** A candidate failing a tier's filter is stamped with
`rejected_at` and skips only the more expensive downstream tiers. It is never
removed, so later reweighting can resurrect it. Gates throttle compute, not
membership.
