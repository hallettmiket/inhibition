# Dance with Inhibition

Four independent computational approaches to finding an **inhibitor** of human
Pin1, and the choreography that organizes them.

A murmurent **choreography** is a problem statement many people attack, each
with their own **approach**. This site documents the four approaches, the shared
substrate they all depend on, and the decisions taken along the way.

!!! quote "The design stance"
    **There is no authoritative cross-approach numeric ranking.** Vina affinity
    (kcal/mol, lower better) and gnina `CNNaffinity` (dimensionless, *higher*
    better) are not the same axis. A non-covalent complex and a covalent adduct
    with the leaving group removed are not the same physical quantity. Forcing
    a merge would be the easiest way to produce a confident wrong answer, so
    the integration phase **presents** rather than merges.

## The four approaches

| | Approach | Seed | Search | Covalent? | Inhibition proxy |
|---|---|---|---|---|---|
| **T_1** | de novo generation | none — the pocket | DiffSBDD diffusion | non-covalent | weak (reversible occupancy) |
| **T_2** | ATRA neighborhood | ATRA | CReM, degree-bounded | non-covalent | weak (reversible occupancy) |
| **T_3** | R-group decoration | sulfopin (core + warhead fixed) | REINVENT 4 `libinvent` | covalent | strong (covalent Cys113) |
| **T_4** | warhead x R-group | sulfopin (core fixed) | combinatorial | covalent | strong (covalent Cys113) |

## Where to start

- **New to the project** — [Framing](overview/framing.md), then [Approaches](approaches/index.md).
- **Want to know why something is the way it is** — [Decisions](decisions/index.md).
- **About to make a judgment call** — [Runbooks](runbooks/index.md).
- **Wondering what is blocked right now** — [Status](overview/status.md).
- **Running something** — [Operations](ops/running.md).

## Source of truth

Files in the repo, always. This site, and the Streamlit GUI, are **views** over
them — see [D0008](decisions/index.md). A published result must be reproducible
by reading the repo, without standing up either.
