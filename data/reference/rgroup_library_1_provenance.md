# R-group library — provenance

`rgroup_library_1.csv`, built by
`approaches/t4_combinatorial/build_rgroup_library.py`.

## How it was derived

**Aryl / heteroaryl groups: by frequency, not by taste.** The 33 ring
systems are the most common small aromatic Murcko ring systems in the ChEMBL
pool pinned in `config/sources.lock.json`
(`chembl_pool.csv`, 20858 molecules). Ring systems were kept when
they had 5-12 heavy atoms and at least
one aromatic atom — small enough to be decoration rather than a second scaffold.

This is deliberate. D0004 rejected inheriting the prior run's 444-member library
because it predated the reference set. Replacing it with substituents chosen by
intuition would be no better grounded, only differently arbitrary. "These are
the substituents medicinal chemists actually use" is a claim a reviewer can
check by re-running this against the same hash-pinned file.

**Linkers: explicit, and a design choice.** 6 linkers connect the
sulfolane nitrogen to the aryl group. This is not a frequency question, and the
verified anchors are informative: Sulfopin uses neopentyl and Reddi 4g uses
cyclohexylmethyl — both CH2-linked, which is why `1C` is included and why direct
attachment and longer chains bracket it.

| id | SMILES | rationale |
|---|---|---|
| `direct` | (none) | aryl bonded straight to N |
| `1C` | `C` | the anchors' own linker |
| `2C` | `CC` | one atom longer |
| `3C` | `CCC` | reach into a deeper sub-pocket |
| `alpha_Me` | `C(C)` | branch at the alpha carbon |
| `gem_diMe` | `C(C)(C)` | quaternary, conformationally restricted |

## Size

198 R-groups = 33 aryls x 6 linkers (minus any that
failed to construct). The enumerated library is this multiplied by the enumerable
warhead classes, so `gates.yaml` keeps `library_size: null` until the
enumeration stage pins the real number.

## Known limitations

- **Frequency is not suitability.** A substituent common in ChEMBL is common
  across all targets; nothing here is Pin1-specific. The pocket is shallow and
  solvent-exposed, and a chemist may well want groups this method will not
  surface.
- **Murcko scaffolds discard substituents on the ring**, so `phenyl` stands in
  for every substituted phenyl. That keeps the set small and orthogonal to the
  linker axis, but it means fluorinated and methylated variants are absent.
  Adding them is a CSV edit.
- The pool was filtered to MW 150-700 when fetched, so ring systems appearing
  only in larger molecules are underrepresented.
