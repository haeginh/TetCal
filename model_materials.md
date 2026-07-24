# Model Material Definitions

This file documents the material assignments used for the additional TET models
under `../models/merged`.

The generated `.material` files are:

- `../models/merged/L_block_vial.material`
- `../models/merged/L_block_shield.material`
- `../models/merged/syringe.material`
- `../models/merged/syringe_assembly.material`
- `../models/merged/vial.material`
- `../models/merged/vial_assembly.material`

All six files intentionally contain the same complete material table so any
merged model can be loaded independently or together with other model fragments.

## Assignments

| ID | CAD group | Material used | Density (g/cm3) | Basis |
|---:|---|---|---:|---|
| -1 | `lead` | elemental lead | 11.35 | NIST elemental Pb |
| -2 | `lead_glass` | NIST `Glass, Lead` | 6.22 | NIST XCOM compound |
| -3 | `not_defined` | polypropylene | 0.90 | assumed plastic for unlabeled syringe assembly part |
| -4 | `source` | liquid water | 1.00 | NIST XCOM water |
| -5 | `stainless_steel` | AISI 304 stainless steel | 8.00 | nominal 304 composition |
| -6 | `syringe` | polypropylene | 0.90 | common disposable syringe barrel material |
| -7 | `tungsten` | elemental tungsten | 19.30 | NIST elemental W |
| -8 | `vial` | NIST `Glass, Borosilicate (Pyrex)` | 2.23 | NIST XCOM compound |

## Composition Notes

Fractions in the `.material` files are mass fractions, matching the existing
MRCP material-file convention.

- Lead, tungsten: pure elements from NIST elemental-material constants.
- Lead glass, borosilicate glass, water: direct NIST XCOM compound/mixture
  mass fractions.
- Polypropylene: repeat unit `(C3H6)n`; mass fractions are H 0.143716 and
  C 0.856284. Density is set to 0.90 g/cm3, within the common commercial PP
  range.
- Stainless steel: AISI 304 nominal composition was used with Fe as balance:
  C 0.0008, Si 0.0100, P 0.00045, S 0.00030, Cr 0.1900, Mn 0.0200,
  Fe 0.68345, Ni 0.0950.

## Sources

- NIST X-Ray Mass Attenuation Coefficients, Table 1: elemental densities and
  material constants.
- NIST X-Ray Mass Attenuation Coefficients, Table 2: compound/mixture densities
  and mass fractions for water, lead glass, and borosilicate glass.
- AISI/ASTM 304 stainless steel nominal composition ranges; midrange Cr/Ni and
  maximum minor-element fractions were used, with Fe as balance.
