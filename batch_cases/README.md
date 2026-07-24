# Table 2 batch cases

These macros implement the Table 2 cases that can be represented with the
currently available `../phantoms` and `../models/merged` tetrahedral models.

Run any case with:

```bash
/private/tmp/tetcal-build/TetCal -m batch_cases/<case>.in -o <case>.out
```

Each case uses:

- source material ID `-4` (`source_water`)
- `/spec/IAEA iaea_spectra <nuclide>`
- `/run/beamOn 10000`

## Implemented

| File | Table 2 case |
|---|---|
| `Preinjection_I-131_1.in` | I-131, Posture 1, L block shield + L-block vial |
| `Preinjection_I-131_2.in` | I-131, Posture 2, Vial Assembly |
| `Preinjection_I-131_3.in` | I-131, Posture 2, Vial |
| `Injection_I-131_1.in` | I-131, Posture 3, Syringe assembly + syringe source |
| `Injection_I-131_2.in` | I-131, Posture 3, syringe |
| `Injection_I-131_3.in` | I-131, Posture 2, Vial assembly |
| `Preinjection_Tc-99m_1.in` | Tc-99m, Posture 1, L block shield + L-block vial |
| `Preinjection_Tc-99m_2.in` | Tc-99m, Posture 2, Vial Assembly |
| `Preinjection_Tc-99m_3.in` | Tc-99m, Posture 2, Vial |
| `Injection_Tc-99m_1.in` | Tc-99m, Posture 3, Syringe assembly + syringe source |
| `Injection_Tc-99m_2.in` | Tc-99m, Posture 3, Syringe |
| `Preinjection_F-18_1.in` | F-18, Posture 1, L block shield + L-block vial |
| `Preinjection_F-18_2.in` | F-18, Posture 2, Vial Assembly |
| `Preinjection_F-18_3.in` | F-18, Posture 2, Vial |
| `Injection_F-18_1.in` | F-18, Posture 3, Syringe assembly + syringe source |
| `Injection_F-18_2.in` | F-18, Posture 3, syringe |

## Not yet implemented

The Table 2 administration-site and post-injection patient cases need additional
geometry/source setup before a defensible macro can be generated:

- IV administration site in arm
- swallowed liquid/capsule source
- patient at 3 ft / 6 ft
- patient at delayed time points

Those cases need a source-region definition and, for patient-distance cases, a
second source/person geometry or an external-source approximation.
