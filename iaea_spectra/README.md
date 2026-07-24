# IAEA LiveChart per-decay spectra

Generated from IAEA LiveChart API v1 on demand for I-131, Tc-99m, and F-18.

Files:
- `*_per_decay_discrete.tsv`: ungrouped discrete gamma, x-ray, Auger/conversion electron, and alpha lines.
- `*_per_decay_discrete_aggregated.tsv`: discrete lines grouped by particle, radiation code, shell, and energy.
- `*_beta_branches.tsv`: beta-/beta+/electron-capture branches from `decay_rads`.
- `*_beta_continuum.tsv`: binned beta spectra from `bin_beta`; `dn_dE` is already per decay.
- `raw/*.csv`: raw API responses.
- `summary.tsv`: yield/integral checks by nuclide and particle.

IAEA units are keV and percent. The generated files also include MeV and yield per decay.
For beta continuum, `dn_dE_per_decay_per_keV` integrates over keV to the beta yield per decay.
Tc-99m is selected from the `99tc` API response by `p_energy = 142.6836 keV`.
