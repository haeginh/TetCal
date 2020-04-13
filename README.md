# TetCal
(Lung installed, for adult(new ID) lung)

* __Under Development__: visualization
* __Developed Functions__: External, Internal (volume/surface), Hadronic, dose arrangement (tallying), effective dose calculation
* __Dose Arrangement__: doses can be organized in [phantom name].dose file (similar to tallying cells in MCNP)
```
[dose ID] [dose name] [organ ID list]
```
* __Effective Dose Calculation__: Implemented in SetEffectiveDose() in RunAction.cc. When the organ IDs were changed, please modify the function

*__MACROS__ _(see sample.in)_
```
/lung/volchk ICRP-AM.lungVol (optional, calculates the lung branch volume and print the data file)
/lung/sampling 10000         (optional, for source sampling #needs to be checked)
/lung/BB-bas 465.345 mm3     (required, set BB-bas volume)
/lung/BB-sec 285.198 mm3
/lung/bb-sec 1477.61 mm3
/run/initialize       (required before any run)
/gun/particle  gamma  (particle)
/gun/energy    1 MeV  (energy)
/external/dir  AP     (idealized external irradiation geometries)
/internal/source "1 2" [or 1] (internal soruce organ IDs, double quotation is required for multiple sources) 
/internal/surface "1 2" [or 1] (internal soruce organ IDs, double quotation is required for multiple sources) 
/run/beamOn    10000  (nps)
```
