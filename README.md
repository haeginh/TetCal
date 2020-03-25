# TetCal
(synchronized with VoxelVersion - except for surface source)

* __Under Development__: visualization
* __Developed Functions__: External, Internal (volume/surface), Hadronic, dose arrangement (tallying), effective dose calculation
* __Dose Arrangement__: doses can be organized in [phantom name].dose file (similar to tallying cells in MCNP)
```
[dose ID] [dose name] [organ ID list]
```
* __Effective Dose Calculation__: Implemented in SetEffectiveDose() in RunAction.cc. When the organ IDs were changed, please modify the function

*__MACROS__ _(see sample.in)_
```
/run/initialize       (required before any run)
/gun/particle  gamma  (particle)
/gun/energy    1 MeV  (energy)
/external/dir  AP     (idealized external irradiation geometries)
/internal/source "1 2" [or 1] (internal soruce organ IDs, double quotation is required for multiple sources) 
/internal/surface "1 2" [or 1] (internal soruce organ IDs, double quotation is required for multiple sources) 
/run/beamOn    10000  (nps)
```
