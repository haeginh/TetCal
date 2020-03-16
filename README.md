# ODDCal
(synchronized with ODDCal_Vox)

* __Developed Functions__: multiple organs, RBM
* __odd.mac__: contains macro to run 1E7 particles for six directions

*__MACROS__ _(see sample.in)_
```
/run/initialize             (required before any run)
/odd/organ     "1 2" [or 1] (organ IDs, double quotation is required for multiple sources) 
/odd/organ     "0 1 2"      (for RBM IDs, put 0 in the beginning) 
/run/beamOn    10000        (nps)
```
