
**************************************************************
 Geant4 version Name: geant4-11-01-patch-02 [MT]   (15-June-2023)
  << in Multi-threaded mode >> 
                       Copyright : Geant4 Collaboration
                      References : NIM A 506 (2003), 250-303
                                 : IEEE-TNS 53 (2006), 270-278
                                 : NIM A 835 (2016), 186-225
                             WWW : http://geant4.org/
**************************************************************

================================================================================
	/Users/hanh4/workspace/phantomData/MRCP/MRCP_AM was implemented in this CODE!!   
================================================================================
  Opening TETGEN node (vertex points: x y z) file '/Users/hanh4/workspace/phantomData/MRCP/MRCP_AM.node'
  Opening TETGEN elements (tetrahedron with node No.) file '/Users/hanh4/workspace/phantomData/MRCP/MRCP_AM.ele'
  Opening material file '/Users/hanh4/workspace/phantomData/MRCP/MRCP_AM.material'

 Organ ID   # of Tet  vol [cm3]  d [g/cm3]   mass [g]	organ/tissue
--------------------------------------------------------------------------------
      100       2682      8.381      1.036      8.683	Adrenal_left
      200       2661      8.381      1.036      8.683	Adrenal_right
      300       4675      0.022      1.031      0.022	ET1(0-8)
      301       4671      0.087      1.031      0.090	ET1(8-40)
      302       4675      0.027      1.031      0.028	ET1(40-50)
      303       5910     10.952      1.031     11.291	ET1(50-Surface)
      400       8615      0.141      1.000      0.141	ET2(-15-0)
      401       8617      0.378      1.031      0.390	ET2(0-40)
      402       8605      0.095      1.031      0.098	ET2(40-50)
      403       8610      0.048      1.031      0.049	ET2(50-55)
      404       8611      0.095      1.031      0.098	ET2(55-65)
      405      11327     27.942      1.031     28.808	ET2(65-Surface)
      500       6493      0.082      1.050      0.086	Oral_mucosa_tongue
      501        768      0.022      1.050      0.024	Oral_mucosa_moth_floor
      600       5884      0.022      1.050      0.023	Oral_mucosa_lips_and_cheeks
      700       2848     10.052      1.031     10.364	Trachea
      800       7423      0.025      1.000      0.025	BB(-11--6)
      801       7419      0.030      1.031      0.031	BB(-6-0)
      802       7427      0.050      1.031      0.052	BB(0-10)
      803       7424      0.126      1.031      0.130	BB(10-35)
      804       7426      0.025      1.031      0.026	BB(35-40)
      805       7419      0.051      1.031      0.052	BB(40-50)
      806       7422      0.051      1.031      0.052	BB(50-60)
      807       7423      0.051      1.031      0.053	BB(60-70)
      808       7469      2.693      1.031      2.777	BB(70-surface)
      810          0      0.000      1.000      0.000	bb(-6--4)
      811          0      0.000      1.031      0.000	bb(-4-0)
      812          0      0.000      1.031      0.000	bb(0-4)
      813          0      0.000      1.031      0.000	bb(4-12)
      814          0      0.000      1.031      0.000	bb(12-20)
      815          0      0.000      1.031      0.000	bb(20-25)
      816          0      0.000      1.031      0.000	bb(25-surface)
      900       8184      1.419      1.060      1.505	Blood_in_large_arteries_head
      910       6763      6.551      1.060      6.944	Blood_in_large_veins_head
     1000     199208    182.248      1.060    193.183	Blood_in_large_arteries_trunk
     1010     194299    418.889      1.060    444.022	Blood_in_large_veins_trunk
     1100      79897     30.628      1.060     32.466	Blood_in_large_arteries_arms
     1110      64756    157.837      1.060    167.307	Blood_in_large_veins_arms
     1200     144062    102.685      1.060    108.846	Blood_in_large_arteries_legs
     1210     129955    367.665      1.060    389.725	Blood_in_large_veins_legs
     1300       4514     83.684      1.904    159.335	Humeri_upper_cortical
     1400       1591    118.158      1.233    145.689	Humeri_upper_spogiosa
     1500       1036     34.907      0.981     34.244	Humeri_upper_medullary_cavity
     1600       4206     55.978      1.904    106.582	Humeri_lower_cortical
     1700       1604     45.888      1.109     50.890	Humeri_lower_spongiosa
     1800        794     38.121      0.981     37.397	Humeri_lower_medullary_cavity
     1900      13312    143.644      1.904    273.498	Ulnae_and_radii_cotical
     2000       4193    139.748      1.109    154.981	Ulnae_and_radii_spongiosa
     2100       3982     23.442      0.981     22.996	Ulnae_and_radii_medullary_cavity
     2200      55485     95.341      1.904    181.529	Wrists_and_hand_bones_cortical
     2300      26623    107.238      1.109    118.927	Wrists_and_hand_bones_spongiosa
     2400       6931     25.342      1.904     48.252	Clavicles_cortical
     2500       3511     38.943      1.157     45.057	Clavicles_spongiosa
     2600     182290    298.566      1.904    568.469	Cranium_cortical
     2700      94869    327.960      1.165    382.073	Cranium_spongiosa
     2800       9093    133.166      1.904    253.548	Femora_upper_cortical
     2900       4287    367.318      1.125    413.232	Femora_upper_spongiosa
     3000       1039     26.550      0.981     26.045	Femora_upper_medullary_cavity
     3100       7647    161.639      1.904    307.761	Femora_lower_cortical
     3200       2935    336.927      1.109    373.652	Femora_lower_spongiosa
     3300       1539     83.771      0.981     82.179	Femora_lower_medullary_cavity
     3400      24141    281.855      1.904    536.651	Tibiae_cortical
     3500      11014    560.332      1.109    621.408	Tibiae_spongiosa
     3600       1807     81.361      0.981     79.815	Tibiae_medullary_cavity
     3700      34425    123.363      1.904    234.882	Ankles_and_foot_cortical
     3800      15057    390.095      1.109    432.615	Ankles_and_foot_spongiosa
     3900      10806     40.377      1.904     76.877	Mandible_cortical
     4000       4795     44.286      1.271     56.287	Mandible_spongiosa
     4100      20889    211.447      1.904    402.595	Pelvis_cortical
     4200      11688    552.785      1.121    619.672	Pelvis_spongiosa
     4300     163309    193.696      1.904    368.797	Ribs_cortical
     4400      84967    390.898      1.170    457.351	Ribs_spongiosa
     4500      43903    117.297      1.904    223.333	Scapulae_cortical
     4600      29007    130.450      1.201    156.670	Scapulae_spongiosa
     4700     113671     54.592      1.904    103.943	Cervical_spine_cortical
     4800      52165     75.229      1.049     78.915	Cervical_spine_spongiosa
     4900     188736    152.017      1.904    289.440	Thoracic_spine_cortical
     5000      83302    322.637      1.070    345.222	Thoracic_spine_spongiosa
     5100      78135     98.764      1.904    188.047	Lumbar_spine_cortical
     5200      35089    263.163      1.108    291.584	Lumbar_spine_spongiosa
     5300      52281     57.941      1.904    110.320	Sacrum_cortical
     5400      25006    186.083      1.033    192.224	Sacrum_spongiosa
     5500      19001      5.247      1.904      9.991	Sternum_cortical
     5600       9944     59.001      1.041     61.420	Sternum_spongiosa
     5700      30691     51.257      1.099     56.331	Cartilage_costal
     5800      41596     74.671      1.099     82.063	Cartilage_discs
     6100      16686   1457.627      1.041   1517.390	Brain
     6200       1178      8.152      0.953      7.769	Breast_left_adipose_tissue
     6300        767      5.074      1.021      5.180	Breast_left_glandular_tissue
     6400       1312      8.152      0.953      7.769	Breast_right_adipose_tissue
     6500        783      5.074      1.021      5.180	Breast_right_glandular_tissue
     6600       2163      0.037      1.060      0.039	Eye_lens_sensitive_left
     6601       3207      0.179      1.060      0.189	Eye_lens_insensitive_left
     6700      24539      1.012      1.100      1.113	Cornea_left
     6701       3227      0.300      1.025      0.308	Aqueous_left
     6702       6392      5.938      1.031      6.122	Vitreous_left
     6800       2175      0.037      1.060      0.039	Eye_lens_sensitive_right
     6801       3206      0.179      1.060      0.189	Eye_lens_insensitive_right
     6900      24536      1.012      1.100      1.113	Cornea_right
     6901       3226      0.300      1.025      0.308	Aqueous_right
     6902       6391      5.938      1.031      6.122	Vitreous_right
     7000       3709     10.052      1.031     10.364	Gall_bladder_wall
     7100       1733     56.311      1.030     58.000	Gall_bladder_contents
     7200      13472      1.720      1.037      1.784	Stomach_wall(0-60)
     7201      13457      1.150      1.037      1.193	Stomach_wall(60-100)
     7202      13447      5.794      1.037      6.008	Stomach_wall(100-300)
     7203      14006    178.675      1.037    185.286	Stomach_wall(300-surface)
     7300       6770    240.385      1.040    250.000	Stomach_contents
     7400     107655     14.028      1.037     14.547	Small_intestine_wall(0-130)
     7401     108211      2.184      1.037      2.264	Small_intestine_wall(130-150)
     7402     107956      5.489      1.037      5.692	Small_intestine_wall(150-200)
     7403     102273    810.122      1.037    840.096	Small_intestine_wall(200-surface)
     7500     105813     51.286      1.040     53.337	Small_intestine_contents(-500-0)
     7501      68911    285.253      1.040    296.663	Small_intestine_contents(centre--500)
     7600       5033      2.962      1.037      3.071	Ascending_colon_wall(0-280)
     7601       5045      0.215      1.037      0.223	Ascending_colon_wall(280-300)
     7602       5752    112.473      1.037    116.634	Ascending_colon_wall(300-surface)
     7700       2556     52.885      1.040     55.000	Ascending_colon_content
     7800       4209      3.850      1.037      3.993	Transverse_colon_wall_right(0-280)
     7801       4192      0.278      1.037      0.289	Transverse_colon_wall_right(280-300)
     7802       5675     72.971      1.037     75.671	Transverse_colon_wall_right(300-surface)
     7900       2076     91.346      1.040     95.000	Transverse_colon_contents_right
     8000       3269      2.723      1.037      2.824	Transverse_colon_wall_left(0-280)
     8001       3235      0.198      1.037      0.205	Transverse_colon_wall_left(280-300)
     8002       4072     74.179      1.037     76.924	Transverse_colon_wall_left(300-surface)
     8100       1611     38.462      1.040     40.000	Transverse_colon_content_left
     8200       5852      2.680      1.037      2.779	Descending_colon_wall(0-280)
     8201       5887      0.196      1.037      0.203	Descending_colon_wall(280-300)
     8202       6284    112.774      1.037    116.946	Descending_colon_wall(300-surface)
     8300       2866     33.654      1.040     35.000	Descending_colon_content
     8400       6071      4.292      1.037      4.451	Sigmoid_colon_wall(0-280)
     8401       6068      0.312      1.037      0.324	Sigmoid_colon_wall(280-300)
     8402       6004     46.784      1.037     48.515	Sigmoid_colon_wall(300-surface)
     8500       3040     72.115      1.040     75.000	Sigmoid_colon_contents
     8600        851     38.550      1.037     39.976	Rectum_wall
     8700      15528    367.116      1.051    385.839	Heart_wall
     8800       8946    481.132      1.060    510.000	Blood_in_heart_chamber
     8900       9947    154.168      1.053    162.338	Kidney_left_cortex
     9000       5835     36.429      1.053     38.359	Kidney_left_medulla
     9100       1079      7.267      1.053      7.652	Kidney_left_pelvis
     9200       7378    158.159      1.053    166.542	Kidney_right_cortex
     9300       3176     37.381      1.053     39.362	Kidney_right_medulla
     9400        786      7.495      1.053      7.892	Kidney_right_pelvis
     9500      12572   2226.415      1.060   2360.000	Liver
     9700       8585   1315.366      0.415    545.877	Lung(AI)_left
     9900       9944   1573.160      0.415    652.861	Lung(AI)_right
    10000       6820     15.455      1.032     15.949	Lymphatic_nodes_ET
    10100       6820     15.455      1.032     15.949	Lymphatic_nodes_thoracic
    10200       2370      5.339      1.032      5.510	Lymphatic_nodes_head
    10300      55761    126.166      1.032    130.204	Lymphatic_nodes_trunk
    10400       4755     10.678      1.032     11.019	Lymphatic_nodes_arms
    10500       4712     10.678      1.032     11.019	Lymphatic_nodes_legs
    10600     131490   1143.640      1.050   1200.822	Muscle_head
    10700     396774  14134.733      1.050  14841.469	Muscle_trunk
    10800     168276   2708.029      1.050   2843.431	Muscle_arms
    10900     254621  10372.285      1.050  10890.899	Muscle_legs
    11000      10079      1.851      1.037      1.919	Oesophagus_wall(0-190)
    11001      10340      0.099      1.037      0.103	Oesophagus_wall(190-200)
    11002       8309     48.007      1.037     49.783	Oesophagus_wall(200-surface)
    11003       5932     21.990      1.040     22.870	Oesophagus_contents
    11300       3149    166.314      1.044    173.631	Pancreas
    11400       2835      0.603      1.031      0.622	Pituitary_gland
    11500       1730     17.088      1.031     17.618	Prostate
    11600     384611   1039.002      0.939    975.623	RST_head
    11700    1428065  11902.982      0.939  11176.900	RST_trunk
    11800     313652   1650.524      0.939   1549.842	RST_arms
    11900     537504   4803.129      0.939   4510.138	RST_legs
    12000       5947     42.721      1.031     44.045	Salivary_glands_left
    12100       6717     42.721      1.031     44.045	Salivary_glandss_right
    12200      86263    238.044      1.089    259.230	Skin_head_insensitive
    12201      43176      7.778      1.089      8.470	Skin_head_sensitive(50-100)
    12300     110988   1167.258      1.089   1271.144	Skin_trunk_insensitive
    12301      55549     35.278      1.089     38.418	Skin_trunk_sensitive(50-100)
    12400     149397    528.657      1.089    575.708	Skin_arms_insensitive
    12401      75136     17.303      1.089     18.843	Skin_arms_sensitive(50-100)
    12500     209125   1156.993      1.089   1259.965	Skin_legs_insensitive
    12501     105071     34.701      1.089     37.790	Skin_legs_sensitive(50-100)
    12600       6127     36.811      1.031     37.952	Spinal_cord
    12700       2516    215.472      1.060    228.400	Spleen
    12800       7296     18.872      2.688     50.727	Teeth
    12801      11936      0.041      1.040      0.043	Teeth_retention_region
    12900       2596     17.884      1.041     18.617	Testis_left
    13000       2424     17.884      1.041     18.617	Testis_right
    13100       5112     25.130      1.031     25.909	Thymus
    13200       9199     22.217      1.051     23.351	Thyroid
    13300       1757     19.993      1.050     20.993	Tongue_upper(food)
    13301       9474     51.954      1.050     54.552	Tongue_lower
    13400       2132      3.016      1.031      3.109	Tonsils
    13500       2853      8.544      1.031      8.809	Ureter_left
    13600       6074      7.539      1.031      7.773	Ureter_right
    13700       6155     47.866      1.040     49.781	Urinary_bladder_wall_insensitive
    13701       2969      1.267      1.040      1.318	Urinary_bladder_wall_sensitive(75-118)
    13800       1438    192.308      1.040    200.000	Urinary_bladder_content
    14000      20910    140.165      0.001      0.140	Air_inside_body
Visualization Manager instantiating with verbosity "warnings (3)"...
Visualization Manager initialising...
Registering graphics systems...

You have successfully registered the following graphics systems.
Registered graphics systems are:
  ASCIITree (ATree)
  DAWNFILE (DAWNFILE)
  G4HepRepFile (HepRepFile)
  RayTracer (RayTracer)
  VRML2FILE (VRML2FILE)
  gMocrenFile (gMocrenFile)
  TOOLSSG_OFFSCREEN (TSG_OFFSCREEN)
  TOOLSSG_OFFSCREEN (TSG_OFFSCREEN, TSG_FILE)
  OpenGLImmediateQt (OGLIQt, OGLI)
  OpenGLStoredQt (OGLSQt, OGL, OGLS)
  Qt3D (Qt3D)
  TOOLSSG_QT_GLES (TSG_QT_GLES, TSGQt, TSG)

Registering model factories...

You have successfully registered the following model factories.
Registered model factories:
  generic
  drawByAttribute
  drawByCharge
  drawByOriginVolume
  drawByParticleID
  drawByEncounteredVolume

Registered models:
  None

Registered filter factories:
  attributeFilter
  chargeFilter
  originVolumeFilter
  particleFilter
  encounteredVolumeFilter

Registered filters:
  None

You have successfully registered the following user vis actions.
Run Duration User Vis Actions: none
End of Event User Vis Actions: none
End of Run User Vis Actions: none

Some /vis commands (optionally) take a string to specify colour.
"/vis/list" to see available colours.

==========================================================================================
G4TaskRunManager :: Using G4ThreadPool...
==========================================================================================

userDetector->Construct() start.

   Phantom name                /Users/hanh4/workspace/phantomData/MRCP/MRCP_AM TET phantom
   Phantom size                557.708 * 291.260 * 1760.047 mm3
   Phantom box position (min)  -286.958 mm, -136.028 mm, -878.602 mm
   Phantom box position (max)  270.750 mm, 155.232 mm, 881.445 mm
   Number of tetrahedrons      8234101

worldLogical is registered to the default region.
physicsList->Construct() start.
physicsList->CheckParticleList() start.
physicsList->setCut() start.
=======================================================================
======                 Electromagnetic Physics Parameters      ========
=======================================================================
LPM effect enabled                                 1
Enable creation and use of sampling tables         0
Apply cuts on all EM processes                     0
Use combined TransportationWithMsc                 Disabled
Use general process                                0
Enable linear polarisation for gamma               0
Enable photoeffect sampling below K-shell          1
Enable sampling of quantum entanglement            0
X-section factor for integral approach             0.80000
Min kinetic energy for tables                      100.00000 eV 
Max kinetic energy for tables                      100.00000 TeV
Number of bins per decade of a table               20
Verbose level                                      1
Verbose level for worker thread                    0
Bremsstrahlung energy threshold above which 
  primary e+- is added to the list of secondary    100.00000 TeV
Bremsstrahlung energy threshold above which primary
  muon/hadron is added to the list of secondary    100.00000 TeV
Lowest triplet kinetic energy                      1.00000 MeV
Enable sampling of gamma linear polarisation       0
5D gamma conversion model type                     0
5D gamma conversion model on isolated ion          0
Livermore data directory                           epics_2017
=======================================================================
======                 Ionisation Parameters                   ========
=======================================================================
Step function for e+-                              (0.20000, 0.01000 mm)
Step function for muons/hadrons                    (0.10000, 0.05000 mm)
Step function for light ions                       (0.10000, 0.02000 mm)
Step function for general ions                     (0.10000, 0.00100 mm)
Lowest e+e- kinetic energy                         100.00000 eV 
Lowest muon/hadron kinetic energy                  1.00000 keV
Use ICRU90 data                                    1
Fluctuations of dE/dx are enabled                  1
Type of fluctuation model for leptons and hadrons  Urban
Use built-in Birks satuaration                     0
Build CSDA range enabled                           0
Use cut as a final range enabled                   0
Enable angular generator interface                 1
Max kinetic energy for CSDA tables                 1.00000 GeV
Max kinetic energy for NIEL computation            1.00000 MeV
Linear loss limit                                  0.01000
Read data from file for e+e- pair production by mu 0
=======================================================================
======                 Multiple Scattering Parameters          ========
=======================================================================
Type of msc step limit algorithm for e+-           2
Type of msc step limit algorithm for muons/hadrons 0
Msc lateral displacement for e+- enabled           1
Msc lateral displacement for muons and hadrons     1
Urban msc model lateral displacement alg96         1
Range factor for msc step limit for e+-            0.08000
Range factor for msc step limit for muons/hadrons  0.20000
Geometry factor for msc step limitation of e+-     2.50000
Safety factor for msc step limit for e+-           0.60000
Skin parameter for msc step limitation of e+-      3.00000
Lambda limit for msc step limit for e+-            1.00000 mm
Use Mott correction for e- scattering              1
Factor used for dynamic computation of angular 
  limit between single and multiple scattering     1.00000
Fixed angular limit between single 
  and multiple scattering                          3.14159 rad
Upper energy limit for e+- multiple scattering     100.00000 MeV
Type of electron single scattering model           0
Type of nuclear form-factor                        1
Screening factor                                   1.00000
=======================================================================
======                 Atomic Deexcitation Parameters          ========
=======================================================================
Fluorescence enabled                               1
Directory in G4LEDATA for fluorescence data files  fluor
Auger electron cascade enabled                     1
PIXE atomic de-excitation enabled                  0
De-excitation module ignores cuts                  1
Type of PIXE cross section for hadrons             Empirical
Type of PIXE cross section for e+-                 Livermore
=======================================================================

### ===  Deexcitation model UAtomDeexcitation is activated for 1 region:
          DefaultRegionForTheWorld  1  1  0
### ===  Auger flag: 1
### ===  Ignore cuts flag:   1

phot:  for gamma SubType=12 BuildTable=0
      LambdaPrime table from 200.000000 keV to 100.000000 TeV in 174 bins 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
 LivermorePhElectric : Emin=0.000000 eV  Emax=100.000000 TeV  SauterGavrila Fluo

compt:  for gamma SubType=13 BuildTable=1
      Lambda table from 100.000000 eV  to 1.000000 MeV, 20 bins/decade, spline: 1
      LambdaPrime table from 1.000000 MeV to 100.000000 TeV in 160 bins 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
    LivermoreCompton : Emin=0.000000 eV  Emax=1.000000 GeV Fluo
        KleinNishina : Emin=1.000000 GeV Emax=100.000000 TeV Fluo

conv:  for gamma SubType=14 BuildTable=1
      Lambda table from 1.021998 MeV to 100.000000 TeV, 20 bins/decade, spline: 1
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
Livermore5DConversion : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedTsai

Rayl:  for gamma SubType=11 BuildTable=1
      Lambda table from 100.000000 eV  to 150.000000 keV, 20 bins/decade, spline: 0
      LambdaPrime table from 150.000000 keV to 100.000000 TeV in 176 bins 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
   LivermoreRayleigh : Emin=0.000000 eV  Emax=100.000000 TeV  CullenGenerator

msc:  for e-  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  GoudsmitSaunderson : Emin=0.000000 eV  Emax=100.000000 MeV Nbins=120 100.000000 eV  - 100.000000 MeV
          StepLim=SafetyPlus Rfact=0.080000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm
        WentzelVIUni : Emin=100.000000 MeV Emax=100.000000 TeV Nbins=120 100.000000 MeV - 100.000000 TeV
          StepLim=SafetyPlus Rfact=0.080000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

eIoni:  for e-  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.200000, 0.010000 mm), integ: 3, fluct: 1, linLossLim= 0.010000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
       LowEnergyIoni : Emin=0.000000 eV  Emax=100.000000 keV  deltaVI
        MollerBhabha : Emin=100.000000 keV Emax=100.000000 TeV  deltaVI

eBrem:  for e-  XStype:4  SubType=3
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      LPM flag: 1 for E > 1.000000 GeV,  VertexHighEnergyTh(GeV)= 100000.000000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
             eBremSB : Emin=0.000000 eV  Emax=1.000000 GeV  AngularGen2BS
            eBremLPM : Emin=1.000000 GeV Emax=100.000000 TeV  AngularGen2BS

ePairProd:  for e-  XStype:1  SubType=4
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 0
      Sampling table 25x1001; from 0.100000 GeV to 100.000000 TeV 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
           ePairProd : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

CoulombScat:  for e- XStype:1 SubType=1 BuildTable=1
      Lambda table from 100.000000 MeV to 100.000000 TeV, 20 bins/decade, spline: 0
      ThetaMin(p) < Theta(degree) < 180; pLimit(GeV^1)= 0.139531
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  eCoulombScattering : Emin=100.000000 MeV Emax=100.000000 TeV

msc:  for e+  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  GoudsmitSaunderson : Emin=0.000000 eV  Emax=100.000000 MeV Nbins=120 100.000000 eV  - 100.000000 MeV
          StepLim=SafetyPlus Rfact=0.080000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm
        WentzelVIUni : Emin=100.000000 MeV Emax=100.000000 TeV Nbins=120 100.000000 MeV - 100.000000 TeV
          StepLim=SafetyPlus Rfact=0.080000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

eIoni:  for e+  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.200000, 0.010000 mm), integ: 3, fluct: 1, linLossLim= 0.010000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
        MollerBhabha : Emin=0.000000 eV  Emax=100.000000 TeV  deltaVI

eBrem:  for e+  XStype:4  SubType=3
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      LPM flag: 1 for E > 1.000000 GeV,  VertexHighEnergyTh(GeV)= 100000.000000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
             eBremSB : Emin=0.000000 eV  Emax=1.000000 GeV  AngularGen2BS
            eBremLPM : Emin=1.000000 GeV Emax=100.000000 TeV  AngularGen2BS

ePairProd:  for e+  XStype:1  SubType=4
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 0
      Sampling table 25x1001; from 0.100000 GeV to 100.000000 TeV 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
           ePairProd : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

annihil:  for e+ XStype:2 SubType=5 BuildTable=0
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
            eplus2gg : Emin=0.000000 eV  Emax=100.000000 TeV

CoulombScat:  for e+ XStype:1 SubType=1 BuildTable=1
      Lambda table from 100.000000 MeV to 100.000000 TeV, 20 bins/decade, spline: 0
      ThetaMin(p) < Theta(degree) < 180; pLimit(GeV^1)= 0.139531
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  eCoulombScattering : Emin=100.000000 MeV Emax=100.000000 TeV

msc:  for proton  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
        WentzelVIUni : Emin=0.000000 eV  Emax=100.000000 TeV Nbins=240 100.000000 eV  - 100.000000 TeV
          StepLim=Minimal Rfact=0.200000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

hIoni:  for proton  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.100000, 0.050000 mm), integ: 3, fluct: 1, linLossLim= 0.010000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
               Bragg : Emin=0.000000 eV  Emax=2.000000 MeV  deltaVI
          BetheBloch : Emin=2.000000 MeV Emax=100.000000 TeV  deltaVI

hBrems:  for proton  XStype:1  SubType=3
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
               hBrem : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

hPairProd:  for proton  XStype:1  SubType=4
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      Sampling table 17x1001; from 7.506176 GeV to 100.000000 TeV 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
           hPairProd : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

CoulombScat:  for proton XStype:1 SubType=1 BuildTable=1
      Lambda table from threshold  to 100.000000 TeV, 20 bins/decade, spline: 0
      ThetaMin(p) < Theta(degree) < 180; pLimit(GeV^1)= 0.139531
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  eCoulombScattering : Emin=0.000000 eV  Emax=100.000000 TeV

nuclearStopping:  for proton SubType=8 BuildTable=0
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
   ICRU49NucStopping : Emin=0.000000 eV  Emax=1.000000 MeV

msc:  for GenericIon  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
            UrbanMsc : Emin=0.000000 eV  Emax=100.000000 TeV
          StepLim=Minimal Rfact=0.200000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

ionIoni:  for GenericIon  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.100000, 0.001000 mm), integ: 3, fluct: 1, linLossLim= 0.020000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
    LindhardSorensen : Emin=0.000000 eV  Emax=100.000000 TeV  deltaVI

nuclearStopping:  for GenericIon SubType=8 BuildTable=0
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
   ICRU49NucStopping : Emin=0.000000 eV  Emax=1.000000 MeV

msc:  for alpha  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
            UrbanMsc : Emin=0.000000 eV  Emax=100.000000 TeV
          StepLim=Minimal Rfact=0.200000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

ionIoni:  for alpha  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.100000, 0.020000 mm), integ: 3, fluct: 1, linLossLim= 0.020000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
            BraggIon : Emin=0.000000 eV  Emax=7.945199 MeV  deltaVI
          BetheBloch : Emin=7.945199 MeV Emax=100.000000 TeV  deltaVI

nuclearStopping:  for alpha SubType=8 BuildTable=0
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
   ICRU49NucStopping : Emin=0.000000 eV  Emax=1.000000 MeV

msc:  for anti_proton  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
        WentzelVIUni : Emin=0.000000 eV  Emax=100.000000 TeV Nbins=240 100.000000 eV  - 100.000000 TeV
          StepLim=Minimal Rfact=0.200000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

hIoni:  for anti_proton  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.100000, 0.050000 mm), integ: 3, fluct: 1, linLossLim= 0.010000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
            ICRU73QO : Emin=0.000000 eV  Emax=2.000000 MeV  deltaVI
          BetheBloch : Emin=2.000000 MeV Emax=100.000000 TeV  deltaVI

hBrems:  for anti_proton  XStype:1  SubType=3
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
               hBrem : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

hPairProd:  for anti_proton  XStype:1  SubType=4
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      Sampling table 17x1001; from 7.506176 GeV to 100.000000 TeV 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
           hPairProd : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

CoulombScat:  for anti_proton XStype:1 SubType=1 BuildTable=1
      Lambda table from threshold  to 100.000000 TeV, 20 bins/decade, spline: 0
      ThetaMin(p) < Theta(degree) < 180; pLimit(GeV^1)= 0.139531
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  eCoulombScattering : Emin=0.000000 eV  Emax=100.000000 TeV

msc:  for kaon+  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
        WentzelVIUni : Emin=0.000000 eV  Emax=100.000000 TeV Nbins=240 100.000000 eV  - 100.000000 TeV
          StepLim=Minimal Rfact=0.200000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

hIoni:  for kaon+  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.100000, 0.050000 mm), integ: 3, fluct: 1, linLossLim= 0.010000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
               Bragg : Emin=0.000000 eV  Emax=1.052311 MeV  deltaVI
          BetheBloch : Emin=1.052311 MeV Emax=100.000000 TeV  deltaVI

hBrems:  for kaon+  XStype:1  SubType=3
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
               hBrem : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

hPairProd:  for kaon+  XStype:1  SubType=4
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      Sampling table 18x1001; from 3.949416 GeV to 100.000000 TeV 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
           hPairProd : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

CoulombScat:  for kaon+ XStype:1 SubType=1 BuildTable=1
      Lambda table from threshold  to 100.000000 TeV, 20 bins/decade, spline: 0
      ThetaMin(p) < Theta(degree) < 180; pLimit(GeV^1)= 0.139531
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  eCoulombScattering : Emin=0.000000 eV  Emax=100.000000 TeV

msc:  for kaon-  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
        WentzelVIUni : Emin=0.000000 eV  Emax=100.000000 TeV Nbins=240 100.000000 eV  - 100.000000 TeV
          StepLim=Minimal Rfact=0.200000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

hIoni:  for kaon-  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.100000, 0.050000 mm), integ: 3, fluct: 1, linLossLim= 0.010000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
            ICRU73QO : Emin=0.000000 eV  Emax=1.052311 MeV  deltaVI
          BetheBloch : Emin=1.052311 MeV Emax=100.000000 TeV  deltaVI

hBrems:  for kaon-  XStype:1  SubType=3
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
               hBrem : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

hPairProd:  for kaon-  XStype:1  SubType=4
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      Sampling table 18x1001; from 3.949416 GeV to 100.000000 TeV 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
           hPairProd : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

CoulombScat:  for kaon- XStype:1 SubType=1 BuildTable=1
      Used Lambda table of kaon+
      ThetaMin(p) < Theta(degree) < 180; pLimit(GeV^1)= 0.139531
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  eCoulombScattering : Emin=0.000000 eV  Emax=100.000000 TeV

msc:  for mu+  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
        WentzelVIUni : Emin=0.000000 eV  Emax=100.000000 TeV Nbins=240 100.000000 eV  - 100.000000 TeV
          StepLim=Minimal Rfact=0.200000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

muIoni:  for mu+  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.100000, 0.050000 mm), integ: 3, fluct: 1, linLossLim= 0.010000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
               Bragg : Emin=0.000000 eV  Emax=200.000000 keV  deltaVI
        MuBetheBloch : Emin=200.000000 keV Emax=100.000000 TeV  deltaVI

muBrems:  for mu+  XStype:1  SubType=3
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
              MuBrem : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

muPairProd:  for mu+  XStype:1  SubType=4
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      Sampling table 21x1001; from 0.850000 GeV to 100.000000 TeV 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
          muPairProd : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

CoulombScat:  for mu+ XStype:1 SubType=1 BuildTable=1
      Lambda table from threshold  to 100.000000 TeV, 20 bins/decade, spline: 0
      ThetaMin(p) < Theta(degree) < 180; pLimit(GeV^1)= 0.139531
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  eCoulombScattering : Emin=0.000000 eV  Emax=100.000000 TeV

msc:  for mu-  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
        WentzelVIUni : Emin=0.000000 eV  Emax=100.000000 TeV Nbins=240 100.000000 eV  - 100.000000 TeV
          StepLim=Minimal Rfact=0.200000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

muIoni:  for mu-  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.100000, 0.050000 mm), integ: 3, fluct: 1, linLossLim= 0.010000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
            ICRU73QO : Emin=0.000000 eV  Emax=200.000000 keV  deltaVI
        MuBetheBloch : Emin=200.000000 keV Emax=100.000000 TeV  deltaVI

muBrems:  for mu-  XStype:1  SubType=3
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
              MuBrem : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

muPairProd:  for mu-  XStype:1  SubType=4
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      Sampling table 21x1001; from 0.850000 GeV to 100.000000 TeV 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
          muPairProd : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

CoulombScat:  for mu- XStype:1 SubType=1 BuildTable=1
      Used Lambda table of mu+
      ThetaMin(p) < Theta(degree) < 180; pLimit(GeV^1)= 0.139531
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  eCoulombScattering : Emin=0.000000 eV  Emax=100.000000 TeV

msc:  for pi+  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
        WentzelVIUni : Emin=0.000000 eV  Emax=100.000000 TeV Nbins=240 100.000000 eV  - 100.000000 TeV
          StepLim=Minimal Rfact=0.200000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

hIoni:  for pi+  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.100000, 0.050000 mm), integ: 3, fluct: 1, linLossLim= 0.010000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
               Bragg : Emin=0.000000 eV  Emax=297.504557 keV  deltaVI
          BetheBloch : Emin=297.504557 keV Emax=100.000000 TeV  deltaVI

hBrems:  for pi+  XStype:1  SubType=3
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
               hBrem : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

hPairProd:  for pi+  XStype:1  SubType=4
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      Sampling table 20x1001; from 1.116561 GeV to 100.000000 TeV 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
           hPairProd : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

CoulombScat:  for pi+ XStype:1 SubType=1 BuildTable=1
      Lambda table from threshold  to 100.000000 TeV, 20 bins/decade, spline: 0
      ThetaMin(p) < Theta(degree) < 180; pLimit(GeV^1)= 0.139531
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  eCoulombScattering : Emin=0.000000 eV  Emax=100.000000 TeV

msc:  for pi-  SubType= 10
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
        WentzelVIUni : Emin=0.000000 eV  Emax=100.000000 TeV Nbins=240 100.000000 eV  - 100.000000 TeV
          StepLim=Minimal Rfact=0.200000 Gfact=2.500000 Sfact=0.600000 DispFlag:1 Skin=3.000000 Llim=1.000000 mm

hIoni:  for pi-  XStype:3  SubType=2
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      StepFunction=(0.100000, 0.050000 mm), integ: 3, fluct: 1, linLossLim= 0.010000
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
            ICRU73QO : Emin=0.000000 eV  Emax=297.504557 keV  deltaVI
          BetheBloch : Emin=297.504557 keV Emax=100.000000 TeV  deltaVI

hBrems:  for pi-  XStype:1  SubType=3
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
               hBrem : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

hPairProd:  for pi-  XStype:1  SubType=4
      dE/dx and range tables from 100.000000 eV  to 100.000000 TeV in 240 bins
      Lambda tables from threshold to 100.000000 TeV, 20 bins/decade, spline: 1
      Sampling table 20x1001; from 1.116561 GeV to 100.000000 TeV 
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
           hPairProd : Emin=0.000000 eV  Emax=100.000000 TeV  ModifiedMephi

CoulombScat:  for pi- XStype:1 SubType=1 BuildTable=1
      Used Lambda table of pi+
      ThetaMin(p) < Theta(degree) < 180; pLimit(GeV^1)= 0.139531
      ===== EM models for the G4Region  DefaultRegionForTheWorld ======
  eCoulombScattering : Emin=0.000000 eV  Emax=100.000000 TeV
Start closing geometry.
G4GeometryManager::ReportVoxelStats -- Voxel Statistics

    Total memory consumed for geometry optimisation:   3289151 kByte
    Total CPU time elapsed for geometry optimisation: 45.71 seconds

    Voxelisation: top CPU users:
    Percent   Total CPU    System CPU       Memory  Volume
    -------   ----------   ----------     --------  ----------
      94.47        43.18         3.45      3289151k phantomLogical

    Voxelisation: top memory users:
    Percent     Memory      Heads    Nodes   Pointers    Total CPU    Volume
    -------   --------     ------   ------   --------   ----------    ----------
     100.00    3289151k    605813 46718609  135247371        43.18    phantomLogical

==============================================================================
--> G4TaskRunManager::ComputeNumberOfTasks() --> 0 tasks with 1 events/task...
==============================================================================


=========================================================================
--> G4TaskRunManager::CreateAndStartWorkers() --> Initializing workers...
=========================================================================


Region <DefaultRegionForTheWorld> --  -- appears in <worldPhysical> world volume
 This region is in the mass world.
 Root logical volume(s) : worldLogical 
 Pointers : G4VUserRegionInformation[0x0], G4UserLimits[0x0], G4FastSimulationManager[0x0], G4UserSteppingAction[0x0]
 Materials : G4_Galactic Skin_legs_insensitive RST_trunk RST_arms Skin_legs_sensitive(50-100) RST_legs BB(60-70) Small_intestine_wall(130-150) ET2(0-40) ET1(0-8) Muscle_arms Air_inside_body Cranium_cortical Cranium_spongiosa RST_head Oral_mucosa_lips_and_cheeks Skin_arms_sensitive(50-100) Wrists_and_hand_bones_cortical Skin_arms_insensitive Heart_wall Femora_lower_spongiosa Tongue_lower Small_intestine_wall(0-130) ET2(65-Surface) Skin_trunk_insensitive Muscle_trunk Oral_mucosa_tongue Teeth Small_intestine_contents(-500-0) Ribs_cortical Blood_in_large_arteries_trunk Thymus Stomach_wall(60-100) Lung(AI)_left Skin_head_insensitive Muscle_legs Teeth_retention_region Blood_in_large_veins_trunk Lung(AI)_right Ureter_right Testis_left Skin_trunk_sensitive(50-100) Cervical_spine_spongiosa Muscle_head Thyroid Small_intestine_wall(200-surface) Blood_in_large_arteries_legs Trachea Thoracic_spine_spongiosa Blood_in_large_veins_arms ET2(55-65) Tibiae_cortical Blood_in_large_veins_legs BB(40-50) Small_intestine_contents(centre--500) Urinary_bladder_wall_insensitive Ulnae_and_radii_medullary_cavity Urinary_bladder_content Ankles_and_foot_cortical Ankles_and_foot_spongiosa Ribs_spongiosa Ulnae_and_radii_cotical Testis_right BB(-6-0) Scapulae_cortical Blood_in_large_arteries_arms Pelvis_cortical Lumbar_spine_cortical Urinary_bladder_wall_sensitive(75-118) Small_intestine_wall(150-200) ET1(8-40) Spinal_cord Lymphatic_nodes_arms Thoracic_spine_cortical Oesophagus_wall(200-surface) Adrenal_left Skin_head_sensitive(50-100) Oesophagus_wall(190-200) Prostate Blood_in_heart_chamber Cervical_spine_cortical Lymphatic_nodes_trunk Sacrum_spongiosa Adrenal_right Transverse_colon_wall_right(0-280) BB(0-10) Brain BB(50-60) Tonsils Sigmoid_colon_wall(300-surface) Scapulae_spongiosa ET2(50-55) Transverse_colon_wall_left(280-300) ET1(50-Surface) ET2(-15-0) Eye_lens_insensitive_right Lymphatic_nodes_ET Femora_lower_cortical BB(-11--6) Tongue_upper(food) Cartilage_discs Clavicles_cortical BB(35-40) Ureter_left Oesophagus_contents BB(70-surface) Salivary_glands_left Blood_in_large_arteries_head Mandible_cortical ET2(40-50) Stomach_wall(100-300) Mandible_spongiosa Lumbar_spine_spongiosa ET1(40-50) Kidney_right_cortex Sacrum_cortical Pelvis_spongiosa Sternum_cortical Humeri_lower_spongiosa Descending_colon_wall(280-300) Breast_right_adipose_tissue Wrists_and_hand_bones_spongiosa Femora_upper_cortical Stomach_wall(0-60) Tibiae_medullary_cavity Transverse_colon_contents_right Ascending_colon_wall(300-surface) Ascending_colon_wall(0-280) Sigmoid_colon_wall(0-280) Tibiae_spongiosa Femora_upper_spongiosa Kidney_left_cortex Gall_bladder_wall Ulnae_and_radii_spongiosa Humeri_upper_cortical Kidney_right_medulla Oesophagus_wall(0-190) Transverse_colon_wall_left(0-280) Cartilage_costal Liver Salivary_glandss_right BB(10-35) Cornea_left Oral_mucosa_moth_floor Transverse_colon_wall_right(300-surface) Descending_colon_wall(0-280) Ascending_colon_content Blood_in_large_veins_head Vitreous_right Cornea_right Pituitary_gland Sigmoid_colon_contents Lymphatic_nodes_legs Lymphatic_nodes_head Clavicles_spongiosa Breast_right_glandular_tissue Sternum_spongiosa Lymphatic_nodes_thoracic Stomach_wall(300-surface) Ascending_colon_wall(280-300) Sigmoid_colon_wall(280-300) Humeri_upper_medullary_cavity Eye_lens_insensitive_left Vitreous_left Spleen Humeri_lower_cortical Aqueous_left Stomach_contents Aqueous_right Transverse_colon_wall_right(280-300) Breast_left_adipose_tissue Rectum_wall Descending_colon_wall(300-surface) Kidney_left_medulla Femora_lower_medullary_cavity Breast_left_glandular_tissue Eye_lens_sensitive_right Femora_upper_medullary_cavity Transverse_colon_wall_left(300-surface) Humeri_lower_medullary_cavity Gall_bladder_contents Transverse_colon_content_left Kidney_right_pelvis Eye_lens_sensitive_left Humeri_upper_spogiosa Pancreas Descending_colon_content Kidney_left_pelvis 
 Production cuts :   gamma 1 mm      e- 1 mm      e+ 1 mm  proton 1 mm 

Region <DefaultRegionForParallelWorld> --  -- is not associated to any world.
 Root logical volume(s) : 
 Pointers : G4VUserRegionInformation[0x0], G4UserLimits[0x0], G4FastSimulationManager[0x0], G4UserSteppingAction[0x0]
 Materials : 
 Production cuts :   gamma 1 mm      e- 1 mm      e+ 1 mm  proton 1 mm 

========= Table of registered couples ============================

Index : 0     used in the geometry : Yes
 Material : G4_Galactic
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  990 eV     e-  990 eV     e+  990 eV  proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 1     used in the geometry : Yes
 Material : Skin_legs_insensitive
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93499 keV    e-  369.615 keV    e+  359.18 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 2     used in the geometry : Yes
 Material : RST_trunk
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.54133 keV    e-  339.86 keV    e+  330.771 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 3     used in the geometry : Yes
 Material : RST_arms
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.54133 keV    e-  339.86 keV    e+  330.771 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 4     used in the geometry : Yes
 Material : Skin_legs_sensitive(50-100)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93499 keV    e-  369.615 keV    e+  359.18 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 5     used in the geometry : Yes
 Material : RST_legs
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.54133 keV    e-  339.86 keV    e+  330.771 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 6     used in the geometry : Yes
 Material : BB(60-70)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 7     used in the geometry : Yes
 Material : Small_intestine_wall(130-150)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 8     used in the geometry : Yes
 Material : ET2(0-40)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 9     used in the geometry : Yes
 Material : ET1(0-8)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 10     used in the geometry : Yes
 Material : Muscle_arms
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93667 keV    e-  360.8 keV    e+  350.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 11     used in the geometry : Yes
 Material : Air_inside_body
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  990 eV     e-  990 eV     e+  990 eV  proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 12     used in the geometry : Yes
 Material : Cranium_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 13     used in the geometry : Yes
 Material : Cranium_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.22162 keV    e-  382.315 keV    e+  371.247 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 14     used in the geometry : Yes
 Material : RST_head
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.54133 keV    e-  339.86 keV    e+  330.771 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 15     used in the geometry : Yes
 Material : Oral_mucosa_lips_and_cheeks
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93667 keV    e-  360.8 keV    e+  350.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 16     used in the geometry : Yes
 Material : Skin_arms_sensitive(50-100)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93499 keV    e-  369.615 keV    e+  359.18 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 17     used in the geometry : Yes
 Material : Wrists_and_hand_bones_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 18     used in the geometry : Yes
 Material : Skin_arms_insensitive
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93499 keV    e-  369.615 keV    e+  359.18 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 19     used in the geometry : Yes
 Material : Heart_wall
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.9379 keV    e-  361.592 keV    e+  351.442 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 20     used in the geometry : Yes
 Material : Femora_lower_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.00206 keV    e-  373.473 keV    e+  362.881 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 21     used in the geometry : Yes
 Material : Tongue_lower
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93667 keV    e-  360.8 keV    e+  350.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 22     used in the geometry : Yes
 Material : Small_intestine_wall(0-130)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 23     used in the geometry : Yes
 Material : ET2(65-Surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 24     used in the geometry : Yes
 Material : Skin_trunk_insensitive
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93499 keV    e-  369.615 keV    e+  359.18 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 25     used in the geometry : Yes
 Material : Muscle_trunk
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93667 keV    e-  360.8 keV    e+  350.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 26     used in the geometry : Yes
 Material : Oral_mucosa_tongue
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93667 keV    e-  360.8 keV    e+  350.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 27     used in the geometry : Yes
 Material : Teeth
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  6.71752 keV    e-  634.21 keV    e+  609.635 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 28     used in the geometry : Yes
 Material : Small_intestine_contents(-500-0)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88376 keV    e-  358.325 keV    e+  348.336 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 29     used in the geometry : Yes
 Material : Ribs_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 30     used in the geometry : Yes
 Material : Blood_in_large_arteries_trunk
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96963 keV    e-  362.885 keV    e+  352.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 31     used in the geometry : Yes
 Material : Thymus
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 32     used in the geometry : Yes
 Material : Stomach_wall(60-100)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 33     used in the geometry : Yes
 Material : Lung(AI)_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.00418 keV    e-  197.94 keV    e+  193.838 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 34     used in the geometry : Yes
 Material : Skin_head_insensitive
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93499 keV    e-  369.615 keV    e+  359.18 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 35     used in the geometry : Yes
 Material : Muscle_legs
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93667 keV    e-  360.8 keV    e+  350.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 36     used in the geometry : Yes
 Material : Teeth_retention_region
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88376 keV    e-  358.325 keV    e+  348.336 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 37     used in the geometry : Yes
 Material : Blood_in_large_veins_trunk
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96963 keV    e-  362.885 keV    e+  352.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 38     used in the geometry : Yes
 Material : Lung(AI)_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.00418 keV    e-  197.94 keV    e+  193.838 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 39     used in the geometry : Yes
 Material : Ureter_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 40     used in the geometry : Yes
 Material : Testis_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.94463 keV    e-  359.64 keV    e+  349.574 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 41     used in the geometry : Yes
 Material : Skin_trunk_sensitive(50-100)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93499 keV    e-  369.615 keV    e+  359.18 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 42     used in the geometry : Yes
 Material : Cervical_spine_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.80594 keV    e-  362.335 keV    e+  352.224 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 43     used in the geometry : Yes
 Material : Muscle_head
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93667 keV    e-  360.8 keV    e+  350.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 44     used in the geometry : Yes
 Material : Thyroid
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.95629 keV    e-  361.469 keV    e+  351.319 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 45     used in the geometry : Yes
 Material : Small_intestine_wall(200-surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 46     used in the geometry : Yes
 Material : Blood_in_large_arteries_legs
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96963 keV    e-  362.885 keV    e+  352.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 47     used in the geometry : Yes
 Material : Trachea
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 48     used in the geometry : Yes
 Material : Thoracic_spine_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88783 keV    e-  365.906 keV    e+  355.639 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 49     used in the geometry : Yes
 Material : Blood_in_large_veins_arms
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96963 keV    e-  362.885 keV    e+  352.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 50     used in the geometry : Yes
 Material : ET2(55-65)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 51     used in the geometry : Yes
 Material : Tibiae_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 52     used in the geometry : Yes
 Material : Blood_in_large_veins_legs
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96963 keV    e-  362.885 keV    e+  352.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 53     used in the geometry : Yes
 Material : BB(40-50)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 54     used in the geometry : Yes
 Material : Small_intestine_contents(centre--500)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88376 keV    e-  358.325 keV    e+  348.336 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 55     used in the geometry : Yes
 Material : Urinary_bladder_wall_insensitive
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.95353 keV    e-  359.023 keV    e+  348.982 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 56     used in the geometry : Yes
 Material : Ulnae_and_radii_medullary_cavity
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.49917 keV    e-  351.594 keV    e+  342.076 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 57     used in the geometry : Yes
 Material : Urinary_bladder_content
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.99933 keV    e-  359.093 keV    e+  349.032 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 58     used in the geometry : Yes
 Material : Ankles_and_foot_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 59     used in the geometry : Yes
 Material : Ankles_and_foot_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.00206 keV    e-  373.473 keV    e+  362.881 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 60     used in the geometry : Yes
 Material : Ribs_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.24433 keV    e-  383.192 keV    e+  372.08 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 61     used in the geometry : Yes
 Material : Ulnae_and_radii_cotical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 62     used in the geometry : Yes
 Material : Testis_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.94463 keV    e-  359.64 keV    e+  349.574 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 63     used in the geometry : Yes
 Material : BB(-6-0)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 64     used in the geometry : Yes
 Material : Scapulae_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 65     used in the geometry : Yes
 Material : Blood_in_large_arteries_arms
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96963 keV    e-  362.885 keV    e+  352.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 66     used in the geometry : Yes
 Material : Pelvis_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 67     used in the geometry : Yes
 Material : Lumbar_spine_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 68     used in the geometry : Yes
 Material : Urinary_bladder_wall_sensitive(75-118)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.95353 keV    e-  359.023 keV    e+  348.982 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 69     used in the geometry : Yes
 Material : Small_intestine_wall(150-200)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 70     used in the geometry : Yes
 Material : ET1(8-40)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 71     used in the geometry : Yes
 Material : Spinal_cord
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 72     used in the geometry : Yes
 Material : Lymphatic_nodes_arms
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96377 keV    e-  357.772 keV    e+  347.793 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 73     used in the geometry : Yes
 Material : Thoracic_spine_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 74     used in the geometry : Yes
 Material : Oesophagus_wall(200-surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.86422 keV    e-  358.912 keV    e+  348.911 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 75     used in the geometry : Yes
 Material : Adrenal_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.85992 keV    e-  358.71 keV    e+  348.721 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 76     used in the geometry : Yes
 Material : Skin_head_sensitive(50-100)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93499 keV    e-  369.615 keV    e+  359.18 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 77     used in the geometry : Yes
 Material : Oesophagus_wall(190-200)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.86422 keV    e-  358.912 keV    e+  348.911 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 78     used in the geometry : Yes
 Material : Prostate
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 79     used in the geometry : Yes
 Material : Blood_in_heart_chamber
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96963 keV    e-  362.885 keV    e+  352.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 80     used in the geometry : Yes
 Material : Cervical_spine_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 81     used in the geometry : Yes
 Material : Lymphatic_nodes_trunk
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96377 keV    e-  357.772 keV    e+  347.793 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 82     used in the geometry : Yes
 Material : Sacrum_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.74882 keV    e-  359.456 keV    e+  349.484 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 83     used in the geometry : Yes
 Material : Adrenal_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.85992 keV    e-  358.71 keV    e+  348.721 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 84     used in the geometry : Yes
 Material : Transverse_colon_wall_right(0-280)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 85     used in the geometry : Yes
 Material : BB(0-10)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 86     used in the geometry : Yes
 Material : Brain
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92639 keV    e-  360.127 keV    e+  350.05 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 87     used in the geometry : Yes
 Material : BB(50-60)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 88     used in the geometry : Yes
 Material : Tonsils
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 89     used in the geometry : Yes
 Material : Sigmoid_colon_wall(300-surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 90     used in the geometry : Yes
 Material : Scapulae_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.34215 keV    e-  388.432 keV    e+  377.091 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 91     used in the geometry : Yes
 Material : ET2(50-55)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 92     used in the geometry : Yes
 Material : Transverse_colon_wall_left(280-300)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 93     used in the geometry : Yes
 Material : ET1(50-Surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 94     used in the geometry : Yes
 Material : ET2(-15-0)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93533 keV    e-  351.27 keV    e+  341.576 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 95     used in the geometry : Yes
 Material : Eye_lens_insensitive_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.89856 keV    e-  361.815 keV    e+  351.651 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 96     used in the geometry : Yes
 Material : Lymphatic_nodes_ET
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96377 keV    e-  357.772 keV    e+  347.793 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 97     used in the geometry : Yes
 Material : Femora_lower_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 98     used in the geometry : Yes
 Material : BB(-11--6)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93533 keV    e-  351.27 keV    e+  341.576 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 99     used in the geometry : Yes
 Material : Tongue_upper(food)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93667 keV    e-  360.8 keV    e+  350.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 100     used in the geometry : Yes
 Material : Cartilage_discs
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.08934 keV    e-  369.292 keV    e+  358.8 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 101     used in the geometry : Yes
 Material : Clavicles_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 102     used in the geometry : Yes
 Material : BB(35-40)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 103     used in the geometry : Yes
 Material : Ureter_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 104     used in the geometry : Yes
 Material : Oesophagus_contents
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88376 keV    e-  358.325 keV    e+  348.336 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 105     used in the geometry : Yes
 Material : BB(70-surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 106     used in the geometry : Yes
 Material : Salivary_glands_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 107     used in the geometry : Yes
 Material : Blood_in_large_arteries_head
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96963 keV    e-  362.885 keV    e+  352.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 108     used in the geometry : Yes
 Material : Mandible_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 109     used in the geometry : Yes
 Material : ET2(40-50)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 110     used in the geometry : Yes
 Material : Stomach_wall(100-300)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 111     used in the geometry : Yes
 Material : Mandible_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.56826 keV    e-  400.299 keV    e+  388.334 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 112     used in the geometry : Yes
 Material : Lumbar_spine_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.02744 keV    e-  372.443 keV    e+  361.878 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 113     used in the geometry : Yes
 Material : ET1(40-50)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 114     used in the geometry : Yes
 Material : Kidney_right_cortex
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.94836 keV    e-  361.684 keV    e+  351.524 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 115     used in the geometry : Yes
 Material : Sacrum_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 116     used in the geometry : Yes
 Material : Pelvis_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.06368 keV    e-  375.042 keV    e+  364.338 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 117     used in the geometry : Yes
 Material : Sternum_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 118     used in the geometry : Yes
 Material : Humeri_lower_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.00206 keV    e-  373.473 keV    e+  362.881 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 119     used in the geometry : Yes
 Material : Descending_colon_wall(280-300)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 120     used in the geometry : Yes
 Material : Breast_right_adipose_tissue
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.50816 keV    e-  344.246 keV    e+  334.965 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 121     used in the geometry : Yes
 Material : Wrists_and_hand_bones_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.00206 keV    e-  373.473 keV    e+  362.881 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 122     used in the geometry : Yes
 Material : Femora_upper_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 123     used in the geometry : Yes
 Material : Stomach_wall(0-60)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 124     used in the geometry : Yes
 Material : Tibiae_medullary_cavity
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.49917 keV    e-  351.594 keV    e+  342.076 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 125     used in the geometry : Yes
 Material : Transverse_colon_contents_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88376 keV    e-  358.325 keV    e+  348.336 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 126     used in the geometry : Yes
 Material : Ascending_colon_wall(300-surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 127     used in the geometry : Yes
 Material : Ascending_colon_wall(0-280)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 128     used in the geometry : Yes
 Material : Sigmoid_colon_wall(0-280)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 129     used in the geometry : Yes
 Material : Tibiae_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.00206 keV    e-  373.473 keV    e+  362.881 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 130     used in the geometry : Yes
 Material : Femora_upper_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.0741 keV    e-  375.684 keV    e+  364.947 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 131     used in the geometry : Yes
 Material : Kidney_left_cortex
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.94836 keV    e-  361.684 keV    e+  351.524 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 132     used in the geometry : Yes
 Material : Gall_bladder_wall
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 133     used in the geometry : Yes
 Material : Ulnae_and_radii_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.00206 keV    e-  373.473 keV    e+  362.881 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 134     used in the geometry : Yes
 Material : Humeri_upper_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 135     used in the geometry : Yes
 Material : Kidney_right_medulla
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.94836 keV    e-  361.684 keV    e+  351.524 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 136     used in the geometry : Yes
 Material : Oesophagus_wall(0-190)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.86422 keV    e-  358.912 keV    e+  348.911 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 137     used in the geometry : Yes
 Material : Transverse_colon_wall_left(0-280)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 138     used in the geometry : Yes
 Material : Cartilage_costal
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.08934 keV    e-  369.292 keV    e+  358.8 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 139     used in the geometry : Yes
 Material : Liver
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.95629 keV    e-  363.011 keV    e+  352.8 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 140     used in the geometry : Yes
 Material : Salivary_glandss_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 141     used in the geometry : Yes
 Material : BB(10-35)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 142     used in the geometry : Yes
 Material : Cornea_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.98825 keV    e-  372.037 keV    e+  361.486 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 143     used in the geometry : Yes
 Material : Oral_mucosa_moth_floor
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.93667 keV    e-  360.8 keV    e+  350.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 144     used in the geometry : Yes
 Material : Transverse_colon_wall_right(300-surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 145     used in the geometry : Yes
 Material : Descending_colon_wall(0-280)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 146     used in the geometry : Yes
 Material : Ascending_colon_content
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88376 keV    e-  358.325 keV    e+  348.336 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 147     used in the geometry : Yes
 Material : Blood_in_large_veins_head
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96963 keV    e-  362.885 keV    e+  352.673 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 148     used in the geometry : Yes
 Material : Vitreous_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.97103 keV    e-  358.611 keV    e+  348.595 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 149     used in the geometry : Yes
 Material : Cornea_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.98825 keV    e-  372.037 keV    e+  361.486 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 150     used in the geometry : Yes
 Material : Pituitary_gland
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.8375 keV    e-  357.986 keV    e+  348.043 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 151     used in the geometry : Yes
 Material : Sigmoid_colon_contents
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88376 keV    e-  358.325 keV    e+  348.336 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 152     used in the geometry : Yes
 Material : Lymphatic_nodes_legs
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96377 keV    e-  357.772 keV    e+  347.793 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 153     used in the geometry : Yes
 Material : Lymphatic_nodes_head
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96377 keV    e-  357.772 keV    e+  347.793 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 154     used in the geometry : Yes
 Material : Clavicles_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.18965 keV    e-  381.04 keV    e+  370.032 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 155     used in the geometry : Yes
 Material : Breast_right_glandular_tissue
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.76089 keV    e-  356.529 keV    e+  346.689 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 156     used in the geometry : Yes
 Material : Sternum_spongiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.77716 keV    e-  360.916 keV    e+  350.864 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 157     used in the geometry : Yes
 Material : Lymphatic_nodes_thoracic
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96377 keV    e-  357.772 keV    e+  347.793 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 158     used in the geometry : Yes
 Material : Stomach_wall(300-surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 159     used in the geometry : Yes
 Material : Ascending_colon_wall(280-300)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 160     used in the geometry : Yes
 Material : Sigmoid_colon_wall(280-300)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 161     used in the geometry : Yes
 Material : Humeri_upper_medullary_cavity
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.49917 keV    e-  351.594 keV    e+  342.076 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 162     used in the geometry : Yes
 Material : Eye_lens_insensitive_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.89856 keV    e-  361.815 keV    e+  351.651 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 163     used in the geometry : Yes
 Material : Vitreous_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.97103 keV    e-  358.611 keV    e+  348.595 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 164     used in the geometry : Yes
 Material : Spleen
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.9643 keV    e-  363.22 keV    e+  353.001 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 165     used in the geometry : Yes
 Material : Humeri_lower_cortical
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  5.21027 keV    e-  503.928 keV    e+  486.537 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 166     used in the geometry : Yes
 Material : Aqueous_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96357 keV    e-  357.196 keV    e+  347.251 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 167     used in the geometry : Yes
 Material : Stomach_contents
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88376 keV    e-  358.325 keV    e+  348.336 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 168     used in the geometry : Yes
 Material : Aqueous_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.96357 keV    e-  357.196 keV    e+  347.251 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 169     used in the geometry : Yes
 Material : Transverse_colon_wall_right(280-300)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 170     used in the geometry : Yes
 Material : Breast_left_adipose_tissue
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.50816 keV    e-  344.246 keV    e+  334.965 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 171     used in the geometry : Yes
 Material : Rectum_wall
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 172     used in the geometry : Yes
 Material : Descending_colon_wall(300-surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 173     used in the geometry : Yes
 Material : Kidney_left_medulla
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.94836 keV    e-  361.684 keV    e+  351.524 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 174     used in the geometry : Yes
 Material : Femora_lower_medullary_cavity
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.49917 keV    e-  351.594 keV    e+  342.076 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 175     used in the geometry : Yes
 Material : Breast_left_glandular_tissue
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.76089 keV    e-  356.529 keV    e+  346.689 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 176     used in the geometry : Yes
 Material : Eye_lens_sensitive_right
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.89856 keV    e-  361.815 keV    e+  351.651 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 177     used in the geometry : Yes
 Material : Femora_upper_medullary_cavity
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.49917 keV    e-  351.594 keV    e+  342.076 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 178     used in the geometry : Yes
 Material : Transverse_colon_wall_left(300-surface)
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.92336 keV    e-  358.567 keV    e+  348.559 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 179     used in the geometry : Yes
 Material : Humeri_lower_medullary_cavity
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.49917 keV    e-  351.594 keV    e+  342.076 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 180     used in the geometry : Yes
 Material : Gall_bladder_contents
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.83297 keV    e-  357.783 keV    e+  347.852 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 181     used in the geometry : Yes
 Material : Transverse_colon_content_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88376 keV    e-  358.325 keV    e+  348.336 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 182     used in the geometry : Yes
 Material : Kidney_right_pelvis
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.94836 keV    e-  361.684 keV    e+  351.524 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 183     used in the geometry : Yes
 Material : Eye_lens_sensitive_left
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.89856 keV    e-  361.815 keV    e+  351.651 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 184     used in the geometry : Yes
 Material : Humeri_upper_spogiosa
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  3.45048 keV    e-  393.96 keV    e+  382.307 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 185     used in the geometry : Yes
 Material : Pancreas
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.90971 keV    e-  360.434 keV    e+  350.342 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 186     used in the geometry : Yes
 Material : Descending_colon_content
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.88376 keV    e-  358.325 keV    e+  348.336 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

Index : 187     used in the geometry : Yes
 Material : Kidney_left_pelvis
 Range cuts        :  gamma  1 mm     e-  1 mm     e+  1 mm  proton 1 mm 
 Energy thresholds :  gamma  2.94836 keV    e-  361.684 keV    e+  351.524 keV proton 100 keV
 Region(s) which use this couple : 
    DefaultRegionForTheWorld

==================================================================

G4VisManager: Using G4TrajectoryDrawByCharge as fallback trajectory model.
See commands in /vis/modeling/trajectories/ for other options.
### Run 0 starts.

========================================================================================
--> G4TaskRunManager::CreateAndStartWorkers() --> Creating 3 tasks with 3 events/task...
========================================================================================

Adding task 0 to task-group...
Adding task 1 to task-group...
Adding task 2 to task-group...
Adding task 3 to task-group...
G4WT0 > ### Run 0 starts on worker thread 0.
G4WT0 > --> Event 0 starts with initial seeds (13049039,61775110).
G4WT0 > 
G4WT0 > *********************************************************************************************************
G4WT0 > * G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
G4WT0 > *********************************************************************************************************
G4WT0 > 
G4WT0 > Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
G4WT0 >     0     2.01     36.5      440      0.01        0        0         0 wholePhantom initStep
G4WT0 >     1     2.01     36.5      440   0.00917  0.00083 0.000184  0.000184 wholePhantom msc
G4WT0 >     2     2.01     36.5      440         0  0.00917  0.00197   0.00215 wholePhantom eIoni
G4WT0 > --> Event 1 starts with initial seeds (99594690,49590222).
G4WT0 > 
G4WT0 > *********************************************************************************************************
G4WT0 > * G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
G4WT0 > *********************************************************************************************************
G4WT0 > 
G4WT0 > Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
G4WT0 >     0     5.01     21.3      396      0.01        0        0         0 wholePhantom initStep
G4WT0 >     1     5.01     21.3      396   0.00957 0.000435 0.000203  0.000203 wholePhantom msc
G4WT0 >     2     5.01     21.3      396         0  0.00957  0.00235   0.00255 wholePhantom eIoni
G4WT0 > --> Event 2 starts with initial seeds (11291749,28987050).
G4WT0 > 
G4WT0 > *********************************************************************************************************
G4WT0 > * G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
G4WT0 > *********************************************************************************************************
G4WT0 > 
G4WT0 > Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
G4WT0 >     0    -3.84     30.6      437      0.01        0        0         0 wholePhantom initStep
G4WT0 >     1    -3.84     30.6      437   0.00917 0.000832 0.000203  0.000203 wholePhantom msc
G4WT0 >     2    -3.84     30.6      437         0  0.00917  0.00218   0.00238 wholePhantom eIoni
G4WT0 > --> Event 3 starts with initial seeds (47304437,83761946).
G4WT0 > 
G4WT0 > *********************************************************************************************************
G4WT0 > * G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
G4WT0 > *********************************************************************************************************
G4WT0 > 
G4WT0 > Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
G4WT0 >     0   -0.293       14      403      0.01        0        0         0 wholePhantom initStep
G4WT0 >     1   -0.293       14      403   0.00903 0.000967 0.000203  0.000203 wholePhantom msc
G4WT0 >     2   -0.293       14      403         0  0.00903  0.00212   0.00232 wholePhantom eIoni
G4WT0 > --> Event 4 starts with initial seeds (35935592,92693821).
G4WT0 > 
G4WT0 > *********************************************************************************************************
G4WT0 > * G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
G4WT0 > *********************************************************************************************************
G4WT0 > 
G4WT0 > Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
G4WT0 >     0      8.9     1.22      598      0.01        0        0         0 wholePhantom initStep
G4WT0 >     1      8.9     1.22      598   0.00972 0.000277 0.000121  0.000121 wholePhantom msc
G4WT0 >     2      8.9     1.22      598         0  0.00972  0.00144   0.00156 wholePhantom eIoni
G4WT0 > --> Event 5 starts with initial seeds (54798860,13721668).
G4WT0 > 
G4WT0 > *********************************************************************************************************
G4WT0 > * G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
G4WT0 > *********************************************************************************************************
G4WT0 > 
G4WT0 > Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
G4WT0 >     0     4.38       19      414      0.01        0        0         0 wholePhantom initStep
G4WT0 >     1     4.38       19      414   0.00971 0.000293  0.00019   0.00019 wholePhantom msc
G4WT0 >     2     4.38       19      414         0  0.00971  0.00226   0.00245 wholePhantom eIoni
G4WT0 > --> Event 6 starts with initial seeds (11900984,65585301).
G4WT0 > 
G4WT0 > *********************************************************************************************************
G4WT0 > * G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
G4WT0 > *********************************************************************************************************
G4WT0 > 
G4WT0 > Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
G4WT0 >     0    -2.22     18.2      400      0.01        0        0         0 wholePhantom initStep
G4WT0 >     1    -2.22     18.2      400   0.00876  0.00124 0.000482  0.000482 wholePhantom msc
G4WT0 >     2    -2.22     18.2      400   0.00851  0.00025 0.000381  0.000863 wholePhantom msc
G4WT0 >     3    -2.22     18.2      400    0.0058   0.0027 0.000362   0.00122 wholePhantom msc
G4WT0 >     4    -2.22     18.2      400         0   0.0058  0.00231   0.00354 wholePhantom eIoni
G4WT0 > --> Event 7 starts with initial seeds (7212557,98037391).
G4WT0 > 
G4WT0 > *********************************************************************************************************
G4WT0 > * G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
G4WT0 > *********************************************************************************************************
G4WT0 > 
G4WT0 > Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
G4WT0 >     0     10.9     5.31      395      0.01        0        0         0 wholePhantom initStep
G4WT0 >     1     10.9     5.31      395   0.00952 0.000485 0.000191  0.000191 wholePhantom msc
G4WT0 >     2     10.9     5.31      395         0  0.00952  0.00219   0.00238 wholePhantom eIoni
G4WT0 > --> Event 8 starts with initial seeds (74052786,96299938).
G4WT0 > 
G4WT0 > *********************************************************************************************************
G4WT0 > * G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
G4WT0 > *********************************************************************************************************
G4WT0 > 
G4WT0 > Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
G4WT0 >     0     8.04     3.77      603      0.01        0        0         0 wholePhantom initStep
G4WT0 >     1     8.04     3.77      603   0.00897  0.00103 0.000203  0.000203 wholePhantom msc
G4WT0 >     2     8.04     3.77      603         0  0.00897   0.0021    0.0023 wholePhantom eIoni
G4WT0 > --> Event 9 starts with initial seeds (685966,37053858).
G4WT0 > 
G4WT0 > *********************************************************************************************************
G4WT0 > * G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
G4WT0 > *********************************************************************************************************
G4WT0 > 
G4WT0 > Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
G4WT0 >     0     15.9    -1.07      604      0.01        0        0         0 wholePhantom initStep
G4WT0 >     1     15.9    -1.07      604   0.00954 0.000465 0.000184  0.000184 wholePhantom msc
G4WT0 >     2     15.9    -1.07      604         0  0.00954  0.00212    0.0023 wholePhantom eIoni
G4WT0 > [thread 0] Thread-local run terminated.
G4WT0 > [thread 0] Run Summary
G4WT0 > [thread 0]   Number of events processed : 10
G4WT0 > [thread 0]   User=0.000000s Real=0.000013s Sys=0.000000s [Cpu=0.0%]
 Run terminated.
Run Summary
  Number of events processed : 10
  User=0.030000s Real=0.295884s Sys=0.040000s [Cpu=23.7%]

=======================================================================
 Run #0 / Number of event processed : 10
=======================================================================
 Init time: 85.9846 s / Run time: 0.296007 s
=======================================================================
                 organ ID|  Organ Mass (g)     SAF (kg-1) Relative Error
                 RBM(DRF)|           0.000      0.000e+00            nan
                  BS(DRF)|           0.000      0.000e+00            nan
                      RBM|           0.000      4.674e-02          0.949
                       BS|           0.000      1.431e-02          0.949
                    Colon|         493.029      0.000e+00            nan
            Colon(target)|           1.244      0.000e+00            nan
                    Lungs|        1198.738      8.342e-02          0.949
                  Stomach|         194.271      0.000e+00            nan
          Stomach(target)|           1.193      0.000e+00            nan
                  Breasts|          25.898      0.000e+00            nan
                   Testes|          37.234      0.000e+00            nan
                       UB|          51.099      0.000e+00            nan
               UB(target)|           1.318      0.000e+00            nan
               Oesophagus|          51.806      1.930e+00          0.949
       Oesophagus(target)|           0.103      0.000e+00            nan
                    Liver|        2360.000      0.000e+00            nan
                  Thyroid|          23.351      0.000e+00            nan
                    Brain|        1517.390      0.000e+00            nan
           SalivaryGlands|          88.091      0.000e+00            nan
                     Skin|        3469.569      0.000e+00            nan
             Skin(target)|         103.521      0.000e+00            nan
                  Adrenal|          17.365      0.000e+00            nan
                      ET1|          11.431      0.000e+00            nan
                      ET2|          29.584      0.000e+00            nan
              ET1(target)|           0.028      0.000e+00            nan
              ET2(target)|           0.098      0.000e+00            nan
              GallBladder|          10.364      0.000e+00            nan
                    Heart|         385.839      0.000e+00            nan
                  Kidneys|         422.146      0.000e+00            nan
               LymphNodes|         189.651      0.000e+00            nan
                   Muscle|       29776.621      0.000e+00            nan
               OralMucosa|           0.133      0.000e+00            nan
                 Pancreas|         173.631      0.000e+00            nan
                 Prostate|          17.618      0.000e+00            nan
                       SI|         862.600      0.000e+00            nan
               SI(target)|           2.264      0.000e+00            nan
                   Spleen|         228.400      0.000e+00            nan
                   Thymus|          25.909      0.000e+00            nan
                  EyeLens|           0.456      0.000e+00            nan
                      RST|       18212.502      2.196e-02          0.387
=======================================================================

Graphics systems deleted.
Visualization Manager deleting...
G4WT0 > Destroying WorkerRunManager (0x4c3c2d920)
G4WT0 > G4 kernel has come to Quit state.
G4WT0 > UserDetectorConstruction deleted.
G4WT0 > UserPhysicsList deleted.
G4WT0 > UserActionInitialization deleted.
G4WT0 > UserWorkerInitialization deleted.
G4WT0 > UserWorkerThreadInitialization deleted.
G4WT0 > UserRunAction deleted.
G4WT0 > UserPrimaryGenerator deleted.
G4WT0 > RunManager is deleting RunManagerKernel.
G4WT0 > G4SDManager deleted.
G4WT0 > EventManager deleted.
G4WT0 > Units table cleared.
G4WT0 > TransportationManager deleted.
G4WT0 > Total navigation history collections cleaned: 5
G4WT0 > ================== Deleting memory pools ===================
G4WT0 > Pool ID '20G4NavigationLevelRep', size : 0.00481 MB
G4WT0 > Pool ID '24G4ReferenceCountedHandleIvE', size : 0.000961 MB
G4WT0 > Pool ID '7G4Event', size : 0.000961 MB
G4WT0 > Pool ID '15G4PrimaryVertex', size : 0.000961 MB
G4WT0 > Pool ID '17G4PrimaryParticle', size : 0.000961 MB
G4WT0 > Pool ID '15G4HCofThisEvent', size : 0.000961 MB
G4WT0 > Pool ID '17G4DynamicParticle', size : 0.000961 MB
G4WT0 > Pool ID '7G4Track', size : 0.000961 MB
G4WT0 > Pool ID '18G4TouchableHistory', size : 0.000961 MB
G4WT0 > Pool ID '15G4CountedObjectIvE', size : 0.000961 MB
G4WT0 > Number of memory pools allocated: 10; of which, static: 0
G4WT0 > Dynamic pools deleted: 10 / Total memory freed: 0.013 MB
G4WT0 > ============================================================
G4WT0 > G4Allocator objects are deleted.
G4WT0 > Thread-local UImanager is to be deleted.
G4WT0 > There should not be any thread-local G4cout/G4cerr hereafter.
UImanager deleted.
StateManager deleted.
RunManagerKernel is deleted. Good bye :)
RunManager is deleted.
G4 kernel has come to Quit state.
UserDetectorConstruction deleted.
UserPhysicsList deleted.
UserActionInitialization deleted.
UserWorkerInitialization deleted.
UserWorkerThreadInitialization deleted.
UserRunAction deleted.
UserPrimaryGenerator deleted.
RunManager is deleting RunManagerKernel.
G4SDManager deleted.
EventManager deleted.
Units table cleared.
TransportationManager deleted.
Total navigation history collections cleaned: 3
G4RNGHelper object is deleted.
================== Deleting memory pools ===================
Pool ID '20G4NavigationLevelRep', size : 0.003 MB
Pool ID '24G4ReferenceCountedHandleIvE', size : 0.001 MB
Pool ID '17G4DynamicParticle', size : 0.001 MB
Pool ID '16G4SmartVoxelNode', size : 1448.742 MB
Pool ID '17G4SmartVoxelProxy', size : 722.117 MB
Pool ID '15G4CountedObjectIvE', size : 0.001 MB
Number of memory pools allocated: 6; of which, static: 0
Dynamic pools deleted: 6 / Total memory freed: 2170.87 MB
============================================================
G4Allocator objects are deleted.
UImanager deleted.
StateManager deleted.
RunManagerKernel is deleted. Good bye :)
RunManager is deleted.
