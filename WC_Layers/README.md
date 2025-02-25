# README for Two Layer Spherical Shell Problem
Describes the basic parameters and order-of-operations for this problem. The boundary conditions and tallies are the same as in the single-layer model.
## 1)
- Monoenergetic 14.0MeV neutrons at the origin, incident upon nested spherical shells of W and C (by default).
- Run OpenMC and obtain neutron flux distribution over unstructured mesh (Mesh.h5).
## 2)
- Use Conversion_h5m.py and PyNE_Lib.py to prepare inputs for R2S Step 1.
- Run R2S Step 1.
## 2)
- Run ALARA using R2S Step 1 outputs.
## 3)
- Run R2S Step 2 using ALARA outputs.
## 4)
- Use Source_Mesh_Reader to extract relevant data from R2S Step 2 outputs.
## 5)
- Run OpenMC photon transport calculation (OpenMC_PhotonTransport, TwoLayers_Geometry, and TwoLayers_Materials).
## 6)
- Use Photon_TallytoVtk to convert OpenMC tally data to vtk format.
