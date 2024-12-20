from pyne.material import Material
from pyne.material_library import MaterialLibrary

W=Material({"W":1})
W.density=19.3
Graphite=Material({"C":1})
Graphite.density=2.62

MatLib = MaterialLibrary()
MatLib["W"]=W
MatLib["Graphite"]=Graphite
MatLib.write_hdf5("PyNE_Lib_WC.h5", h5_overwrite=True)
