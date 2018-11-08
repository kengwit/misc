from mechanics.mechanics import Mechanics
from mechanics.mesh_reader import MeshReadYAML

infile = "./examples/tet2.yaml"

mec = Mechanics()
MeshReadYAML(mec,infile)


mec.mech_implicit_solve()


