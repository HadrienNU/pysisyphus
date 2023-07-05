import mdi
import numpy as np
from mpi4py import MPI

from pysisyphus.calculators.Calculator import Calculator


class LAMMPS(Calculator):
    def __init__(self, sigma=1.8897261251, epsilon=1, rc=None):
        super().__init__()
        # méthode TCP pour découpler les processus?
        mdiarg = "-role DRIVER -name sequence -method MPI"
        # Command pour lancer Lammps A vérifier

        # Not plugin mode
        mdi.MDI_Init(mdiarg)
        # LAMMPS engine is a stand-alone code
        # world = MPI communicator for just this driver
        # invoke perform_tasks() directly
        self.world = mdi.MDI_MPI_get_world_comm()
        self.mdicomm = mdi.MDI_Accept_Communicator()

    def get_energy(self, atoms, coords):
        energy, _ = self.calculate(coords.reshape(-1, 3))
        return {"energy": energy}

    def get_forces(self, atoms, coords):
        energy, forces = self.calculate(coords.reshape(-1, 3))
        return {
            "energy": energy,
            "forces": forces.flatten(),
        }


me = self.world.Get_rank()
nprocs = self.world.Get_size()

# allocate vectors for per-atom types, coords, vels, forces
natoms = nx * ny * nz
atypes = np.zeros(natoms, dtype=np.int)
coords = np.zeros(3 * natoms, dtype=np.float64)
vels = np.zeros(3 * natoms, dtype=np.float64)
forces = np.zeros(3 * natoms, dtype=np.float64)

atypes[:] = 1


# define simulation box

onerho = rho + (random.random() - 0.5) * rhodelta
sigma = pow(1.0 / onerho, 1.0 / 3.0)

xlo = ylo = zlo = 0.0
xhi = nx * sigma
yhi = ny * sigma
zhi = nz * sigma

# send simulation box to engine

vec = [xhi - xlo, 0.0, 0.0] + [0.0, yhi - ylo, 0.0] + [0.0, 0.0, zhi - zlo]
mdi.MDI_Send_command(">CELL", mdicomm)
mdi.MDI_Send(vec, 9, mdi.MDI_DOUBLE, mdicomm)

for m in range(3 * natoms):
    coords[m] += 2.0 * random.random() * delta - delta

# define initial velocities

for m in range(3 * natoms):
    vels[m] = 0.0

# send atoms and their properties to engine

mdi.MDI_Send_command(">NATOMS", mdicomm)
mdi.MDI_Send(natoms, 1, mdi.MDI_INT, mdicomm)
mdi.MDI_Send_command(">TYPES", mdicomm)
mdi.MDI_Send(atypes, natoms, mdi.MDI_INT, mdicomm)
mdi.MDI_Send_command(">COORDS", mdicomm)
mdi.MDI_Send(coords, 3 * natoms, mdi.MDI_DOUBLE, mdicomm)
mdi.MDI_Send_command(">VELOCITIES", mdicomm)
mdi.MDI_Send(vels, 3 * natoms, mdi.MDI_DOUBLE, mdicomm)
# request potential energy

mdi.MDI_Send_command("<PE", mdicomm)
pe = mdi.MDI_Recv(1, mdi.MDI_DOUBLE, mdicomm)
pe = world.bcast(pe, root=0)

# request forces

mdi.MDI_Send_command("<FORCES", mdicomm)
mdi.MDI_Recv(3 * natoms, mdi.MDI_DOUBLE, mdicomm, buf=forces)
world.Bcast(forces, root=0)

# final output from each calculation
# pressure = trace of virial tensor, no kinetic component

aveeng = pe / natoms

m = 0
fx = fy = fz = 0.0
for i in range(natoms):
    fx += forces[m]
    fy += forces[m + 1]
    fz += forces[m + 2]
    m += 3

fx /= natoms
fy /= natoms
fz /= natoms
