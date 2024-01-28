from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV, EVANG2AUBOHR
from pysisyphus.Geometry import Geometry


class ASECalc(Calculator):
    def __init__(self, ase_calc):
        super().__init__()
        self.calc = ase_calc
        # Dummies
        self.mult = 1
        self.charge = 0

    def get_energy(self, atoms, coords):
        try:
            import ase
        except ImportError:
            print("Please install the 'ase' package!")
            return None
        ase_atoms = ase.Atoms(symbols=atoms, positions=coords.reshape(-1, 3) * BOHR2ANG)
        energy = self.calc.get_potential_energy(ase_atoms) / AU2EV
        return {
            "energy": energy,
        }

    def get_forces(self, atoms, coords):
        try:
            import ase
        except ImportError:
            print("Please install the 'ase' package!")
            return None
        ase_atoms = ase.Atoms(symbols=atoms, positions=coords.reshape(-1, 3) * BOHR2ANG)
        forces = self.calc.get_forces(ase_atoms) * EVANG2AUBOHR  # Do Units conversion
        results = {
            "forces": forces,
        }
        results.update({"energy": self.calc.get_potential_energy(atoms) / AU2EV})
        return results

    @classmethod
    def get_geom_from_ase(cls, ase_atoms, calc=None, geom_kwargs=None):
        if calc is None:
            calc = ase_atoms.calc
        if geom_kwargs is None:
            geom_kwargs = dict()

        geom = Geometry(ase_atoms.species, ase_atoms.positions * ANG2BOHR, **geom_kwargs)
        geom.set_calculator(cls(calc))
        return geom
