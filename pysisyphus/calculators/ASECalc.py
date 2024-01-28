from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG


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
        energy = self.calc.get_potential_energy(ase_atoms)
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
        forces = self.calc.get_forces(ase_atoms)
        results = {
            "forces": forces,
        }
        results.update(self.calc.get_potential_energy(atoms))
        return results
