from rdkit import Chem
from rdkit.Chem.rdchem import BondType, RWMol
from src.utils import *
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import rdFingerprintGenerator

class RediscoverySubstructure:
    def __init__(self,target_smi,highsymm=False):
        self.highsymm = highsymm
        self.is_bond_fixed = False
        try:
            self.target_smi = Chem.MolToSmiles(Chem.MolFromSmiles(target_smi))
            self.target_mol = fix_fixed_bonds(Chem.MolFromSmiles(target_smi))
            self.target_mol_skeleton = mol2skeleton(self.target_mol)
            self.target_generation = sum(num_unit_method(target_smi))
        except:
            print("Error! Please change the target molecule.")
            exit()
    
    def score(self,smi):
        if smi == "[H][H]":
            return 0.0
        try:
            mol = Chem.MolFromSmiles(smi)
            mol_skeleton = mol2skeleton(mol)
        except:
            return -1.0
        if self.highsymm:
            target_atom_idxs = self.target_mol_skeleton.GetSubstructMatch(mol_skeleton)
            if len(target_atom_idxs) == 0:
                return -1.0
            if substructure_search_single(mol,self.target_mol,target_atom_idxs):
                return sum(num_unit_method(smi))/self.target_generation
            else:
                return -1.0
        else:
            for mapping in self.target_mol_skeleton.GetSubstructMatches(mol_skeleton,uniquify=False):
                if substructure_search_single(mol,self.target_mol,mapping):
                    return sum(num_unit_method(smi))/self.target_generation
            return -1.0

def substructure_search_single(mol,target_mol,mapping):
    # substructure search start
    mapping_set = set(mapping)
    target_rwmol = RWMol(target_mol)
    for atom_idx, atom in enumerate(mol.GetAtoms()):
        target_atom_idx = mapping[atom_idx]
        for neighbor_atom in atom.GetNeighbors():
            neighbor_atom_idx = neighbor_atom.GetIdx()
            target_rwmol.RemoveBond(mapping[atom_idx], mapping[neighbor_atom_idx])
        for target_neighbor_atom in target_rwmol.GetAtomWithIdx(target_atom_idx).GetNeighbors():
            target_neighbor_atom_idx = target_neighbor_atom.GetIdx()
            if target_rwmol.GetBondBetweenAtoms(target_atom_idx, target_neighbor_atom_idx).GetBondType() in {BondType.DOUBLE,BondType.TRIPLE}:
                return False
            target_rwmol.RemoveBond(target_atom_idx, target_neighbor_atom_idx)
            if target_neighbor_atom_idx not in mapping_set:
                target_neighbor_atom.SetNumExplicitHs(target_neighbor_atom.GetNumExplicitHs() + 1)
    target_atom_list = sorted(list(mapping_set), reverse=True)
    for target_atom_idx in target_atom_list:
        target_rwmol.RemoveAtom(target_atom_idx)
    if not propagate_mol(target_rwmol) is None:
        return True
    else:
        return False

class RediscoveryTanimoto:
    def __init__(self,target_smi):
        try:
            self.target_smi = target_smi
            self.fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
            self.target_fp = self.fpgen.GetSparseCountFingerprint(Chem.MolFromSmiles(self.target_smi))
            if self.target_fp is None:
                print("Error! Please change the target molecule.")
                exit()
        except:
            print("Error! Please change the target molecule.")
            exit()
    
    def score(self,smi):
        try:
            fp = self.fpgen.GetSparseCountFingerprint(Chem.MolFromSmiles(smi))
            return TanimotoSimilarity(self.target_fp,fp)
        except:
            return -1.0
        
class RediscoverySubstructureEnsemble:
    def __init__(self,target_smi,highsymm=False):
        self.highsymm = highsymm
        self.is_bond_fixed = False
        try:
            self.target_smi = Chem.MolToSmiles(Chem.MolFromSmiles(target_smi))
            self.target_mol = fix_fixed_bonds(Chem.MolFromSmiles(target_smi))
            self.target_mol_skeleton = mol2skeleton(self.target_mol)
            self.target_generation = sum(num_unit_method(target_smi))
        except:
            print("Error! Please change the target molecule.")
            exit()
    
    def score(self,smi):
        if smi == "[H][H]":
            return 0.0
        try:
            smi = ".".join([part for part in smi.split(".") if not part in {"C=C","C#C"}])
            mol = Chem.MolFromSmiles(smi)
            mol_skeleton = mol2skeleton(mol)
            if self.highsymm:
                target_atom_idxs = self.target_mol_skeleton.GetSubstructMatch(mol_skeleton)
                if len(target_atom_idxs) == 0:
                    return -1.0
                if substructure_search_single(mol,self.target_mol,target_atom_idxs):
                    return sum(num_unit_method(smi))/self.target_generation
                else:
                    return -1.0
            else:
                for mapping in self.target_mol_skeleton.GetSubstructMatches(mol_skeleton):
                    if substructure_search_single(mol,self.target_mol,mapping):
                        return sum(num_unit_method(smi))/self.target_generation
                return -1.0
        except:
            return -1.0

if __name__ == "__main__":
    rds = RediscoverySubstructure("C1=C2C=c3ccccc3=C2c2ccccc21")
    print(rds.score("C=CC=C1C(=C)C=C2C=c3ccccc3=C21"))

