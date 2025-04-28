from rdkit import Chem
from utils import *

class Rediscovery:
    def __init__(self,target_smi,weak=False):
        self.weak = weak
        try:
            self.target_smi = Chem.MolToSmiles(Chem.MolFromSmiles(target_smi))
            self.target_mol = fix_fixed_bonds(target_smi)
            self.target_mol_skeleton = mol2skeleton(self.target_mol)
            self.target_generation = sum(num_unit_method(target_smi))
        except:
            print("Error! Please change the target molecule.")
            exit()
    
    def score(self,smi):
        try:
            mol = Chem.MolFromSmiles(smi)
            mol_skeleton = mol2skeleton(mol)
            if self.weak:
                atoms = self.target_mol_skeleton.GetSubstructMatch(mol_skeleton)
                if len(atoms) == 0:
                    return 0
                rwmol = RWMol(mol)
                bond_list = []
                for i in range(len(atoms)-1):
                    for j in range(i+1,len(atoms)):
                        bond = rwmol.GetBondBetweenAtoms(i,j)
                        if not bond is None:
                            bond_list.append((i,j,atoms[i],atoms[j]))
                for i,j,i_target,j_target in bond_list:
                    rwmol.GetBondBetweenAtoms(i,j).SetBondType(self.target_mol.GetBondBetweenAtoms(i_target,j_target).GetBondType())
                rwmol = Kekulize_aromatic(rwmol)
                if rwmol is None:
                    return 0
                target_rwmol = RWMol(self.target_mol)
                for i,j,i_target,j_target in bond_list:
                    target_rwmol.GetBondBetweenAtoms(i_target,j_target).SetBondType(rwmol.GetBondBetweenAtoms(i,j).GetBondType())
                target_rwmol = Kekulize_aromatic(target_rwmol)
                if not target_rwmol is None:
                    return (sum(num_unit_method(smi))+1)/(self.target_generation+1)
                return 0
            else:
                for atoms in self.target_mol_skeleton.GetSubstructMatches(mol_skeleton):
                    rwmol = RWMol(mol)
                    bond_list = []
                    for i in range(len(atoms)-1):
                        for j in range(i+1,len(atoms)):
                            bond = rwmol.GetBondBetweenAtoms(i,j)
                            if not bond is None:
                                bond_list.append((i,j,atoms[i],atoms[j]))
                    for i,j,i_target,j_target in bond_list:
                        rwmol.GetBondBetweenAtoms(i,j).SetBondType(self.target_mol.GetBondBetweenAtoms(i_target,j_target).GetBondType())
                    rwmol = Kekulize_aromatic(rwmol)
                    if rwmol is None:
                        continue
                    target_rwmol = RWMol(self.target_mol)
                    for i,j,i_target,j_target in bond_list:
                        target_rwmol.GetBondBetweenAtoms(i_target,j_target).SetBondType(rwmol.GetBondBetweenAtoms(i,j).GetBondType())
                    target_rwmol = Kekulize_aromatic(target_rwmol)
                    if not target_rwmol is None:
                        return (sum(num_unit_method(smi))+1)/(self.target_generation+1)
                return 0
        except:
            print("structure construction failed!")
            return 0

if __name__ == "__main__":
    rds = Rediscovery("c1c2/C=C\c3cc4cc5ccccc5cc4cc3/C=C\c2cc2cc3ccccc3cc21") #AnFLAP_EZ
    print(rds.score("C=C/C=C1/C=C/c2cc(=C)c(/C=C/C=C)c/c2=C/C=C/1C=C"))