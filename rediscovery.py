from rdkit import Chem
from rdkit.Chem.rdchem import BondStereo
from utils import *
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

class Rediscovery:
    def __init__(self,target_smi,weak=False,is_EZ=False):
        self.weak = weak
        self.is_EZ = is_EZ
        self.is_bond_fixed = False
        try:
            self.target_smi = Chem.MolToSmiles(Chem.MolFromSmiles(target_smi))
            self.target_mol = fix_fixed_bonds(target_smi)
            for bond in Chem.MolFromSmiles(target_smi).GetBonds():
                if bond.IsInRing() and bond.GetBondType() == BondType.DOUBLE:
                    self.is_bond_fixed = True
                    break
            self.target_mol_skeleton = mol2skeleton(self.target_mol)
            self.target_generation = sum(num_unit_method(target_smi))
        except:
            print("Error! Please change the target molecule.")
            exit()
    
    def score(self,smi):
        if smi == "[H][H]":
            return (sum(num_unit_method(smi))+1)/(self.target_generation+1)
        try:
            if self.is_bond_fixed:
                mol = fix_fixed_bonds(smi)
            else:
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
                            if not rwmol.GetBondBetweenAtoms(i,j) is None:
                                bond_list.append((i,j,atoms[i],atoms[j]))
                    # stereo
                    if self.is_EZ:
                        is_stereo_correct = True
                        for i,j,i_target,j_target in bond_list:
                            target_stereo = self.target_mol.GetBondBetweenAtoms(i_target,j_target).GetStereo()
                            parent_stereo = rwmol.GetBondBetweenAtoms(i,j).GetStereo()
                            if all([bond in (BondStereo.STEREOE,BondStereo.STEREOZ) for bond in (target_stereo, parent_stereo)]):
                                k,l = list(rwmol.GetBondBetweenAtoms(i,j).GetStereoAtoms())
                                if 1^(i == k)^(j == k)^(i == l)^(j == l)^(target_stereo == parent_stereo):
                                    is_stereo_correct = False
                                    break
                        if not is_stereo_correct:
                            continue
                    
                    for i,j,i_target,j_target in bond_list:
                        bond = rwmol.GetBondBetweenAtoms(i,j)
                        target_bond = self.target_mol.GetBondBetweenAtoms(i_target,j_target)
                        if bond.GetBondType() != BondType.AROMATIC and target_bond.GetBondType() != BondType.AROMATIC and bond.GetBondType() != target_bond.GetBondType():
                            rwmol = None
                            break
                        bond.SetBondType(target_bond.GetBondType())
                    rwmol = Kekulize_aromatic(rwmol)
                    if rwmol is None:
                        continue

                    target_rwmol = RWMol(self.target_mol)
                    for i,j,i_target,j_target in bond_list:
                        target_rwmol.GetBondBetweenAtoms(i_target,j_target).SetBondType(rwmol.GetBondBetweenAtoms(i,j).GetBondType())
                    if not Kekulize_aromatic(target_rwmol) is None:
                        return (sum(num_unit_method(smi))+1)/(self.target_generation+1)
                return 0
        except:
            print("structure construction failed!")
            return 0
        

class RediscoveryTanimoto:
    def __init__(self,target_smi):
        try:
            self.target_smi = target_smi
            self.target_fp = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(self.target_smi), 2)
            if self.target_fp is None:
                print("Error! Please change the target molecule.")
                exit()
        except:
            print("Error! Please change the target molecule.")
            exit()
    
    def score(self,smi):
        try:
            fp = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi), 2)
            return TanimotoSimilarity(self.target_fp,fp)
        except:
            print("structure construction failed!")
            return 0

if __name__ == "__main__":
    rds = Rediscovery("C1(C=C2C=CC(C=C2)=CC3=CC=C4C=C3)=CC=C(C=C5C=CC(C=C5)=CC6=CC=C(C=C6)C=C(C=C7)C=CC7=C4)C=C1",is_EZ=True) #AnFLAP_EZ
    print(rds.score("C=C/C=C/C=C\\C=C\\C=C/C=C\\C=C/C=C/C=C/C=C\\C=C\\C=C"))
