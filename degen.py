from collections import defaultdict
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from utils import *

rxn_smas = {
    "vinyl": "[C,c:1]-[C;h1]=[C;h2]>>[C,c;h:1]",
    "ethynyl": "[C,c:1]-C#[C;h1]>>[C,c;h:1]",
    "phenyl": "[C,c:1]-c1[c;h1][c;h1][c;h1][c;h1][c;h1]1>>[C,c;h:1]",
    "annulate_pl2": "[c:1]1[c:2][c:3][c:4][c;h1][c;h1]1>>[C:1]=[C:2][C:3]=[C:4]",
    "annulate_pl4": "[c:1]1[c:2][c;h1][c;h1][c;h1][c;h1]1>>[C:1]=[C:2]",
}
rxns = {m: AllChem.ReactionFromSmarts(rxn_smas[m]) for m in rxn_smas}

def degenerate_rxn(mol, method_name):
    rs = rxns[method_name]
    products = set()
    for m in rs.RunReactants([mol]):
        rxmol = Kekulize_aromatic(m[0])
        if not rxmol is None:
            products.add(Chem.MolToSmiles(rxmol))
    return sorted(list(products))


def degen_single(smi, evolmethods=("connect","vinyl","ethynyl")):
    parents = []
    if smi == "[H][H]":
        return []
    if smi == "C=C" and "vinyl" in evolmethods:
        return [("[H][H]",("C=C","vinyl"))]
    if smi == "C#C" and "ethynyl" in evolmethods:
        return [("[H][H]",("C#C","ethynyl"))]
    if smi == "c1ccccc1" and "phenyl" in evolmethods:
        parents.append(("[H][H]",("c1ccccc1","phenyl")))

    if "connect" in evolmethods:
        parents_set = set()
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return parents
        
        mol = RWMol(mol)
        for atom in mol.GetAtoms():
            atom.SetIsAromatic(False)

        # if a bond is located in a ring and is single/aromatic, it is cuttable
        for bond in mol.GetBonds():
            if (not bond.IsInRing()) or (not bond.GetBondType() in {Chem.BondType.SINGLE,Chem.BondType.AROMATIC}):
                continue
            aidx1 = bond.GetBeginAtomIdx()
            aidx2 = bond.GetEndAtomIdx()
            rwmol = Chem.RWMol(mol)
            if bond.GetBondType() == Chem.BondType.AROMATIC:
                # Aromatic -> Single
                rwmol.GetBondBetweenAtoms(aidx1, aidx2).SetBondType(Chem.BondType.SINGLE)
            rwmol = Kekulize_aromatic(rwmol)
            if rwmol is None:
                continue
            rwmol.RemoveBond(aidx1, aidx2)
            try:
                Chem.SanitizeMol(rwmol)
                parents_set.add(Chem.MolToSmiles(rwmol))
            except:
                pass

        for parent_smi in sorted(list(parents_set)):
            parents.append((parent_smi,(smi,"connect")))

    mol = Chem.MolFromSmiles(smi)

    if "vinyl" in evolmethods:
        for parent_smi in degenerate_rxn(mol, "vinyl"):
            parents.append((parent_smi,(smi,"vinyl")))
    if "ethynyl" in evolmethods:
        for parent_smi in degenerate_rxn(mol, "ethynyl"):
            parents.append((parent_smi,(smi,"ethynyl")))
    if "annulate_pl2" in evolmethods:
        for parent_smi in degenerate_rxn(mol, "annulate_pl2"):
            parents.append((parent_smi,(smi,"annulate_pl2")))
    if "annulate_pl4" in evolmethods:
        for parent_smi in degenerate_rxn(mol, "annulate_pl4"):
            parents.append((parent_smi,(smi,"annulate_pl4")))
    if "phenyl" in evolmethods:
        for parent_smi in degenerate_rxn(mol, "phenyl"):
            parents.append((parent_smi,(smi,"phenyl")))

    return parents

if __name__ == "__main__":
    print(degen_single("C=CC1=CC(C#C)=CC2=CC#CC=C12",evolmethods=("connect","vinyl","ethynyl","annulate_pl2","annulate_pl4","phenyl")))






    









    

