from collections import defaultdict
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from src.utils import *

DEGEN_QUERIES = {
    "vinyl": {
        "query": Chem.MolFromSmarts("[!#1]-[#6X3H1!a]=[#6X3H2!a]"),
        "roots": (0,),
        "remove_bonds": ((0,1),(1,2)),
        "remove_atoms": (1,2),
    },
    "ethynyl": {
        "query": Chem.MolFromSmarts("[!#1]-[#6X2H0!a]#[#6X2H1!a]"),
        "roots": (0,),
        "remove_bonds": ((0,1),(1,2)),
        "remove_atoms": (1,2),
    },
    "connect": {
        "query": Chem.MolFromSmarts("[!#1]-@,:@[!#1]"),
        "roots": (0,1),
        "remove_bonds": ((0,1),),
        "remove_atoms": (),
    },
    "annulate_pl2": {
        "query": Chem.MolFromSmarts("[#6X3a]1:[#6X3H1a]:[#6X3H1a]:[#6X3a]:[#6X3a]:[#6X3a]:1"),
        "roots": (0,3),
        "remove_bonds": ((0,1),(1,2),(2,3)),
        "remove_atoms": (1,2),
    },
    "annulate_pl4": {
        "query": Chem.MolFromSmarts("[#6X3a]1:[#6X3H1a]:[#6X3H1a]:[#6X3H1a]:[#6X3H1a]:[#6X3a]:1"),
        "roots": (0,5),
        "remove_bonds": ((0,1),(1,2),(2,3),(3,4),(4,5)),
        "remove_atoms": (1,2,3,4)
    },
    "phenyl": {
        "query": Chem.MolFromSmarts("[!#1]-[#6X3a]1:[#6X3H1a]:[#6X3H1a]:[#6X3H1a]:[#6X3H1a]:[#6X3H1a]:1"),
        "roots": (0,),
        "remove_bonds": ((0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,1)),
        "remove_atoms": (1,2,3,4,5,6)
    }
}

def degenerate_rxn(mol, method_name):
    degen_flow = DEGEN_QUERIES[method_name]
    products = set()
    matches = mol.GetSubstructMatches(degen_flow["query"], uniquify=False)
    for match in matches:
        rwmol = RWMol(mol)
        # remove bonds
        bonds_to_remove = [(match[bond[0]], match[bond[1]]) for bond in degen_flow["remove_bonds"]]
        for atom1, atom2 in bonds_to_remove:
            rwmol.RemoveBond(atom1, atom2)
        # increase hydrogens to avoid valence issues
        for idx in degen_flow["roots"]:
            atom = rwmol.GetAtomWithIdx(match[idx])
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
        # remove unnecessary atoms
        atoms_to_remove = sorted([match[idx] for idx in degen_flow["remove_atoms"]], reverse=True)
        for idx in atoms_to_remove:
            rwmol.RemoveAtom(idx)
        rwmol = propagate_mol(rwmol)
        if not rwmol is None:
            for atom in rwmol.GetAtoms():
                if all([bond.GetBondType() != BondType.AROMATIC for bond in atom.GetBonds()]):
                    atom.SetIsAromatic(False)
            for bond in rwmol.GetBonds():
                if bond.GetBondType() != BondType.AROMATIC:
                    bond.SetIsAromatic(False)
            try:
                products.add(Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(rwmol))))
            except:
                pass
    return sorted(list(products))

def degen_single(smi, operations=("connect","vinyl","ethynyl")):
    original_mol = Chem.MolFromSmiles(smi)
    if any([operation in {"connect", "annulate_pl2", "annulate_pl4"} for operation in operations]):
        original_mol = fix_fixed_bonds(original_mol)

    parents = {k:set() for k in operations}
    if smi == "[H][H]":
        return parents
    if smi == "C=C" and "vinyl" in operations:
        parents["vinyl"].add("[H][H]")
    if smi == "C#C" and "ethynyl" in operations:
        parents["ethynyl"].add("[H][H]")
    if smi == "c1ccccc1" and "phenyl" in operations:
        parents["phenyl"].add("[H][H]")

    for operation in operations:
        for smi in degenerate_rxn(original_mol, operation):
            parents[operation].add(smi)

    return {k:sorted(list(v)) for k,v in parents.items()}

if __name__ == "__main__":
    print(degen_single("C=Cc1c(C#C)cccc1",operations=("connect","vinyl","ethynyl","annulate_pl2","annulate_pl4","phenyl")))








    









    

