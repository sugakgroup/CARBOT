from rdkit import Chem
from rdkit.Chem.rdchem import RWMol, BondType

def Kekulize_aromatic(mol,start=0):
    rwmol = propagate_structure(mol)
    if rwmol is None:
        return None
    for i in range(start,len(rwmol.GetBonds())):
        bond = rwmol.GetBondWithIdx(i)
        if bond.GetBondType() == BondType.AROMATIC:
            for bondtype in (BondType.SINGLE, BondType.DOUBLE):
                bond.SetBondType(bondtype)
                res_mol = Kekulize_aromatic(rwmol,i)
                if not res_mol is None:
                    return res_mol
            return None
    return rwmol

def propagate_structure(mol):
    rwmol = RWMol(mol)
    atoms = set()
    for bond in rwmol.GetBonds():
        if bond.GetBondType() in {BondType.SINGLE,BondType.DOUBLE,BondType.TRIPLE}:
            atoms.add(bond.GetBeginAtom())
            atoms.add(bond.GetEndAtom())
    determined_atoms = set()
    while atoms:
        next_atoms = set()
        for atom in atoms:
            aromatic = []
            single = []
            multiple = []
            for bond in atom.GetBonds():
                if bond.GetBondType() == BondType.AROMATIC:
                    aromatic.append(bond)
                elif bond.GetBondType() == BondType.SINGLE:
                    single.append(bond)
                else:
                    multiple.append(bond)
            if len(multiple) > 1:
                return None
            elif len(multiple) == 1:
                determined_atoms.add(atom)
                for bond in aromatic:
                    bond.SetBondType(BondType.SINGLE)
                    if not bond.GetEndAtom() in determined_atoms:
                        next_atoms.add(bond.GetEndAtom())
                    if not bond.GetBeginAtom() in determined_atoms:
                        next_atoms.add(bond.GetBeginAtom())
            elif len(aromatic) == 1:
                determined_atoms.add(atom)
                bond = aromatic[0]
                bond.SetBondType(BondType.DOUBLE)
                if not bond.GetEndAtom() in determined_atoms:
                    next_atoms.add(bond.GetEndAtom())
                if not bond.GetBeginAtom() in determined_atoms:
                    next_atoms.add(bond.GetBeginAtom())
            elif len(aromatic) == 0:
                return None
        atoms = next_atoms.copy()
    return rwmol


def is_benzene_constructed(smi_parent, pos):
    mol = Chem.MolFromSmiles(smi_parent)
    begin_atom = mol.GetAtomWithIdx(pos[0])
    end_atom = mol.GetAtomWithIdx(pos[1])

    def dfs_find_len5path(now_atom, now_depth, now_visited):
        if now_depth == 5:
            if now_atom.GetIdx() == end_atom.GetIdx():
                return [now_visited]
            else:
                return None
        paths = []
        for a in now_atom.GetNeighbors():
            if a.GetIdx() in now_visited or mol.GetBondBetweenAtoms(a.GetIdx(),now_atom.GetIdx()).GetBondType() == BondType.TRIPLE:
                continue
            res = dfs_find_len5path(a, now_depth+1, now_visited+[a.GetIdx()])
            if res:
                paths += res
        return paths

    paths = dfs_find_len5path(begin_atom, 0, [begin_atom.GetIdx()])

    for p in paths:
        rwmol = Chem.RWMol(mol)
        bonds = [rwmol.GetBondBetweenAtoms(p[i], p[i + 1]) for i in range(5)]
        for i, b in enumerate(bonds):
            if i % 2 == 0:
                b.SetBondType(Chem.BondType.DOUBLE)
            else:
                b.SetBondType(Chem.BondType.SINGLE)
        kekulized_mol = Kekulize_aromatic(rwmol)
        if not kekulized_mol is None:
            return True
    return False

def num_unit_method(smi):
    mol = Chem.MolFromSmiles(smi)
    unsat = (mol.GetNumAtoms()*3 - Chem.AddHs(mol).GetNumAtoms())//2 + 1
    n_ace = 0
    for b in mol.GetBonds():
        if b.GetBondType() == Chem.BondType.TRIPLE:
            n_ace += 1
    n_eth = mol.GetNumAtoms()//2 - n_ace
    n_con = unsat - n_eth - n_ace*2
    return n_eth, n_ace, n_con

def molfromsmiles_skeleton(smi):
    mol = Chem.MolFromSmiles(smi)
    # Chem.Kekulize(mol)
    # return mol
    rwmol = RWMol(mol)
    for bond in rwmol.GetBonds():
        if bond.GetBondType() != BondType.TRIPLE:
            bond.SetBondType(BondType.SINGLE)
    for atom in rwmol.GetAtoms():
        atom.SetIsAromatic(False)
    return rwmol

def mol2skeleton(mol):
    rwmol = RWMol(mol)
    for bond in rwmol.GetBonds():
        if bond.GetBondType() != BondType.TRIPLE:
            bond.SetBondType(BondType.SINGLE)
    for atom in rwmol.GetAtoms():
        atom.SetIsAromatic(False)
    return rwmol


def fix_fixed_bonds(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = RWMol(propagate_structure(mol))
    for i in range(len(mol.GetBonds())):
        rwmol = RWMol(mol)
        bond = rwmol.GetBondWithIdx(i)
        if rwmol.GetBondWithIdx(i).GetBondType() == BondType.AROMATIC:
            bond.SetBondType(BondType.DOUBLE)
            modmol = propagate_structure(rwmol)
            if modmol is None:
                mol.GetBondWithIdx(i).SetBondType(BondType.SINGLE)
    return mol