from rdkit import Chem
from rdkit.Chem.rdchem import RWMol, BondType

def fix_fixed_bonds(mol):
    modified_mol = RWMol(mol)
    decided_bond_order = []
    for bond in modified_mol.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type == BondType.AROMATIC:
            decided_bond_order.append(0)
        else:
            decided_bond_order.append(3)
    for bondidx, bond_order in enumerate(decided_bond_order):
        if bond_order&1 == 0:
            checkmol = RWMol(modified_mol)
            checkmol.GetBondWithIdx(bondidx).SetBondType(BondType.SINGLE)
            propagated_mol = propagate_mol(checkmol)
            if not propagated_mol is None:
                for i, propagated_bond in enumerate(propagated_mol.GetBonds()):
                    if propagated_bond.GetBondType() == BondType.SINGLE:
                        decided_bond_order[i] |= 1
                    elif propagated_bond.GetBondType() == BondType.DOUBLE:
                        decided_bond_order[i] |= 2
        if bond_order&2 == 0:
            checkmol = RWMol(modified_mol)
            checkmol.GetBondWithIdx(bondidx).SetBondType(BondType.DOUBLE)
            propagated_mol = propagate_mol(checkmol)
            if not propagated_mol is None:
                for i, propagated_bond in enumerate(propagated_mol.GetBonds()):
                    if propagated_bond.GetBondType() == BondType.SINGLE:
                        decided_bond_order[i] |= 1
                    elif propagated_bond.GetBondType() == BondType.DOUBLE:
                        decided_bond_order[i] |= 2
        if decided_bond_order[bondidx] == 1:
            modified_mol.GetBondWithIdx(bondidx).SetBondType(BondType.SINGLE)
        elif decided_bond_order[bondidx] == 2:
            modified_mol.GetBondWithIdx(bondidx).SetBondType(BondType.DOUBLE)
    return modified_mol

def propagate_mol(mol):
    modified_mol = RWMol(mol)
    atom_idxs_tobe_checked = set([i for i in range(len(modified_mol.GetAtoms()))])
    determined_atoms = set()
    while atom_idxs_tobe_checked:
        next_atom_idxs = set()
        for atom_idx in atom_idxs_tobe_checked:
            atom = modified_mol.GetAtomWithIdx(atom_idx)
            total_valence = atom.GetNumImplicitHs()+atom.GetNumExplicitHs()
            aromatic_neighbors = []
            for neighbor_atom in atom.GetNeighbors():
                neighbor_atom_idx = neighbor_atom.GetIdx()
                bondType = modified_mol.GetBondBetweenAtoms(atom_idx, neighbor_atom_idx).GetBondType()
                if bondType == BondType.AROMATIC:
                    aromatic_neighbors.append(neighbor_atom_idx)
                elif bondType == BondType.SINGLE:
                    total_valence += 1
                elif bondType == BondType.DOUBLE:
                    total_valence += 2
                elif bondType == BondType.TRIPLE:
                    total_valence += 3
            if total_valence+len(aromatic_neighbors) > 4: # exceeded valence (NG)
                return None
            elif total_valence+len(aromatic_neighbors) == 4: # determine all aromatic bonds to single bonds (OK)
                for neighbor_atom_idx in aromatic_neighbors:
                    modified_mol.GetBondBetweenAtoms(atom_idx, neighbor_atom_idx).SetIsAromatic(False)
                    modified_mol.GetBondBetweenAtoms(atom_idx, neighbor_atom_idx).SetBondType(BondType.SINGLE)
                    if not neighbor_atom_idx in determined_atoms:
                        next_atom_idxs.add(neighbor_atom_idx)
            elif len(aromatic_neighbors) == 0 and total_valence != 4: # valence not satisfied (NG)
                return None
            elif len(aromatic_neighbors) == 1: 
                if total_valence == 2: # determine the aromatic bond to double bond (OK)
                    modified_mol.GetBondBetweenAtoms(atom_idx, aromatic_neighbors[0]).SetIsAromatic(False)
                    modified_mol.GetBondBetweenAtoms(atom_idx, aromatic_neighbors[0]).SetBondType(BondType.DOUBLE)
                    if not aromatic_neighbors[0] in determined_atoms:
                        next_atom_idxs.add(aromatic_neighbors[0])
                else: # cannot determine the aromatic bond (NG)
                    return None
            elif len(aromatic_neighbors) > 1: # cannot determine the aromatic bonds (pass)
                continue
            determined_atoms.add(atom_idx)
        atom_idxs_tobe_checked = next_atom_idxs
    return modified_mol

def num_unit_method(smi):
    if smi == "[H][H]":
        return 0,0,0
    mol = Chem.MolFromSmiles(smi)
    unsat = (mol.GetNumAtoms()*3 - Chem.AddHs(mol).GetNumAtoms())//2 + 1
    n_eth = len([0 for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.TRIPLE])
    n_vin = mol.GetNumAtoms()//2 - n_eth
    n_con = unsat - n_vin - n_eth*2
    return n_con, n_vin, n_eth

def num_unit_method_ensemble(smi_ensemble):
    n_con_total = 0
    n_vin_total = 0
    n_eth_total = 0
    for smi in smi_ensemble.split("."):
        n_con, n_vin, n_eth = num_unit_method(smi)
        n_con_total += n_con
        n_vin_total += n_vin
        n_eth_total += n_eth
    return n_con_total, n_vin_total, n_eth_total, -len(smi_ensemble.split("."))

def mol2skeleton(mol):
    rwmol = RWMol(mol)
    for bond in rwmol.GetBonds():
        if bond.GetBondType() != BondType.TRIPLE:
            bond.SetBondType(BondType.SINGLE)
    for atom in rwmol.GetAtoms():
        atom.SetIsAromatic(False)
    return rwmol

if __name__ == "__main__":
    mol = Chem.MolFromSmiles("c1cc2ccc3c4cc3c2cc14")
    mol.Debug()
    fixed_mol = fix_fixed_bonds(mol)
    fixed_mol.Debug()

