from collections import defaultdict

from rdkit import Chem
from rdkit.Chem.rdchem import RWMol, BondType
from rdkit.Chem.rdmolops import AddHs, RemoveHs

from src.utils import *

def enumerate_EZ(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return []
    Chem.FindPotentialStereoBonds(mol)
    [bond.SetStereo(Chem.BondStereo.STEREOANY) for bond in mol.GetBonds() if (bond.GetBondType() == Chem.BondType.DOUBLE and bond.IsInRing() and bond.GetStereo() == Chem.BondStereo.STEREONONE)]
    must_setstereo_bonds = []
    for bond in mol.GetBonds():
        if bond.GetStereo() == Chem.BondStereo.STEREOANY:
            bidx = bond.GetIdx()
            atoms = (bond.GetBeginAtom(), bond.GetEndAtom())
            neighbors = []
            for a in atoms:
                for n in a.GetNeighbors():
                    if n not in atoms:
                        neighbors.append(n.GetIdx())
                        break
            must_setstereo_bonds.append(tuple([bidx])+tuple(neighbors))

    maxcount = 1 << len(must_setstereo_bonds)
    smis = set()
    for i_ez in range(maxcount):
        rwmol = Chem.RWMol(mol)
        for i, (bidx, n0, n1) in enumerate(must_setstereo_bonds):
            bond = rwmol.GetBondWithIdx(bidx)
            bond.SetStereoAtoms(n0,n1)
            if i_ez & (1 << i):
                bond.SetStereo(Chem.BondStereo.STEREOE)
            else:
                bond.SetStereo(Chem.BondStereo.STEREOZ)
        try:
            smis.add(Chem.MolToSmiles(rwmol))
        except:
            pass
    smis = list(smis)
    smis.sort()
    return smis

def evol_single(smi, operations=("connect","vinyl","ethynyl"), is_EZ=False):
    gen_smis = {k:set() for k in operations}

    original_mol = Chem.MolFromSmiles(smi)
    if original_mol is None:
        return gen_smis
    
    # handle [H][H]
    if smi == "[H][H]":
        if "vinyl" in operations:
            gen_smis["vinyl"].add("C=C")
        if "ethynyl" in operations:
            gen_smis["ethynyl"].add("C#C")
        if "phenyl" in operations:
            gen_smis["phenyl"].add("c1ccccc1")
        return {k:sorted(list(v)) for k,v in gen_smis.items()}

    # pre-calculations
    # search carbon atom with X-H bond
    atoms_XHn = []
    atoms_XHn_set = set()
    for original_atom in original_mol.GetAtoms():
        if original_atom.GetNumImplicitHs() > 0 or original_atom.GetNumExplicitHs() > 0:
            atomidx = original_atom.GetIdx()
            atoms_XHn.append(atomidx)
            atoms_XHn_set.add(atomidx)

    # find connectable pairs for annulation
    pl2_pairs = set()
    pl4_pairs = set()
    if "annulate_pl2" in operations or "annulate_pl4" in operations:
        fixed_mol = fix_fixed_bonds(original_mol)
        if "annulate_pl2" in operations:
            for bond in fixed_mol.GetBonds():
                if bond.GetBondType() in {BondType.SINGLE,BondType.AROMATIC}:
                    atom1_idx = bond.GetBeginAtomIdx()
                    atom2_idx = bond.GetEndAtomIdx()
                    atom1_neighbor_idxs = []
                    for nbr in fixed_mol.GetAtomWithIdx(atom1_idx).GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx != atom2_idx:
                            atom1_neighbor_idxs.append(nbr_idx)
                    atom2_neighbor_idxs = []
                    for nbr in fixed_mol.GetAtomWithIdx(atom2_idx).GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx != atom1_idx:
                            atom2_neighbor_idxs.append(nbr_idx)
                    for nbr1_idx in atom1_neighbor_idxs:
                        for nbr2_idx in atom2_neighbor_idxs:
                            if nbr1_idx == nbr2_idx or (not nbr1_idx in atoms_XHn_set) or (not nbr2_idx in atoms_XHn_set):
                                continue
                            check_mol = RWMol(fixed_mol)
                            check_mol.GetBondBetweenAtoms(atom1_idx, nbr1_idx).SetBondType(BondType.DOUBLE)
                            check_mol.GetBondBetweenAtoms(atom1_idx, atom2_idx).SetBondType(BondType.SINGLE)
                            check_mol.GetBondBetweenAtoms(atom2_idx, nbr2_idx).SetBondType(BondType.DOUBLE)
                            if propagate_mol(check_mol) is None:
                                continue
                            else:
                                pl2_pairs.add((nbr1_idx,atom1_idx,atom2_idx,nbr2_idx))
        if "annulate_pl4" in operations:
            for bond in fixed_mol.GetBonds():
                atom1,atom2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if atom1 in atoms_XHn_set and atom2 in atoms_XHn_set and bond.GetBondType() in {BondType.DOUBLE,BondType.AROMATIC}:
                    pl4_pairs.add((bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()))

    # operations application
    if "vinyl" in operations:
        for atom_idx in atoms_XHn:
            rwmol = RWMol(original_mol)
            new_atom_idx = [rwmol.AddAtom(Chem.Atom(6)),rwmol.AddAtom(Chem.Atom(6))]
            rwmol.AddBond(new_atom_idx[0],new_atom_idx[1],BondType.DOUBLE)
            rwmol.AddBond(atom_idx,new_atom_idx[0],BondType.SINGLE)
            rwmol.GetAtomWithIdx(atom_idx).SetNumExplicitHs(max(0,rwmol.GetAtomWithIdx(atom_idx).GetNumExplicitHs()-1))
            try:
                gen_smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(rwmol)))
            except:
                continue
            if is_EZ:
                for gen_smi_EZ in enumerate_EZ(gen_smi):                              
                    gen_smis["vinyl"].add(gen_smi_EZ)
            else:
                gen_smis["vinyl"].add(gen_smi)     

    if "ethynyl" in operations:
        for atom_idx in atoms_XHn:
            rwmol = RWMol(original_mol)
            new_atom_idx = [rwmol.AddAtom(Chem.Atom(6)),rwmol.AddAtom(Chem.Atom(6))]
            rwmol.AddBond(new_atom_idx[0],new_atom_idx[1],BondType.TRIPLE)
            rwmol.AddBond(atom_idx,new_atom_idx[0],BondType.SINGLE)
            rwmol.GetAtomWithIdx(atom_idx).SetNumExplicitHs(max(0,rwmol.GetAtomWithIdx(atom_idx).GetNumExplicitHs()-1))
            try:
                gen_smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(rwmol)))
            except:
                continue
            if is_EZ:
                for gen_smi_EZ in enumerate_EZ(gen_smi):                              
                    gen_smis["ethynyl"].add(gen_smi_EZ)
            else:
                gen_smis["ethynyl"].add(gen_smi)     

    if "connect" in operations:
        num_XH = len(atoms_XHn)
        for i in range(num_XH):
            for j in range(i+1,num_XH):
                if original_mol.GetBondBetweenAtoms(atoms_XHn[i],atoms_XHn[j]) is None:
                    rwmol = RWMol(original_mol)
                    rwmol.AddBond(atoms_XHn[i],atoms_XHn[j],BondType.SINGLE)
                    rwmol.GetAtomWithIdx(atoms_XHn[i]).SetNumExplicitHs(max(0,rwmol.GetAtomWithIdx(atoms_XHn[i]).GetNumExplicitHs()-1))
                    rwmol.GetAtomWithIdx(atoms_XHn[j]).SetNumExplicitHs(max(0,rwmol.GetAtomWithIdx(atoms_XHn[j]).GetNumExplicitHs()-1))
                    try:
                        gen_smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(rwmol)))
                    except:
                        continue
                    if is_EZ:
                        for gen_smi_EZ in enumerate_EZ(gen_smi):                              
                            gen_smis["connect"].add(gen_smi_EZ)
                    else:
                        gen_smis["connect"].add(gen_smi)

    if "annulate_pl2" in operations:
        for atom_pair in pl2_pairs:
            rwmol = RWMol(original_mol)
            new_atom_idx = [rwmol.AddAtom(Chem.Atom(6)) for _ in range(2)]+list(atom_pair)
            [rwmol.GetAtomWithIdx(idx).SetIsAromatic(True) for idx in new_atom_idx]
            [rwmol.AddBond(new_atom_idx[i-1],new_atom_idx[i],BondType.AROMATIC) for i in range(3)]
            [rwmol.GetBondBetweenAtoms(new_atom_idx[i-1],new_atom_idx[i]).SetBondType(BondType.AROMATIC) for i in range(3,6)]
            [rwmol.GetBondBetweenAtoms(new_atom_idx[i-1],new_atom_idx[i]).SetIsAromatic(True) for i in range(6)]
            [rwmol.GetAtomWithIdx(atom_idx).SetNumExplicitHs(max(0,rwmol.GetAtomWithIdx(atom_idx).GetNumExplicitHs()-1)) for atom_idx in {atom_pair[0],atom_pair[-1]}]
            try:
                gen_smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(rwmol)))
            except:
                continue
            if is_EZ:
                for gen_smi_EZ in enumerate_EZ(gen_smi):                              
                    gen_smis["annulate_pl2"].add(gen_smi_EZ)
            else:
                gen_smis["annulate_pl2"].add(gen_smi)

    if "annulate_pl4" in operations:
        for atom_pair in pl4_pairs:
            rwmol = RWMol(original_mol)
            new_atom_idx = [rwmol.AddAtom(Chem.Atom(6)) for _ in range(4)]+list(atom_pair)
            [rwmol.GetAtomWithIdx(idx).SetIsAromatic(True) for idx in new_atom_idx]
            [rwmol.AddBond(new_atom_idx[i-1],new_atom_idx[i],BondType.AROMATIC) for i in range(5)]
            [rwmol.GetBondBetweenAtoms(new_atom_idx[i-1],new_atom_idx[i]).SetBondType(BondType.AROMATIC) for i in range(5,6)]
            [rwmol.GetBondBetweenAtoms(new_atom_idx[i-1],new_atom_idx[i]).SetIsAromatic(True) for i in range(6)]
            [rwmol.GetAtomWithIdx(atom_idx).SetNumExplicitHs(max(0,rwmol.GetAtomWithIdx(atom_idx).GetNumExplicitHs()-1)) for atom_idx in {atom_pair[0],atom_pair[-1]}]
            try:
                gen_smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(rwmol)))
            except:
                continue
            if is_EZ:
                for gen_smi_EZ in enumerate_EZ(gen_smi):                              
                    gen_smis["annulate_pl4"].add(gen_smi_EZ)
            else:
                gen_smis["annulate_pl4"].add(gen_smi)

    if "phenyl" in operations:
        for atom_idx in atoms_XHn:
            rwmol = RWMol(original_mol)
            new_atom_idx = [rwmol.AddAtom(Chem.Atom(6)) for _ in range(6)]
            [rwmol.GetAtomWithIdx(idx).SetIsAromatic(True) for idx in new_atom_idx]
            rwmol.AddBond(atom_idx,new_atom_idx[0],BondType.SINGLE)
            [rwmol.AddBond(new_atom_idx[i-1],new_atom_idx[i],BondType.AROMATIC) for i in range(6)]
            [rwmol.GetBondBetweenAtoms(new_atom_idx[i-1],new_atom_idx[i]).SetIsAromatic(True) for i in range(6)]
            rwmol.GetAtomWithIdx(atom_idx).SetNumExplicitHs(max(0,rwmol.GetAtomWithIdx(atom_idx).GetNumExplicitHs()-1))
            try:
                gen_smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(rwmol)))
            except:
                continue
            if is_EZ:
                for gen_smi_EZ in enumerate_EZ(gen_smi):                              
                    gen_smis["phenyl"].add(gen_smi_EZ)
            else:
                gen_smis["phenyl"].add(gen_smi)
    
    return {k:sorted(list(v)) for k,v in gen_smis.items()}

if __name__ == "__main__":
    print(evol_single("C=C1C=CC=CC=CC=CC=CC=Cc2cc1ccc2",operations=("connect",),is_EZ=False))



    






    
    





