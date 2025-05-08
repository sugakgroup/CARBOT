from rdkit import Chem
from rdkit.Chem.rdchem import RWMol, BondType
from rdkit.Chem.rdmolops import AddHs, RemoveHs

from utils import *

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
            [a.SetAtomMapNum(0) for a in rwmol.GetAtoms()]
            smis.add(Chem.MolToSmiles(rwmol))
        except:
            pass

    smis = list(smis)
    smis.sort()
    return smis


def evol_single(smi, evolmethods=("connect","vinyl","ethynyl"), is_EZ=True):
    try:
        mother_mol = AddHs(Chem.MolFromSmiles(smi))
    except:
        return []
    
    generation_information = []

    # mapping of atom indeces
    for i, atom in enumerate(mother_mol.GetAtoms(),start=1):
        atom.SetAtomMapNum(i)

    # search carbon atom with C-H bond
    atoms_HH = []
    atoms_CH = []
    atoms_CH2 = []
    for mother_atom in mother_mol.GetAtoms():
        if mother_atom.GetAtomicNum() == 6:
            mother_atom_index = mother_atom.GetIdx()
            mother_atom_H_index = [-1,-1]
            for bond in mother_atom.GetBonds():
                if bond.GetEndAtom().GetAtomicNum() == 1:
                    if mother_atom_H_index[0] == -1:
                        mother_atom_H_index[0] = bond.GetEndAtomIdx()
                    else:
                        mother_atom_H_index[1] = bond.GetEndAtomIdx()
            if mother_atom_H_index[0] != -1:
                atoms_CH.append((mother_atom_index,mother_atom_H_index[0]))
            if mother_atom_H_index[1] != -1:
                atoms_CH2.append((mother_atom_index,mother_atom_H_index[0],mother_atom_H_index[1]))

    atoms_HH = [(0,1)] if smi == "[H][H]" else []

    # connect two C-H carbons
    if "connect" in evolmethods or "annulate_pl2" in evolmethods or "annulate_pl4" in evolmethods:
        # search C-C bonds
        ccbonds = set()
        for bond in mother_mol.GetBonds():
            atom_a = bond.GetBeginAtom()
            atom_b = bond.GetEndAtom()
            if atom_a.GetAtomicNum() != 6 or atom_b.GetAtomicNum() != 6:
                continue
            bond_atoms = [atom_a.GetIdx(),atom_b.GetIdx()]
            bond_atoms.sort()
            ccbonds.add(tuple(bond_atoms))
        # connect C-C w/o bond
        if "connect" in evolmethods:
            num_CH = len(atoms_CH)
            for i in range(num_CH):
                for j in range(i+1,num_CH):
                    begin_atom = atoms_CH[i]
                    end_atom = atoms_CH[j]
                    if (begin_atom[0], end_atom[0]) in ccbonds:
                        continue
                    rw_mol = RWMol(mother_mol)
                    rw_mol.AddBond(begin_atom[0],end_atom[0],BondType.SINGLE)
                    rw_mol.RemoveBond(begin_atom[0],begin_atom[1])
                    rw_mol.RemoveBond(end_atom[0],end_atom[1])
                    if end_atom[1] > begin_atom[1]:
                        rw_mol.RemoveAtom(end_atom[1])
                        rw_mol.RemoveAtom(begin_atom[1])
                    else:
                        rw_mol.RemoveAtom(begin_atom[1])
                        rw_mol.RemoveAtom(end_atom[1])
                    [a.SetAtomMapNum(0) for a in rw_mol.GetAtoms()]
                    try:
                        smi_before = Chem.MolToSmiles(RemoveHs(rw_mol))
                    except Exception as e:
                        continue
                    if is_EZ:
                        for gen_smi in enumerate_EZ(smi_before):                                   
                            generation_information.append((gen_smi,(smi,"connect",begin_atom[0],end_atom[0])))
                    else:
                        generation_information.append((smi_before,(smi,"connect",begin_atom[0],end_atom[0])))

    
    # elongate conjugation from C-H carbons
    if "vinyl" in evolmethods:
        for atom_idx, atom_H_idx in atoms_CH+atoms_HH:
            rw_mol = RWMol(mother_mol)
            new_atom_idx = [rw_mol.AddAtom(Chem.Atom(6)),rw_mol.AddAtom(Chem.Atom(6))]
            rw_mol.AddBond(new_atom_idx[0],new_atom_idx[1],BondType.DOUBLE)
            rw_mol.AddBond(atom_idx,new_atom_idx[0],BondType.SINGLE)
            rw_mol.RemoveBond(atom_idx,atom_H_idx)
            rw_mol.RemoveAtom(atom_H_idx)
            [a.SetAtomMapNum(0) for a in rw_mol.GetAtoms()]
            try:
                smi_before = Chem.MolToSmiles(RemoveHs(rw_mol))
            except Exception as e:
                continue
            if is_EZ:
                for gen_smi in enumerate_EZ(smi_before):                       
                    generation_information.append((gen_smi,(smi,"vinyl",atom_idx)))
            else:
                generation_information.append((smi_before,(smi,"vinyl",atom_idx)))

    if "ethynyl" in evolmethods:
        for atom_idx, atom_H_idx in atoms_CH+atoms_HH:
            rw_mol = RWMol(mother_mol)
            new_atom_idx = [rw_mol.AddAtom(Chem.Atom(6)),rw_mol.AddAtom(Chem.Atom(6))]
            rw_mol.AddBond(new_atom_idx[0],new_atom_idx[1], BondType.TRIPLE)
            rw_mol.AddBond(atom_idx,new_atom_idx[0], BondType.SINGLE)
            rw_mol.RemoveBond(atom_idx,atom_H_idx)
            rw_mol.RemoveAtom(atom_H_idx)
            [a.SetAtomMapNum(0) for a in rw_mol.GetAtoms()]
            try:
                smi_before = Chem.MolToSmiles(RemoveHs(rw_mol))
            except Exception as e:
                continue
            if is_EZ:
                for gen_smi in enumerate_EZ(smi_before):                             
                    generation_information.append((gen_smi,(smi,"ethynyl",atom_idx)))
            else:
                generation_information.append((smi_before,(smi,"ethynyl",atom_idx)))

    # composite paths
    # search triple bonds
    triples = set()
    for bond in mother_mol.GetBonds():
        if str(bond.GetBondType()) == "TRIPLE":
            atom_a = bond.GetBeginAtom()
            atom_b = bond.GetEndAtom()
            bond_atoms = [atom_a.GetIdx(),atom_b.GetIdx()]
            bond_atoms.sort()
            triples.add(tuple(bond_atoms))
    
    BONDTYPE = (BondType.DOUBLE,BondType.SINGLE)
    
    if "annulate_pl2" in evolmethods or "annulate_pl4" in evolmethods:
        if "annulate_pl4" in evolmethods:
            cycle = 1
            connectable_pairs_pl4 = []
        if "annulate_pl2" in evolmethods:
            cycle = 3
            connectable_pairs_pl2 = []
        acceptable_atoms = set([x[0] for x in atoms_CH])
        for mother_atom in mother_mol.GetAtoms():
            radius = []
            seen_atoms = set()
            mother_atom_index = mother_atom.GetIdx()
            if not mother_atom_index in acceptable_atoms:
                continue
            radius.append([[mother_atom_index]])
            seen_atoms.add(mother_atom_index)
            for i in range(cycle):
                radius.append([])
                for atom_idxs in radius[-2]:
                    atom_index = atom_idxs[-1]
                    atom = mother_mol.GetAtomWithIdx(atom_index)
                    for bond in atom.GetBonds():
                        if bond.GetEndAtom().GetIdx() == atom_index:
                            next_atom = bond.GetBeginAtom()
                        else:
                            next_atom = bond.GetEndAtom()
                        next_atom_index = next_atom.GetIdx()
                        if next_atom_index in seen_atoms:
                            continue
                        seen_atoms.add(next_atom_index)
                        if (i%2 == 0 and (not bond.GetBondType() in {Chem.BondType.DOUBLE, Chem.BondType.AROMATIC})) or (i%2 == 1 and (not bond.GetBondType() in {Chem.BondType.SINGLE, Chem.BondType.AROMATIC})):
                            continue
                        if (not (atom_index,next_atom_index) in triples) and (not (next_atom_index,atom_index) in triples) and (((atom_index,next_atom_index) in ccbonds) or ((next_atom_index,atom_index) in ccbonds)):
                            radius[-1].append(atom_idxs+[next_atom_index])
            
            if "annulate_pl2" in evolmethods:
                for atom_idxs in radius[3]:
                    atom_index = atom_idxs[-1]
                    if (not atom_index in acceptable_atoms) or mother_atom_index > atom_index:
                        continue
                    rwmol = RWMol(Chem.MolFromSmiles(smi))
                    for i in range(3):
                        rwmol.GetBondBetweenAtoms(atom_idxs[i],atom_idxs[i+1]).SetBondType(BONDTYPE[i%2])
                    if not Kekulize_aromatic(rwmol) is None:
                        connectable_pairs_pl2.append(tuple([mother_atom_index,atom_index]))
            if "annulate_pl4" in evolmethods:
                for atom_idxs in radius[1]:
                    atom_index = atom_idxs[-1]
                    if (not atom_index in acceptable_atoms) or mother_atom_index > atom_index:
                        continue
                    rwmol = RWMol(Chem.MolFromSmiles(smi))
                    rwmol.GetBondBetweenAtoms(atom_idxs[0],atom_idxs[1]).SetBondType(BONDTYPE[0])
                    if not Kekulize_aromatic(rwmol) is None:               
                        connectable_pairs_pl4.append(tuple([mother_atom_index,atom_index]))

        atoms_CH_dict = dict()
        for (mother_atom_index,mother_atom_H_index) in atoms_CH:
            atoms_CH_dict[mother_atom_index] = mother_atom_H_index
        
        if "annulate_pl2" in evolmethods:
            for atom_pair in connectable_pairs_pl2:
                rw_mol = RWMol(mother_mol)
                new_bond_begin_atom_index = rw_mol.AddAtom(Chem.Atom(6))
                new_bond_end_atom_index = new_bond_begin_atom_index
                bond_types = (BondType.DOUBLE, BondType.SINGLE)
                for j in range(1):
                    prev_atom_idx = new_bond_end_atom_index
                    new_bond_end_atom_index = rw_mol.AddAtom(Chem.Atom(6))
                    rw_mol.AddBond(prev_atom_idx,new_bond_end_atom_index,bond_types[j%2])
                rw_mol.AddBond(atom_pair[0], new_bond_begin_atom_index, BondType.SINGLE)
                rw_mol.AddBond(new_bond_end_atom_index, atom_pair[1], BondType.SINGLE)
                
                H_pair = (atoms_CH_dict[atom_pair[0]],atoms_CH_dict[atom_pair[1]])
                rw_mol.RemoveBond(atom_pair[0],H_pair[0])
                rw_mol.RemoveBond(atom_pair[1],H_pair[1])
                if H_pair[0] > H_pair[1]:
                    rw_mol.RemoveAtom(H_pair[0])
                    rw_mol.RemoveAtom(H_pair[1])
                else:
                    rw_mol.RemoveAtom(H_pair[1])
                    rw_mol.RemoveAtom(H_pair[0])
                [a.SetAtomMapNum(0) for a in rw_mol.GetAtoms()]
                try:
                    smi_before = Chem.MolToSmiles(RemoveHs(rw_mol))
                except Exception as e:
                    continue
                if is_EZ:
                    for gen_smi in enumerate_EZ(smi_before):                              
                        generation_information.append((gen_smi,(smi,"annulate_pl2",atom_pair[0],atom_pair[1])))
                else:
                    generation_information.append((smi_before,(smi,"annulate_pl2",atom_pair[0],atom_pair[1])))
        
        if "annulate_pl4" in evolmethods:
            for atom_pair in connectable_pairs_pl4:
                rw_mol = RWMol(mother_mol)
                new_bond_begin_atom_index = rw_mol.AddAtom(Chem.Atom(6))
                new_bond_end_atom_index = new_bond_begin_atom_index
                bond_types = (BondType.DOUBLE, BondType.SINGLE)
                for j in range(3):
                    prev_atom_idx = new_bond_end_atom_index
                    new_bond_end_atom_index = rw_mol.AddAtom(Chem.Atom(6))
                    rw_mol.AddBond(prev_atom_idx,new_bond_end_atom_index,bond_types[j%2])
                rw_mol.AddBond(atom_pair[0], new_bond_begin_atom_index, BondType.SINGLE)
                rw_mol.AddBond(new_bond_end_atom_index, atom_pair[1], BondType.SINGLE)
                
                H_pair = (atoms_CH_dict[atom_pair[0]],atoms_CH_dict[atom_pair[1]])
                rw_mol.RemoveBond(atom_pair[0],H_pair[0])
                rw_mol.RemoveBond(atom_pair[1],H_pair[1])
                if H_pair[0] > H_pair[1]:
                    rw_mol.RemoveAtom(H_pair[0])
                    rw_mol.RemoveAtom(H_pair[1])
                else:
                    rw_mol.RemoveAtom(H_pair[1])
                    rw_mol.RemoveAtom(H_pair[0])
                [a.SetAtomMapNum(0) for a in rw_mol.GetAtoms()]
                try:
                    smi_before = Chem.MolToSmiles(RemoveHs(rw_mol))
                except Exception as e:
                    continue
                if is_EZ:
                    for gen_smi in enumerate_EZ(smi_before):                               
                        generation_information.append((gen_smi,(smi,"annulate_pl4",atom_pair[0],atom_pair[1])))
                else:
                    generation_information.append((smi_before,(smi,"annulate_pl4",atom_pair[0],atom_pair[1])))

    if "phenyl" in evolmethods:
        for atom_idx, atom_H_idx in atoms_CH+atoms_HH:
            rw_mol = RWMol(mother_mol)
            new_atom_idx = [rw_mol.AddAtom(Chem.Atom(6)) for _ in range(6)]
            [rw_mol.AddBond(new_atom_idx[i-1],new_atom_idx[i],BondType.AROMATIC) for i in range(6)]
            rw_mol.AddBond(atom_idx,new_atom_idx[0],BondType.SINGLE)
            rw_mol.RemoveBond(atom_idx,atom_H_idx)
            rw_mol.RemoveAtom(atom_H_idx)
            [a.SetAtomMapNum(0) for a in rw_mol.GetAtoms()]
            try:
                smi_before = Chem.MolToSmiles(RemoveHs(rw_mol))
            except Exception as e:
                continue
            if is_EZ:
                for gen_smi in enumerate_EZ(smi_before):                             
                    generation_information.append((gen_smi,(smi,"phenyl",atom_idx)))
            else:
                generation_information.append((smi_before,(smi,"phenyl",atom_idx)))

    return generation_information

if __name__ == "__main__":
    print(evol_single("C=C/C=C/C=C\C=C",evolmethods=("annulate_pl4")))





