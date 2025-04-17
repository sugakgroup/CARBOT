from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, QED, rdchem
from rdkit.Chem import rdFingerprintGenerator as fpg
from rdkit.Chem import DataStructs
from rdkit.Chem.rdchem import RWMol
from SA_Score import sascorer
import time
import random

from evolve import *

SA_MEAN = 3.0509305172950056
SA_STD = 0.8381024508503762
LOGP_MEAN = 2.4713477640000012
LOGP_STD = 1.4449891243959943
CYCLE_MEAN = 0.0504
CYCLE_STD = 0.2902754553867757

def calc_Jscore(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        logP = Descriptors.MolLogP(mol)
        sascore = sascorer.calculateScore(mol)
        ri = mol.GetRingInfo()
        max_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
        ring_size_penalty = max(max_ring_size-6,0)
        # print(f"smi:{smi}, logp:{(logP-LOGP_MEAN)/LOGP_STD}, sa:{(sascore-SA_MEAN)/SA_STD}, cycle:{(ring_size_penalty-CYCLE_MEAN)/CYCLE_STD}")
        return logP - sascore - ring_size_penalty
    except:
        return None

def calc_Jscore_std(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        logP = Descriptors.MolLogP(mol)
        sascore = sascorer.calculateScore(mol)
        ri = mol.GetRingInfo()
        max_ring_size = max((len(r) for r in ri.AtomRings()), default=0)
        ring_size_penalty = max(max_ring_size-6,0)
        # print(f"smi:{smi}, logp:{(logP-LOGP_MEAN)/LOGP_STD}, sa:{(sascore-SA_MEAN)/SA_STD}, cycle:{(ring_size_penalty-CYCLE_MEAN)/CYCLE_STD}")
        return (logP-LOGP_MEAN)/LOGP_STD - (sascore-SA_MEAN)/SA_STD - (ring_size_penalty-CYCLE_MEAN)/CYCLE_STD
    except:
        return None
    
def calc_QED(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        qed = QED.qed(mol)
        return qed
    except:
        return None

def anneal_mol(mol):
    rwmol = Chem.RWMol(mol)
    for bond in mol.GetBonds():
        if bond.IsInRing() and (not bond.GetBondType() == Chem.BondType.TRIPLE):
            rwmol.GetBondWithIdx(bond.GetIdx()).SetBondType(Chem.BondType.AROMATIC)
    return Chem.MolFromSmiles(Chem.MolToSmiles(rwmol))

rdkfpgen = fpg.GetRDKitFPGenerator(
    minPath=1,
    maxPath=7,
    useHs=True,
    branchedPaths=True,
    useBondOrder=True,
    countSimulation=False,
    countBounds=None,
    fpSize=2024,
    numBitsPerFeature=2,
    atomInvariantsGenerator=None,
)

def calc_fp(smi):
    try:
        mol = Chem.AddHs(anneal_mol(Chem.MolFromSmiles(smi)))
    except:
        # mol = Chem.AddHs(Chem.MolFromSmiles(""))
        return None
    return rdkfpgen.GetCountFingerprint(mol)

def reproduce(smi,target_smi="C1(C2=C3C4=CC=CC3=CC=C2)=C(C4=CC=C5)C5=CC=C1"):
    try:
        target_smi = Chem.MolToSmiles(Chem.MolFromSmiles(target_smi))
        fp_smi = calc_fp(smi)
        fp_target_smi = calc_fp(target_smi)
        tanimoto = DataStructs.TanimotoSimilarity(fp_smi,fp_target_smi)
        return tanimoto
    except:
        return None

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
        if bond.GetBondType() != rdchem.BondType.TRIPLE:
            bond.SetBondType(rdchem.BondType.SINGLE)
    for atom in rwmol.GetAtoms():
        atom.SetIsAromatic(False)
    return rwmol

TARGET_SMI = "C12=CC=CC(C3=C(C4=CC=C5)C5=CC=C3)=C1C4=CC=C2"
# 'C12=C(C3=C(C4=C5C6=C78)C9=C2C%10=C(C%11=C%12%13)C%14=C9C4=C7C%15=C%14C%11=C%16C%17=C%15C8=C%18C%19=C%20%17)C(C(C3=C5C(C6=C%21%18)=C%22%23)=C%22C%24=C%25%26)=C%26C%27=C1C%10=C%12C%28=C%27C%25=C%29C%30=C%28C%13=C%16C%20=C%30C%31=C%19C%21=C%23C%24=C%31%29'
TARGET_MOL = molfromsmiles_skeleton(TARGET_SMI)

def reproduce2(smi):
    try:
        target_mol = TARGET_MOL 
        mol = molfromsmiles_skeleton(smi)
        if target_mol.HasSubstructMatch(mol):
            return (sum(num_unit_method(smi))+1)/(sum(num_unit_method(TARGET_SMI))+1)
        else:
            return 0
    except:
        return 0

class HillClimbing:
    def __init__(self,
                 root="[H][H]",
                 func=calc_Jscore_std,
                 is_composite_edge=False,
                 num_eval_each_round=(2,20,200),
                 random_seed=42):
        time0 = time.time()
        self.root = root
        self.func = func
        self.num_eval_each_round = num_eval_each_round
        random.seed(random_seed)
        self.log = []
        self.used_smi = set([root])
        self.evolmethods = ("connect","ethylene","acetylene","annulate_pl2","annulate_pl4","phenyl") if is_composite_edge else ("connect","ethylene","acetylene")
        self.score_list = []

        self.score_list.append((func(root),root))
        self.round = 1
        self.log.append((self.round,(func(root),root),1,time.time()-time0))
    
    def single(self):
        time0 = time.time()
        self.round += 1

        num_gen_smi = 0
        while num_gen_smi == 0:
            candidate = set()
            for gen_smi, _ in evol_single(self.score_list[-1][1],self.evolmethods):
                if gen_smi in self.used_smi:
                    continue
                candidate.add((random.random(),gen_smi))
            
            candidate = sorted(list(candidate))

            num_gen_smi = len(candidate)
            if num_gen_smi == 0:
                self.score_list.pop(-1)

        # select the best structure
        num_eval = 0
        i = 0
        num_eval_each_round = random.sample(self.num_eval_each_round,1)[0]
        while i < num_gen_smi and num_eval < num_eval_each_round:
            smi = candidate[i][1]
            i += 1
            if smi in self.used_smi:
                continue
            else:
                score = self.func(smi)
                num_eval += 1
                if score is None:
                    continue
                self.used_smi.add(smi)
                self.score_list.append((score,smi))
        self.score_list.sort()
        
        # if self.score_list[-1][1].lower().count("c") > 38:
        #     print("exceed heavy atom limit")
        #     return False

        print(self.round,self.score_list[-1],num_gen_smi,time.time()-time0)
        self.log.append((self.round,self.score_list[-1],num_gen_smi,time.time()-time0))

        if self.score_list[-1][0] == 1.0:
            print("reproduced!")
            return False
        
        return True
      

    def multiple(self,num):
        for _ in range(num):
            result = self.single()
            if not result:
                return None

if __name__ == "__main__":
    RDLogger.DisableLog('rdApp.*')
    import pandas as pd

    for k in [True,False]:
        for j in [[10],[100],[1000]]:
            num_eval_each_round = tuple(j)

            for i in range(10):
                hc = HillClimbing(is_composite_edge=k,func=reproduce2,num_eval_each_round=num_eval_each_round,random_seed=42+i)
                hc.multiple(1000)
                df = pd.DataFrame([(a,b,c,d,e) for a,(b,c),d,e in hc.log], columns=["round","score","smiles","num_children","time"])
                df.to_csv(f'data/hillclimb_reproduce2-perylene_{"comp" if k else "unit"}_{num_eval_each_round}_{i+42}.csv', index=False)

    