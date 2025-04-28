from rdkit import Chem, RDLogger
from rdkit.Chem import rdchem

from rdkit.Chem import DataStructs
from rdkit.Chem.rdchem import RWMol

import time
import random

from evolve import *
from rediscovery import *
from utils import *

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
                 func=reproduce2,
                 is_composite_edge=False,
                 is_EZ=False,
                 random_seed=42):
        time0 = time.time()
        self.root = root
        self.func = func
        random.seed(random_seed)
        self.log = []
        self.used_smi = dict()
        self.used_smi[root] = func(root)
        self.evolmethods = ("connect","vinyl","ethynyl","annulate_pl2","annulate_pl4","phenyl") if is_composite_edge else ("connect","vinyl","ethynyl")
        self.is_EZ = is_EZ
        self.history = []
        self.history.append((func(root),root))
        self.round = 1
        self.log.append((self.round,1,(func(root),root),1,time.time()-time0))
    
    def single(self):
        time0 = time.time()
        self.round += 1

        num_gen_smi = 0
        while num_gen_smi == 0:
            candidate = set()
            for gen_smi, _ in evol_single(self.history[-1][1],self.evolmethods,is_EZ=self.is_EZ):
                if gen_smi in self.used_smi.keys():
                    continue
                candidate.add((random.random(),gen_smi))
            
            candidate = sorted(list(candidate))

            num_gen_smi = len(candidate)
            if num_gen_smi == 0:
                self.history.pop(-1)

        # select the best structure
        num_eval = 0
        i = 0
        while i < num_gen_smi:
            smi = candidate[i][1]
            i += 1
            if smi in self.used_smi.keys():
                continue
            else:
                score = self.func(smi)
                num_eval += 1
                if score is None:
                    continue
                self.used_smi[smi] = score
                print(self.round,self.log[-1][1]+1,(score,smi),num_gen_smi-i+1,time.time()-time0)
                self.log.append((self.round,self.log[-1][1]+1,(score,smi),num_gen_smi-i+1,time.time()-time0))
                if score > self.history[-1][0]:
                    self.history.append((score,smi))
                    break

        if self.history[-1][0] == 1.0:
            print("reproduced!")
            return False
        else:
            print("updated!")
        
        return True
      

    def multiple(self,num):
        for _ in range(num):
            result = self.single()
            if not result:
                return None

if __name__ == "__main__":
    RDLogger.DisableLog('rdApp.*')
    import pandas as pd
    # rds = Rediscovery("C12=C(C3=C(C4=C5C6=C78)C9=C2C%10=C(C%11=C%12%13)C%14=C9C4=C7C%15=C%14C%11=C%16C%17=C%15C8=C%18C%19=C%20%17)C(C(C3=C5C(C6=C%21%18)=C%22%23)=C%22C%24=C%25%26)=C%26C%27=C1C%10=C%12C%28=C%27C%25=C%29C%30=C%28C%13=C%16C%20=C%30C%31=C%19C%21=C%23C%24=C%31%29") # fullerene
    # rds = Rediscovery("C12=CC=CC(C3=C(C4=CC=C5)C5=CC=C3)=C1C4=CC=C2") # perylene
    # rds = Rediscovery("C12=CC=CC1=CC=CC=C2") # azulene
    # rds = Rediscovery("C1(C2=CC=C(C(C=C3)=CC=C3C4=CC=C(C5=CC=C6C=C5)C=C4)C=C2)=CC=C(C7=CC=C6C=C7)C=C1") # CPP6
    # rds = Rediscovery("c1c2/C=C\c3cc4cc5ccccc5cc4cc3/C=C\c2cc2cc3ccccc3cc21") #AnFLAP_EZ
    # rds = Rediscovery("c1c2C=Cc3cc4cc5ccccc5cc4cc3C=Cc2cc2cc3ccccc3cc21") #AnFLAP_noEZ
    # rds = Rediscovery("C12=CC=CC=C1C=CC3=C2C(C(C(C(C(C=CC=C4)=C4C=C5)=C5C=C6)=C6C=C7)=C7C=C8)=C8C=C3") # helicene8
    # rds = Rediscovery("C1(C=C2)=CC=C3C4=C1C5=C2C=CC(C=C6)=C5C4=C6C=C3") # corannulene
    # rds = Rediscovery("C1(C(C2=CC=C3C=C4)=C3C5=C4C=CC(C=C6)=C57)=C(C=C2)C=CC8=CC=C6C7=C18") # circulene7
    # rds = Rediscovery("C12=C3C=CC(C=C4)=C1C(C=C5)=C4C6=C5C=CC(C=C7)=C6C8=C7C=CC(C=C9)=C8C%10=C9C=CC%11=C%10C=CC(C=C%12)=C%11C%13=C%12C=CC(C=C3)=C2%13") # infinitene
    # rds = Rediscovery("C1(C=CC2=CC=CC=C2)=CC=CC=C1") # stilbene
    # rds = Rediscovery("C1(/C=C\C2=CC=CC=C2)=CC=CC=C1") # stilbeneZ
    # rds = Rediscovery("C1(/C=C/C2=CC=CC=C2)=CC=CC=C1") # stilbeneE
    # rds = Rediscovery("C1(C2=CC3=CC4=CC5=CC6=CC7=CC8=CC9=CC%10=CC%11=CC%12=CC%13=C2)=CC3=CC4=CC5=CC6=CC7=CC8=CC9=CC%10=CC%11=CC%12=CC%13=C1") # cyclacene12
    # rds = Rediscovery("C1(C=C2C3=C4)=C4C(C=C(C=CC5=C6C=C7C(C(C=C(C=CC8=C9C=C%10C(C(C=C(C=C3)C2=C%11)=C%11C=C%10)=C8)C9=C%12)=C%12C=C7)=C5)C6=C%13)=C%13C=C1") # CNB6,6
    # rds = Rediscovery("c1ccccc1") # benzene
    # rds = Rediscovery("C1(C#CC2=CC=C(C#CC3=CC=C(C#CC4=CC=C5C=C4)C=C3)C=C2)=CC=C(C#CC6=CC=C(C#CC7=CC=C(C#C5)C=C7)C=C6)C=C1") # CPPA6
    # rds = Rediscovery("C1(C(C2=CC=CC=C2)=C3C=CC(C=C3)=C4C=CC(C=C4)=C(C5=CC=CC=C5)C6=CC=CC=C6)=CC=CC=C1") # Chichibabin
    # rds = Rediscovery("C1(C=C2C=CC(C=C2)=CC3=CC=C4C=C3)=CC=C(C=C5C=CC(C=C5)=CC6=CC=C(C=C6)C=C(C=C7)C=CC7=C4)C=C1") # CPPM6
    # rds = Rediscovery("C1(C=CC(C=C1)=CC2=CC=C(C=C(C=C3)C=CC3=C4)C=C2)=CC5=CC=C(C=C(C=C6)C=CC6=CC7=CC=C(C=C(C=C8)C=CC8=CC9=CC=C4C=C9)C=C7)C=C5") # CPPM8
    # rds = Rediscovery("C1(C=CC=C2)=C2C=CC(C#CC#CC(C=CC3=C4C=CC=C3)=C4C5=C6C=CC7=C5C=CC=C7)=C1C8=C(C#CC#CC9=C(C%10=C(C#CC#C6)C=CC%11=C%10C=CC=C%11)C(C=CC=C%12)=C%12C=C9)C=CC%13=C8C=CC=C%13") # MöbiusAnnulene
    rds = Rediscovery("C12=C(C3=C(C4=C5C6=C78)C9=C2C%10=C(C%11=C%12%13)C%14=C9C4=C7C%15=C%14C%11=C%16C%17=C%15C8=C%18C%19=C%20%17)C(C(C3=C5C(C6=C%21%18)=C%22%23)=C%22C%24=C%25%26)=C%26C%27=C1C%10=C%12C%28=C%27C%25=C%29C%30=C%28C%13=C%16C%20=C%30C%31=C%19C%21=C%23C%24=C%31%29",weak=True) # fullerene


    for k in [True,False]:
        for j in [[10]]:
            num_eval_each_round = tuple(j)

            for i in range(10):
                hc = HillClimbing(is_composite_edge=k,func=rds.score,random_seed=42+i,is_EZ=False)
                hc.multiple(10000)
                df = pd.DataFrame([(a,b,c,d,e,f) for a,b,(c,d),e,f in hc.log], columns=["round","num_eval","score","smiles","num_children","time"])
                df.to_csv(f'data/rediscovery/hillclimb_rediscovery-fullerene_weak_{"comp" if k else "unit"}_{i+42}.csv', index=False)

    
