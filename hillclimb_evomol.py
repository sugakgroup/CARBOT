from rdkit import Chem, RDLogger
from rdkit.Chem import rdchem

from rdkit.Chem import DataStructs
from rdkit.Chem.rdchem import RWMol

import time
import random

from evolve import *
from degen import *
from rediscovery import *
from utils import *

class ImproveFinder:
    def __init__(self,smi,evolmethods,is_degen,func,used_smi):
        self.smi = smi
        self.evolmethods = evolmethods
        self.is_degen = is_degen
        self.candidate = set()
        self.func = func
        self.used_smi = used_smi
        for gen_smi, _ in evol_single(self.smi,self.evolmethods,is_EZ=False):
            self.candidate.add((random.random(),gen_smi))
        if is_degen:
            for degen_smi, _ in degen_single(self.smi,self.evolmethods):
                self.candidate.add((random.random(),degen_smi))
        self.candidate = sorted(list(self.candidate))
        self.score_max = (-1,None)
        self.next_idx = 0

    def single(self):
        if self.next_idx >= len(self.candidate):
            return None, None
        smi = self.candidate[self.next_idx][1]
        if smi in self.used_smi:
            score = self.used_smi[smi]
        else:
            score = self.func(smi)
        self.next_idx += 1
        return score, smi

class HillClimbingEvoMol:
    def __init__(self,
                 func,
                 root="[H][H]",
                 is_composite_edge=False,
                 is_degen=True,
                 random_seed=42,
                 max_num_search_improver=50):
        time0 = time.time()
        self.root = root
        self.func = func
        random.seed(random_seed)
        self.log = []
        self.parents = dict()
        if is_degen:
            self.children = dict()
        self.used_smi = dict()
        self.used_smi[root] = func(root)
        self.replaced_smi = set([root])
        self.evolmethods = ("connect","vinyl","ethynyl","annulate_pl2","annulate_pl4","phenyl") if is_composite_edge else ("connect","vinyl","ethynyl")
        self.is_degen = is_degen
        self.pop = [(-1,"") for _ in range(1000)]
        self.pop[0] = ((func(root),root))
        self.pop_set = set([root])
        self.round = 1
        self.num_replace = 10
        self.num_eval_total = 1
        self.max_num_search_improver = max_num_search_improver
        print((self.round,self.num_eval_total,(func(root),root),time.time()-time0))
        self.log.append((self.round,self.num_eval_total,(func(root),root),time.time()-time0))
    
    def single(self):
        is_improved = False
        time0 = time.time()
        self.round += 1

        to_be_replaced = list(range(10))
        
        num_update = 0
        for score,smi in self.pop:
            if score < 0:
                break
            else:
                candidate = set()
                if not smi in self.parents.keys():
                    self.parents[smi] = tuple([gen_smi for gen_smi, _ in evol_single(smi,self.evolmethods,is_EZ=False)])
                parents = self.parents[smi]
                for gen_smi in parents:
                    if not gen_smi in self.replaced_smi:
                        candidate.add(gen_smi)
                if self.is_degen:
                    if not smi in self.children.keys():
                        self.children[smi] = tuple([degen_smi for degen_smi, _ in degen_single(smi,self.evolmethods)])
                    children = self.children[smi]
                    for degen_smi in children:
                        if not degen_smi in self.replaced_smi:
                            candidate.add(degen_smi)
                candidate = sorted(list(candidate))
                candidate = [(random.random(),cand) for cand in candidate]
                candidate.sort()
                
                max_score = -1
                max_smi = ""
                num_tries = 0
                
                for _, gen_smi in candidate:
                    if gen_smi in self.replaced_smi:
                        continue
                    if not gen_smi in self.used_smi.keys():
                        self.used_smi[gen_smi] = self.func(gen_smi)
                        self.num_eval_total += 1
                    score = self.used_smi[gen_smi]
                    if max_score < score:
                        max_score = score
                        max_smi = gen_smi
                    num_tries += 1
                    if num_tries >= self.max_num_search_improver:
                        break
                if max_score == -1:
                    continue

                for i in to_be_replaced:
                    if self.pop[-self.num_replace+i][0] <= max_score:
                        self.pop[-self.num_replace+i] = (max_score,max_smi)
                        self.replaced_smi.add(max_smi)
                        to_be_replaced.remove(i)
                        is_improved = True
                        num_update += 1
                        print((self.round,self.num_eval_total,(max_score,max_smi),time.time()-time0))
                        self.log.append((self.round,self.num_eval_total,(max_score,max_smi),time.time()-time0))
                        time0 = time.time()
                        break
                
                if len(to_be_replaced) == 0:
                    break
        
        if not is_improved:
            print("No longer improved!")
            print(f"Best:{self.pop[0]}")
            print(f"Worst:{self.pop[-1]}")
            return False

        self.pop.sort(reverse=True,key=lambda x: x[0])
        # print(self.pop)

        if self.pop[0][0] == 1.0:
            print("reproduced!")
            return False
        else:
            print(f"#{num_update} updated!")
            print(f"Best:{self.pop[0]}")
            print(f"Worst:{self.pop[-1]}")
        
        return True

    def multiple(self,num):
        for _ in range(num):
            result = self.single()
            if not result:
                return None

if __name__ == "__main__":
    RDLogger.DisableLog('rdApp.*')
    import pandas as pd
    rds = RediscoveryTanimoto("C12=C(C3=C(C4=C5C6=C78)C9=C2C%10=C(C%11=C%12%13)C%14=C9C4=C7C%15=C%14C%11=C%16C%17=C%15C8=C%18C%19=C%20%17)C(C(C3=C5C(C6=C%21%18)=C%22%23)=C%22C%24=C%25%26)=C%26C%27=C1C%10=C%12C%28=C%27C%25=C%29C%30=C%28C%13=C%16C%20=C%30C%31=C%19C%21=C%23C%24=C%31%29") # fullerene
    # rds = RediscoveryTanimoto("C12=CC=CC(C3=C(C4=CC=C5)C5=CC=C3)=C1C4=CC=C2") # perylene
    # rds = RediscoveryTanimoto("C12=CC=CC1=CC=CC=C2") # azulene
    # rds = Rediscovery("C1(C2=CC=C(C(C=C3)=CC=C3C4=CC=C(C5=CC=C6C=C5)C=C4)C=C2)=CC=C(C7=CC=C6C=C7)C=C1") # CPP6
    # rds = Rediscovery("c1c2/C=C\c3cc4cc5ccccc5cc4cc3/C=C\c2cc2cc3ccccc3cc21") #AnFLAP_EZ
    # rds = Rediscovery("c1c2C=Cc3cc4cc5ccccc5cc4cc3C=Cc2cc2cc3ccccc3cc21") #AnFLAP_noEZ
    # rds = Rediscovery("C12=CC=CC=C1C=CC3=C2C(C(C(C(C(C=CC=C4)=C4C=C5)=C5C=C6)=C6C=C7)=C7C=C8)=C8C=C3") # helicene8
    # rds = Rediscovery("C1(C=C2)=CC=C3C4=C1C5=C2C=CC(C=C6)=C5C4=C6C=C3") # corannulene
    # rds = Rediscovery("C1(C(C2=CC=C3C=C4)=C3C5=C4C=CC(C=C6)=C57)=C(C=C2)C=CC8=CC=C6C7=C18") # circulene7
    # rds = Rediscovery("C12=C3C=CC(C=C4)=C1C(C=C5)=C4C6=C5C=CC(C=C7)=C6C8=C7C=CC(C=C9)=C8C%10=C9C=CC%11=C%10C=CC(C=C%12)=C%11C%13=C%12C=CC(C=C3)=C2%13") # infinitene
    # rds = RediscoveryTanimoto("C1(C=CC2=CC=CC=C2)=CC=CC=C1") # stilbene
    # rds = Rediscovery("C1(/C=C\C2=CC=CC=C2)=CC=CC=C1") # stilbeneZ
    # rds = Rediscovery("C1(/C=C/C2=CC=CC=C2)=CC=CC=C1") # stilbeneE
    # rds = Rediscovery("C1(C2=CC3=CC4=CC5=CC6=CC7=CC8=CC9=CC%10=CC%11=CC%12=CC%13=C2)=CC3=CC4=CC5=CC6=CC7=CC8=CC9=CC%10=CC%11=CC%12=CC%13=C1") # cyclacene12
    # rds = Rediscovery("C1(C=C2C3=C4)=C4C(C=C(C=CC5=C6C=C7C(C(C=C(C=CC8=C9C=C%10C(C(C=C(C=C3)C2=C%11)=C%11C=C%10)=C8)C9=C%12)=C%12C=C7)=C5)C6=C%13)=C%13C=C1") # CNB6,6
    # rds = RediscoveryTanimoto("c1ccccc1") # benzene
    # rds = Rediscovery("C1(C#CC2=CC=C(C#CC3=CC=C(C#CC4=CC=C5C=C4)C=C3)C=C2)=CC=C(C#CC6=CC=C(C#CC7=CC=C(C#C5)C=C7)C=C6)C=C1") # CPPA6
    rds = Rediscovery("C1(C(C2=CC=CC=C2)=C3C=CC(C=C3)=C4C=CC(C=C4)=C(C5=CC=CC=C5)C6=CC=CC=C6)=CC=CC=C1") # Chichibabin
    # rds = Rediscovery("C1(C=C2C=CC(C=C2)=CC3=CC=C4C=C3)=CC=C(C=C5C=CC(C=C5)=CC6=CC=C(C=C6)C=C(C=C7)C=CC7=C4)C=C1") # CPPM6
    # rds = Rediscovery("C1(C=CC(C=C1)=CC2=CC=C(C=C(C=C3)C=CC3=C4)C=C2)=CC5=CC=C(C=C(C=C6)C=CC6=CC7=CC=C(C=C(C=C8)C=CC8=CC9=CC=C4C=C9)C=C7)C=C5") # CPPM8
    # rds = Rediscovery("C1(C=CC=C2)=C2C=CC(C#CC#CC(C=CC3=C4C=CC=C3)=C4C5=C6C=CC7=C5C=CC=C7)=C1C8=C(C#CC#CC9=C(C%10=C(C#CC#C6)C=CC%11=C%10C=CC=C%11)C(C=CC=C%12)=C%12C=C9)C=CC%13=C8C=CC=C%13") # MöbiusAnnulene
    # rds = Rediscovery("C12=C(C3=C(C4=C5C6=C78)C9=C2C%10=C(C%11=C%12%13)C%14=C9C4=C7C%15=C%14C%11=C%16C%17=C%15C8=C%18C%19=C%20%17)C(C(C3=C5C(C6=C%21%18)=C%22%23)=C%22C%24=C%25%26)=C%26C%27=C1C%10=C%12C%28=C%27C%25=C%29C%30=C%28C%13=C%16C%20=C%30C%31=C%19C%21=C%23C%24=C%31%29",weak=True) # fullerene
    # rds = Rediscovery("C1(C(C(C2=C(C=CC3=C2C(C(C4=C(C=CC5=C4C(C(C6=C(C=CC7=C6C(C(C=CC=C8)=C8C=C9)=C9C=C7)C=C%10)=C%10C=C%11)=C%11C=C5)C=C%12)=C%12C=C%13)=C%13C=C3)C=C%14)=C%14C=C%15)=C%15C=C%16)=C%16C=CC=C1") # helicene16
    # rds = RediscoveryTanimoto("C12=CC=CC=C1C=CC=C2") # naphthalene
    # rds = RediscoveryTanimoto("C1(C2=CC=CC=C2)=CC=CC=C1") # biphenyl
    # rds = RediscoveryTanimoto("C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C") # tetracosadodecaene
    rds = RediscoveryTanimoto("C=CC(C(C(C(C(C(C(C(C(C(C=C)=C)=C)=C)=C)=C)=C)=C)=C)=C)=C") # dendralene12

    for k in [True,False]:
        for j in [[10]]:
            num_eval_each_round = tuple(j)

            for i in range(10):
                hc = HillClimbingEvoMol(is_composite_edge=k,func=rds.score,random_seed=42+i,is_degen=True)
                hc.multiple(3000)
                df = pd.DataFrame([(a,b,c,d,e) for a,b,(c,d),e in hc.log], columns=["round","num_eval","score","smiles","time"])
                df.to_csv(f'data/rediscovery/hillclimbEvoMol_rediscovery-dendralene12_{"comp" if k else "unit"}_{i+42}.csv', index=False)


    
