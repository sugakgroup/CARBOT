from rdkit import Chem, RDLogger
from rdkit.Chem import rdchem

from rdkit.Chem import DataStructs
from rdkit.Chem.rdchem import RWMol

import time
import random

from evolve import *
from rediscovery import *
from utils import *

class HillClimbing:
    def __init__(self,
                 func,
                 root="[H][H]",
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
    # rds = Rediscovery("C12=C(C3=C(C4=C5C6=C78)C9=C2C%10=C(C%11=C%12%13)C%14=C9C4=C7C%15=C%14C%11=C%16C%17=C%15C8=C%18C%19=C%20%17)C(C(C3=C5C(C6=C%21%18)=C%22%23)=C%22C%24=C%25%26)=C%26C%27=C1C%10=C%12C%28=C%27C%25=C%29C%30=C%28C%13=C%16C%20=C%30C%31=C%19C%21=C%23C%24=C%31%29",weak=True) # fullerene
    # rds = Rediscovery("C1(C(C(C2=C(C=CC3=C2C(C(C4=C(C=CC5=C4C(C(C6=C(C=CC7=C6C(C(C=CC=C8)=C8C=C9)=C9C=C7)C=C%10)=C%10C=C%11)=C%11C=C5)C=C%12)=C%12C=C%13)=C%13C=C3)C=C%14)=C%14C=C%15)=C%15C=C%16)=C%16C=CC=C1") # helicene16
    # rds = Rediscovery("C12=CC=CC=C1C=CC=C2") # naphthalene
    # rds = Rediscovery("C1(C2=CC=CC=C2)=CC=CC=C1") # biphenyl
    # rds = Rediscovery("C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C/C=C") # tetracosadodecaene
    # rds = Rediscovery("C=CC(C(C(C(C(C(C(C(C(C(C=C)=C)=C)=C)=C)=C)=C)=C)=C)=C)=C") # dendralene12
    # rds = Rediscovery("C12=CC(C=CC3=C4C=C5C(C(C=C6C(C=C7)=C8)=C8C=C5)=C3)=C4C=C1C9=CC(C(C%10=C%11)=CC%12=C%11C=CC%13=CC%14=C(C=C%13%12)C=CC%15=C%14C=C(C%16=CC(C(C%17=C%18)=CC%19=C%18C=CC%20=CC%21=C(C=C%20%19)C=CC%22=C%21C=C(C%23=CC6=C7C=C%23C=C%24)C%24=C%22)=C(C=C%17)C=C%16C=C%25)C%25=C%15)=C(C=C%10)C=C9C=C2") # MöbiusCNB3
    # rds = Rediscovery("C12=C(C3=C4C=CC=C35)C(C6=C5C=CC=C67)=C(C8=C7C(C9=C%10C(C%11=C(C%12=C%10C=CC=C%12C%13=C%14C%15=CC=C%13)C%14=C(C%16=C%15C=CC=C%16%17)C(C%18=C%17C=CC=C%18%19)=C%11C%20=C%19C=CC=C%20%21)=C%21C%22=C%23C(C%24=C(C%25=C%23C=CC=C%25%26)C(C%27=C%26C=CC=C%27%28)=C(C%29=C%28C=CC=C%29%30)C(C%31=C%30C=CC=C%31%32)=C%24C%33=C%32C=CC=C%33%34)=C%34C%35=C%229)=C(C%35=C%36%37)C(C%38=C%36C%39=C%40C(C%41=C%42C=CC=C%41C%43=C%40C%38=CC=C%43)=C(C%44=C%45C=CC=C%44%42)C%46=C%39C%47=C%37C=CC=C%47C%48=C%46C%45=CC=C%48)=C8%49)C(C%50=C%49C=CC=C%50%51)=C1C%52=C%51C=CC=C%52C%53=C2C4=CC=C%53") # helicalNG
    # rds = Rediscovery("C1(C=C2C(C3=CC(C=C4)=C(C=C3C=C2)C(C4=C5)=CC6=C5C7=CC(C=C8)=C(C=C7C=C6)C(C8=C9)=CC%10=C9C(C=C(C%11=C%12)C=CC(C%11=C%13)=CC%14=C%13C=CC%15=CC%16=C(C=C%15%14)C=CC(C%16=C%17)=CC%18=C%17C=CC%19=CC%20=C(C=C%19%18)C=CC%21=C%20C=C%22C%23=C%21)=C%12C=C%10)=C%24)=C%24C=CC%25=C1C=C%26C(C(C=C(C=CC%27=C%28C=C%29C(C(C=C(C=CC(C%30=C%31)=CC%32=C%31C=CC%33=C%32C=C%34C(C(C=C(C=CC%35=C%36C=C%37C(C(C=C(C=CC(C%38=C%39)=CC%40=C%39C=CC%41=C%40C=C%42C(C(C=C(C=CC%43=C%44C=C%45C(C(C=C(C=C%22)C%23=C%46)=C%46C=C%45)=C%43)C%44=C%47)=C%47C=C%42)=C%41)C%38=C%48)=C%48C=C%37)=C%35)C%36=C%49)=C%49C=C%34)=C%33)C%30=C%50)=C%50C=C%29)=C%27)C%28=C%51)=C%51C=C%26)=C%25") # MöbiusCNB1
    # rds = Rediscovery("C12=CC(C=CC3=C4C=CC5=C3C=CC(C6=CC(C=CC7=C8C=CC9=C7C=CC(C%10=CC(C=CC%11=C%12C=CC%13=C%11C=CC(C%14=CC(C=CC%15=C%16C=CC%17=C%15C=CC2=C%17)=C%16C=C%14)=C%13)=C%12C=C%10)=C9)=C8C=C6)=C5)=C4C=C1") # cyclochrysenylene4
    # rds = Rediscovery("C12=C(C#CC3=CC=CC=C3C#C4)C(C#CC5=C4C=CC=C5)=C(C#CC6=CC=CC=C6C#CC(C=CC=C7)=C7C#C8)C8=C1C#CC9=CC=CC=C9C#CC(C=CC=C%10)=C%10C#C2") # graphyne
    rds = Rediscovery("C1(C=C2)=C(C=C3)C=CC3=C4C=CC(C=C4)=C(C=C5)C=CC5=C(C=C6)C=CC6=C(C=C7)C=CC7=C2C=C1") # CPP6quinoid


    is_EZ = False
    rds.is_EZ = is_EZ

    for k in [True,False]:
        for j in [[10]]:
            num_eval_each_round = tuple(j)

            for i in range(10):
                hc = HillClimbing(is_composite_edge=k,func=rds.score,random_seed=42+i,is_EZ=is_EZ)
                hc.multiple(10000)
                df = pd.DataFrame([(a,b,c,d,e,f) for a,b,(c,d),e,f in hc.log], columns=["round","num_eval","score","smiles","num_children","time"])
                df.to_csv(f'data/rediscovery/hillclimb_rediscovery-CPP6quinoid{"_isEZ" if is_EZ else ""}_{"comp" if k else "unit"}_{i+42}.csv', index=False)


    
