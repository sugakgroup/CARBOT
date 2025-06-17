from rdkit import RDLogger
from rdkit.Chem import rdFMCS
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager, Lock
import time

from rediscovery import *
from evolve import *
from utils import *

class FindMCS:
    def __init__(self,smi1,smi2,is_EZ=False):
        self.rds1 = Rediscovery(smi1,is_EZ=is_EZ)
        self.rds2 = Rediscovery(smi2,is_EZ=is_EZ)
        self.is_EZ = is_EZ
        # print(f"[set] {smi1} and {smi2} are set")

        self.substructures = [set(["[H][H]"])]

    
    def run_task(self,item):
        smi,checked,lock = item
        sbsts = set()
        with lock:
            gen_smis = sorted(list(set([gen_smi for gen_smi, _ in evol_single(smi,evolmethods=("connect","vinyl"),is_EZ=self.is_EZ) if not gen_smi in checked])))
            for gen_smi in gen_smis:
                checked[gen_smi] = True
        for gen_smi in gen_smis:
            if self.rds1.score(gen_smi) != 0 and self.rds2.score(gen_smi) != 0:
                sbsts.add(gen_smi)
        return sbsts
    
    def run_single(self):
        print(f"[start] for generation #{len(self.substructures)} from structures: #{len(self.substructures[-1])}")
        self.substructures.append(set())

        with Manager() as manager:
            checked = manager.dict()
            lock = manager.Lock()

            with ProcessPoolExecutor(max_workers=6) as executor:
                results = list(executor.map(self.run_task, [(smi,checked,lock) for smi in self.substructures[-2]]))

            for result in results:
                for gen_smi in result:
                    self.substructures[-1].add(gen_smi)

        if len(self.substructures[-1]) == 0:
            return False
        else:
            return True
    
    def find(self):
        if sum(num_unit_method(self.rds1.target_smi)) > sum(num_unit_method(self.rds2.target_smi)) and self.rds1.score(self.rds2.target_smi) != 0:
            print(f"MCS found between {self.rds1.target_smi} and {self.rds2.target_smi}: {set([self.rds2.target_smi])} in gen #{sum(num_unit_method(self.rds2.target_smi))}!")
            return [self.rds2.target_smi]
        elif sum(num_unit_method(self.rds1.target_smi)) < sum(num_unit_method(self.rds2.target_smi)) and self.rds2.score(self.rds1.target_smi) != 0:
            print(f"MCS found between {self.rds1.target_smi} and {self.rds2.target_smi}: {set([self.rds1.target_smi])} in gen #{sum(num_unit_method(self.rds1.target_smi))}!")
            return [self.rds1.target_smi]
        
        # # use of rdFMCS
        # res = rdFMCS.FindMCS([self.rds1.target_mol,self.rds2.target_mol])
        # mcs_smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(Kekulize_aromatic(Chem.MolFromSmarts(res.smartsString)))))
        # print(mcs_smi)
        # if self.rds1.score(mcs_smi) != 0 and self.rds2.score(mcs_smi) != 0:
        #     return set(mcs_smi)

        # brute-force
        flag = True
        while flag:
            flag = self.run_single()
        print(f"MCS found between {self.rds1.target_smi} and {self.rds2.target_smi}: {self.substructures[-2]} in gen #{len(self.substructures)-1}")
        return sorted(list(self.substructures[-2]))

if __name__ == "__main__":
    RDLogger.DisableLog('rdApp.*')
    smi1 = "c1(ccc2ccc3ccc4ccc5cc6)c7c2c3c4c5c7c6cc1" # perylene
    smi2 = "c1(c2ccc(c3ccc(c4ccccc4)cc3)cc2)ccccc1" # pyrene
    findmcs = FindMCS(smi1,smi2)
    findmcs.find()
        