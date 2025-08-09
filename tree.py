from collections import defaultdict
from maxsubstructure import *
from utils import *
from rdkit import RDLogger

import csv
import pickle
from concurrent.futures import ProcessPoolExecutor

class TreeConstruction:
    def __init__(self,original_smi=("c1ccccc1","c1ccccc1-c1ccccc1","C12=CC=CC=C1C=CC=C2")):
        self.original_smi = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in original_smi]
        self.nodes = ["[H][H]"]+self.original_smi
        self.edges = []
        self.parents = defaultdict(set)
    
    def find_mcs_flow(self,pair):
        smi1,smi2,i,j = pair
        findmcs = FindMCS(smi1,smi2,is_EZ=False)
        ans = findmcs.find()
        print(i,j,"finish")
        return ans

    def find_mcs_flow_bactch(self,pairs):
        result = dict()
        for smi1,smi2,i,j in pairs:
            findmcs = FindMCS(smi1,smi2,is_EZ=False)
            ans = findmcs.find()
            # print(i,j,ans)
            result[(i,j)] = ans
        return result

    def nodes_construction(self):
        i = 2
        while i < len(self.nodes):
            print(i,self.nodes[i])
            batches = [[] for _ in range(30)]
            for j in range(1,i):
                batches[j%30].append((self.nodes[i],self.nodes[j],i,j))
            batches = [tuple(b) for b in batches]
            with ProcessPoolExecutor(max_workers=30) as executor:
                # results = list(executor.map(self.find_mcs_flow, [(self.nodes[i],self.nodes[j],i,j) for j in range(1,i)]))
                results = list(executor.map(self.find_mcs_flow_bactch, batches))
            
            for ans in results:
                for (i,j), smis in ans.items():
                    # print(i,j,smi)
                    for smi in smis:
            # for j in range(1,i):
            #     for smi in results[j-1]:
                # for smi in self.find_mcs_flow((self.nodes[i],self.nodes[j])):
                        if not smi in self.nodes:
                            self.nodes.append(smi)
                        if smi != self.nodes[i]:
                            self.parents[self.nodes[i]].add((sum(num_unit_method(smi)),smi))
                        if smi != self.nodes[j]:
                            self.parents[self.nodes[j]].add((sum(num_unit_method(smi)),smi))
            i += 1
            print(f"total #{len(self.nodes)} nodes generated now (i = {i})")
        for k in self.nodes:
            if k != "[H][H]":
                self.parents[k].add((0,"[H][H]"))
        for k in self.parents.keys():
            self.parents[k] = sorted(list(self.parents[k]),reverse=True)
                    
    def edges_construction(self):
        for k,v in self.parents.items():
            child_gen = sum(num_unit_method(k))
            for i,(l,smi2) in enumerate(v):
                f = True
                for (_,smi1) in v[:i]:
                    rds = Rediscovery(smi1)
                    f = (rds.score(smi2) == 0)
                    if f == False:
                        break
                if f:
                   self.edges.append((smi2,k,child_gen-l)) 

if __name__ == "__main__":
    RDLogger.DisableLog('rdApp.*')
    with open("data/tree/compounds.csv") as f:
        reader = csv.reader(f)
        compounds = [tuple(row) for row in reader]
    for n, smi, _ in compounds:
        print(n,smi,sum(num_unit_method(smi)))
    tc = TreeConstruction(original_smi=tuple([i for _, i, _ in compounds]))
    tc.nodes_construction()
    tc.edges_construction()
    with open(f'data/tree/tree.pkl', "wb") as f:
        pickle.dump(tc, f)
