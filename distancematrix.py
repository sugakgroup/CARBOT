from collections import defaultdict
from maxsubstructure import *
from utils import *
from rdkit import RDLogger

import csv
import pickle
from concurrent.futures import ProcessPoolExecutor

class DistanceMatrix:
    def __init__(self,original_smi=("c1ccccc1","c1ccccc1-c1ccccc1","C12=CC=CC=C1C=CC=C2")):
        self.original_smi = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in original_smi]
        self.D = [[-1 if i != j else 0 for j in range(len(self.original_smi))] for i in range(len(self.original_smi))]

    def find_mcs_flow_bactch(self,pairs):
        results = []
        for smi1,smi2,i,j in pairs:
            findmcs = FindMCS(smi1,smi2)
            ans = findmcs.find()
            print((smi1,smi2,ans[0],i,j,distance_std(smi1)+distance_std(smi2)-2*distance_std(ans[0])))
            results.append((i,j,distance_std(smi1)+distance_std(smi2)-2*distance_std(ans[0])))
        return results

    def construction(self):
        batches = [[] for _ in range(32)]
        n = 0
        for i in range(len(self.original_smi)):
            for j in range(i):
                batches[n%32].append((self.original_smi[i],self.original_smi[j],i,j))
                n += 1
        batches = [tuple(b) for b in batches]
        with ProcessPoolExecutor(max_workers=32) as executor:
            # results = list(executor.map(self.find_mcs_flow, [(self.nodes[i],self.nodes[j],i,j) for j in range(1,i)]))
            results = list(executor.map(self.find_mcs_flow_bactch, batches))
        for ans in results:
            for i,j,d in ans:
                self.D[i][j] = d
                self.D[j][i] = d
        print(self.D)

if __name__ == "__main__":
    RDLogger.DisableLog('rdApp.*')
    with open("data/tree/compounds.csv") as f:
        reader = csv.reader(f)
        compounds = [tuple(row) for row in reader]
    for n, smi, _ in compounds:
        print(n,smi,sum(num_unit_method(smi)))
    tc = DistanceMatrix(original_smi=tuple([i for _, i, _ in compounds]))
    tc.construction()
    with open(f'data/tree/distancematrix.pkl', "wb") as f:
        pickle.dump(tc, f)
