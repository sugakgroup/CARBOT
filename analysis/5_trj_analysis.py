import pandas as pd

import csv
from pathlib import Path
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

proj_root = (Path.cwd().parent).resolve()
sys.path.insert(0, str(proj_root))

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

from src.evolve import evol_single
from src.rediscovery import *
from src.utils import *

def worker(batch,i):
    for idx,name,target_smi,seed in batch:
        rds = RediscoverySubstructure(target_smi=target_smi)
        for op_name in ["cve","cve24p"]:
            df = pd.read_csv(f"../data/2_rediscovery/hillclimb/substructure_{idx}_{name}_{op_name}_{seed}.csv")
            with open(f'../data/4_trajectory/{idx}_{name}_{op_name}_{seed}_trajanalysis.csv', "w", encoding='utf-8-sig', newline='') as f_out:
                writer = csv.writer(f_out)
                writer.writerow(["SMILES","Step","NumChild","NumOKChild","AcceptRatio","num_connect","num_vinyl","num_ethynyl"])
                for smi,score in zip(df["SMILES"],df["Score"]):
                    if score is not None and score > -0.5:
                        steps = num_unit_method(smi)
                        children = set()
                        ok_children = set()
                        operations = ("connect","vinyl","ethynyl") if op_name == "cve" else ("connect","vinyl","ethynyl","annulate_pl2","annulate_pl4","phenyl")
                        for _,gen_smis in evol_single(smi,operations=operations,is_EZ=False).items():
                            for gen_smi in gen_smis:
                                children.add(gen_smi)
                                if rds.score(gen_smi) > -0.5:
                                    ok_children.add(gen_smi)
                        writer.writerow([smi,sum(steps),len(children),len(ok_children),len(ok_children)/len(children) if len(children) > 0 else 0,steps[0],steps[1],steps[2]])
        print(idx,name,seed,"done!")
    return i

def main():
    with open('../data/target_molecules.csv', encoding='utf-8-sig', newline='') as f:
        items = f.read().split("\r\n")
        items = [tuple(item.split(",")) for item in items[:25]]

    jobs = []

    for idx,name,target_smi in items:
        for seed in range(42,52):
            jobs.append((idx,name,target_smi,seed))
    
    threads = 12
    with ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(worker,batch=[job for j,job in enumerate(jobs) if j % threads == i],i=i) for i in range(threads)]
        for future in as_completed(futures):
            print(f"{future.result()} done!")

if __name__ == "__main__":
    RDLogger.DisableLog('rdApp.*')
    main()
