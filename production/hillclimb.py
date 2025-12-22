from rdkit import RDLogger

from concurrent.futures import ProcessPoolExecutor, as_completed

import csv
import random
import sys
import time
from pathlib import Path

proj_root = (Path.cwd().parent).resolve()
sys.path.insert(0, str(proj_root))


from src.evolve import *
from src.rediscovery import *
from src.utils import *

class HillClimbing:
    def __init__(self,
                 title,
                 score_func,
                 root="[H][H]",
                 operations=("connect","vinyl","ethynyl"),
                 random_seed=42,
                 base_dir=None,
                 save_dir="2_rediscovery"):
        if base_dir is None:
            base_dir = Path.cwd().parent
        base_dir = Path(base_dir)
        self.out_dir = base_dir / "data" / save_dir / "hillclimb"
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.title = title

        self.time0 = time.time()
        self.root = root
        self.score_func = score_func
        self.blacklist = {root}
        random.seed(random_seed)
        self.operations = operations
        self.round = 0
        self.total_eval = 1
        score = score_func(root)
        self.history = [(score, root)]

        with open(self.out_dir / f"{self.title}.out", "w", encoding="utf-8") as f_log:
            f_log.write(f"[{time.time()-self.time0}] Start running HillClimbing from {root}...\n")

        with open(self.out_dir / f"{self.title}.csv", "w", newline="", encoding="utf-8") as f:
            writer_nodes = csv.writer(f)
            writer_nodes.writerow(["TotalEval","Round","NumEvalRound","SMILES","Score","TimeElapsed"])
            writer_nodes.writerow([self.total_eval,self.round,"1_1",root,score,time.time()-self.time0])

    def single(self):
        self.round += 1
        candidate = set()
        for _, smis in evol_single(self.history[-1][1],self.operations).items():
            for smi in smis:
                if smi in self.blacklist:
                    continue
                candidate.add(smi)
        candidate = sorted(list(candidate))
        random.shuffle(candidate)

        num_candidate = len(candidate)

        with open(self.out_dir / f"{self.title}.csv", "a", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["",self.round,f"0_{num_candidate}", "CandidateEnumeration", "", time.time()-self.time0])

        # select the best structure
        is_updated = False
        for i, smi in enumerate(candidate):
            score = self.score_func(smi)
            self.total_eval += 1
            with open(self.out_dir / f"{self.title}.csv", "a", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow([self.total_eval,self.round,f"{i+1}_{num_candidate}", smi, score, time.time()-self.time0])
            if score >= self.history[-1][0]:
                self.history.append((score,smi))
                is_updated = True
                self.blacklist.add(smi)
                break

        if not is_updated:
            self.history.pop(-1)
            if len(self.history) == 0:
                with open(self.out_dir / f"{self.title}.out", "a", encoding="utf-8") as f_log:
                    f_log.write(f"[{time.time()-self.time0}] No more candidate found in round {self.round}. Terminating...\n")
                return False
            else:
                with open(self.out_dir / f"{self.title}.out", "a", encoding="utf-8") as f_log:
                    f_log.write(f"[{time.time()-self.time0}] No improvement found in round {self.round}. Roll back to {self.history[-1]}\n")
        
        else:
            with open(self.out_dir / f"{self.title}.out", "a", encoding="utf-8") as f_log:
                f_log.write(f"[{time.time()-self.time0}] Updated in round {self.round}! max:{self.history[-1]}\n")
        if self.history[-1][0] == 1.0:
            with open(self.out_dir / f"{self.title}.out", "a", encoding="utf-8") as f_log:
                f_log.write(f"[{time.time()-self.time0}] Reproduced in round {self.round}!!\n")
            return False
        return True
      
    def multiple(self,num):
        for _ in range(num):
            result = self.single()
            if not result:
                return None

def run_task_pool(task_pool=[], method="substructure", threads=16):
    task_pool = [tuple([item for j,item in enumerate(task_pool) if j%threads == i]) for i in range(threads)]

    with ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(run_task_batch, task_batch=task_batch, method=method) for task_batch in task_pool]
        for future in as_completed(futures):
            print(future.result())
    print("DONE!")
    return True

def run_task_batch(task_batch,method):
    for idx,name,smi,op_name,seed in task_batch:
        operations = {"cv": ("connect","vinyl","ethynyl"),
                      "cve": ("connect","vinyl","ethynyl"),
                      "cv24p": ("connect","vinyl","annulate_pl2","annulate_pl4","phenyl"),
                      "cve24p": ("connect","vinyl","ethynyl","annulate_pl2","annulate_pl4","phenyl")}[op_name]
        if method == "tanimoto":
            rds = RediscoveryTanimoto(smi)
        elif method == "substructure":
            rds = RediscoverySubstructure(smi,highsymm=(name in ["fullerene"]))
        hc = HillClimbing(title=f"{method}_{idx}_{name}_{op_name}_{seed}", score_func=rds.score, root="[H][H]", operations=operations, random_seed=seed)
        hc.multiple(2000)
        print(f"Finished: {idx}_{name}_{op_name}_{seed}")
    return True

if __name__ == "__main__":
    RDLogger.DisableLog('rdApp.*')

    with open('../data/target_molecules.csv', encoding='utf-8-sig', newline='') as f:
        items = f.read().split("\r\n")
        # items = [tuple(item.split(",")) for item in items[:10]]
        items = [tuple(item.split(",")) for item in items[:27]]
    
    threads = 12
    task_pool = []
    for i,(idx,name,smi) in enumerate(items):
        if i in []:
            continue
        for op_name in (["cve","cve24p"] if "#" in smi else ["cv","cve","cv24p","cve24p"]):
            for seed in range(42,52):
                task_pool.append((idx,name,smi,op_name,seed))

    # run_task_pool(task_pool=task_pool, method="tanimoto", threads=threads)
    run_task_pool(task_pool=task_pool, method="substructure", threads=threads)



    
