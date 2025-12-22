from logging import root
from rdkit import RDLogger

from concurrent.futures import ProcessPoolExecutor, as_completed

import csv
import random
import sys
import time
from pathlib import Path
import heapq

proj_root = (Path.cwd().parent).resolve()
sys.path.insert(0, str(proj_root))

from src.evolve import *
from src.degen import *
from src.rediscovery import *
from src.utils import *

class HillClimbingEvoMol:
    def __init__(self,
                 title,
                 score_func,
                 root="[H][H]",
                 operations=("connect","vinyl","ethynyl"),
                 is_degen=True,
                 random_seed=42,
                 max_num_search_improver=50,
                 base_dir=None,
                 save_dir="2_rediscovery"):
        if base_dir is None:
            base_dir = Path.cwd().parent
        base_dir = Path(base_dir)
        self.out_dir = base_dir / "data" / save_dir / "evomol"
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.title = title

        self.time0 = time.time()
        self.root = root
        self.score_func = score_func
        random.seed(random_seed)
        self.children = dict()
        if is_degen:
            self.parents = dict()
        self.used_smi = set()
        self.operations = operations
        self.is_degen = is_degen
        self.population_size = 1000
        self.pop = [(-99999,"") for _ in range(self.population_size)]
        score = score_func(root)
        self.pop[0] = ((score,root))
        self.pop_set = set([root])
        self.round = 0
        self.num_replace = 10
        self.total_eval = 1
        self.max_num_search_improver = max_num_search_improver

        with open(self.out_dir / f"{self.title}.out", "w", encoding="utf-8") as f_log:
            f_log.write(f"[{time.time()-self.time0}]Start running EvoMol from {root}...\n")

        with open(self.out_dir / f"{self.title}.csv", "w", newline="", encoding="utf-8") as f:
            writer_nodes = csv.writer(f)
            writer_nodes.writerow(["TotalEval","Round","NumEvalRound","SMILES","Score","TimeElapsed"])
            writer_nodes.writerow([self.total_eval,self.round,1,root,score,time.time()-self.time0])
    
        with open(self.out_dir / f"{self.title}_pop.csv", "w", newline="", encoding="utf-8") as f:
            writer_nodes = csv.writer(f)
            writer_nodes.writerow(["Round","SMILES","Score","Action"])
            writer_nodes.writerow([self.round,root,score,"Add"])
    
    def single(self):
        is_improved = False
        self.round += 1
        eval_round = 0

        to_be_replaced = [item+(1,) for item in self.pop[-self.num_replace:]]
        heapq.heapify(to_be_replaced)

        num_update = 0
        for score,smi in self.pop:
            if score < -1.0:
                break
            else:
                candidate = set()
                if not smi in self.children.keys():
                    self.children[smi] = []
                    for _, gen_smis in evol_single(smi,self.operations).items():
                        self.children[smi].extend(gen_smis)
                    self.children[smi] = tuple(sorted(list(self.children[smi])))
                for gen_smi in self.children[smi]:
                    if not gen_smi in self.used_smi:
                        candidate.add(gen_smi)
                if self.is_degen:
                    if not smi in self.parents.keys():
                        self.parents[smi] = []
                        for _, degen_smis in degen_single(smi,self.operations).items():
                            self.parents[smi].extend(degen_smis)
                        self.parents[smi] = tuple(sorted(list(self.parents[smi])))
                    for degen_smi in self.parents[smi]:
                        if not degen_smi in self.used_smi:
                            candidate.add(degen_smi)

                candidate = sorted(list(candidate))
                random.shuffle(candidate)
                
                num_tries = 0
                for gen_smi in candidate:
                    num_tries += 1
                    if not gen_smi in self.used_smi:
                        score_gensmi = self.score_func(gen_smi)
                        self.used_smi.add(gen_smi)
                        self.total_eval += 1
                        eval_round += 1
                        with open(self.out_dir / f"{self.title}.csv", "a", newline="", encoding="utf-8") as f:
                            writer_nodes = csv.writer(f)
                            writer_nodes.writerow([self.total_eval,self.round,eval_round,gen_smi,score_gensmi,time.time()-self.time0])
                        if score_gensmi > to_be_replaced[0][0]:
                            if to_be_replaced[0][2] == 1:
                                num_update += 1
                                is_improved = True
                                num_tries = 0
                            heapq.heappop(to_be_replaced)
                            heapq.heappush(to_be_replaced, (score_gensmi,gen_smi,0))
                        if score_gensmi == 1.0:
                            with open(self.out_dir / f"{self.title}_pop.csv", "a", newline="", encoding="utf-8") as f:
                                writer_nodes = csv.writer(f)
                                for i in range(self.num_replace):
                                    if self.pop[-i-1][0] != -99999:
                                        writer_nodes.writerow([self.round,self.pop[-i-1][1],self.pop[-i-1][0],"Remove"])
                                    item = heapq.heappop(to_be_replaced)
                                    if item[0] != -99999:
                                        writer_nodes.writerow([self.round,item[1],item[0],"Add"])
                                    self.pop[-i-1] = (item[0], item[1])
                            with open(self.out_dir / f"{self.title}.out", "a", encoding="utf-8") as f_log:
                                f_log.write(f"[{time.time()-self.time0}] Reproduced in round {self.round}!!\n")
                            return False
                    if num_tries >= self.max_num_search_improver or num_update >= self.num_replace:
                        break
                
                if num_update >= self.num_replace:
                    break
        
        if not is_improved:
            with open(self.out_dir / f"{self.title}.out", "a", encoding="utf-8") as f_log:
                f_log.write(f"[{time.time()-self.time0}] No improvement found in round {self.round}. Terminating...\n")
            return False
        
        with open(self.out_dir / f"{self.title}_pop.csv", "a", newline="", encoding="utf-8") as f:
            writer_nodes = csv.writer(f)
            for i in range(self.num_replace):
                if self.pop[-i-1][0] != -99999:
                    writer_nodes.writerow([self.round,self.pop[-i-1][1],self.pop[-i-1][0],"Remove"])
                item = heapq.heappop(to_be_replaced)
                if item[0] != -99999:
                    writer_nodes.writerow([self.round,item[1],item[0],"Add"])
                self.pop[-i-1] = (item[0], item[1])

        self.pop.sort(reverse=True,key=lambda x: x[0])
        with open(self.out_dir / f"{self.title}.out", "a", encoding="utf-8") as f_log:
            f_log.write(f"[{time.time()-self.time0}] Updated in round {self.round}! max:{self.pop[0]}, min:{self.pop[-1]}\n")
        
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
        elif method in {"substructure","substructureforward"}:
            rds = RediscoverySubstructure(smi,highsymm=(name in ["fullerene"]))
        if method == "substructureforward":
            evml = HillClimbingEvoMol(title=f"{method}_{idx}_{name}_{op_name}_{seed}", score_func=rds.score, root="[H][H]", operations=operations, random_seed=seed, is_degen=False)
        else:
            evml = HillClimbingEvoMol(title=f"{method}_{idx}_{name}_{op_name}_{seed}", score_func=rds.score, root="[H][H]", operations=operations, random_seed=seed)
        evml.multiple(2000)
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
    # run_task_pool(task_pool=task_pool, method="substructure", threads=threads)
    run_task_pool(task_pool=task_pool, method="substructureforward", threads=threads)
