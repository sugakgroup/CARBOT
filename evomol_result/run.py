from evomol import run_model
from concurrent.futures import ProcessPoolExecutor, as_completed
import random
import numpy as np
import os
import contextlib

def main():
    with open(os.devnull, "w") as fnull, contextlib.redirect_stdout(fnull), contextlib.redirect_stderr(fnull):
        with open('target_molecules.csv', encoding='utf-8-sig', newline='') as f:
            items = f.read().split("\r\n")
        items = [tuple(item.split(",")) for item in items[10:13]]

        threads = 12
        task_pool = [[] for _ in range(threads*2)]

        now_i = 0

        for idx,name,smi in items:
            for seed in range(42,52):
                task_pool[now_i%(threads*2)].append((idx,name,smi,seed))
                now_i += 1
        task_pool = [tuple(task_batch) for task_batch in task_pool if len(task_batch) > 0]

        with ProcessPoolExecutor(max_workers=threads) as ex:
            futures = [ex.submit(run_task_batch, task_batch) for task_batch in task_pool]
            for future in as_completed(futures):
                pass

def run_task_batch(task_batch):
    for idx,name,smi,seed in task_batch:
        random.seed(seed)
        np.random.seed(seed)
        run_model({
            "obj_function": f"rediscovery_{smi}",
            "optimization_parameters": {
                "max_steps": 2000,
                "stop_kth_score_value": (1, 1.0, 5)
            },
            "action_space_parameters": {
                "atoms": "C"
            },
            "io_parameters": {
                "save_n_steps": 1,
                "print_n_steps": 100000,
                "model_path": f"examples/{idx}_{name}/rediscovery_seed{seed}"
            }
        })

if __name__ == "__main__":
    main()
