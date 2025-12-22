from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from pathlib import Path
import sys
import time

proj_root = (Path.cwd().parent).resolve()
sys.path.insert(0, str(proj_root))

from rdkit import RDLogger

from src.evolve import evol_single
from src.utils import *
from src.rediscovery import *

class Enumeration:
    def __init__(self,
                 title,
                 initial_structures=("[H][H]",),
                 operations=("connect","vinyl","ethynyl"),
                 is_EZ=False,
                 edge_save=False,
                 threads=16,
                 base_dir=None,
                 save_dir="1_enumeration"):
        if base_dir is None:
            base_dir = Path.cwd().parent
        base_dir = Path(base_dir)
        self.out_dir = base_dir / "data" / save_dir
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.title = title
        
        self.nodeinfo = defaultdict(dict)
        self.generated_structures = set()
        self.structures_each_generation = []
        self.num_generation = 0
        self.operations = operations
        self.is_EZ = is_EZ
        self.threads = threads
        self.edge_save = edge_save
                    
        with open(self.out_dir / f"{self.title}.out", "w", encoding="utf-8") as f_log:
            f_log.write("---------- Generation #0 ----------\n")
            f_log.write(f"Structures {initial_structures} are loaded as the original structures.\n")
        
        with open(self.out_dir / f"{self.title}_nodes.csv", "w", newline="", encoding="utf-8") as f_nodes:
            writer_nodes = csv.writer(f_nodes)
            writer_nodes.writerow(["SMILES","Generation"])
        
        if self.edge_save:
            for operation in operations:
                with open(self.out_dir / f"{self.title}_edges_{operation}.csv", "w", newline="", encoding="utf-8") as f_edges:
                    writer_edges = csv.writer(f_edges)
                    writer_edges.writerow(["Parent_SMILES","Child_SMILES","Generation"])

        self.structures_each_generation.append([])
        with open(self.out_dir / f"{self.title}_nodes.csv", "a", newline="", encoding="utf-8") as f_nodes:
            writer_nodes = csv.writer(f_nodes)
            for smi in initial_structures:
                self.generated_structures.add(smi)
                self.structures_each_generation[-1].append(smi)
                writer_nodes.writerow([smi, self.num_generation])
        
        with open(self.out_dir / f"{self.title}.out", "a", encoding="utf-8") as f_log:
            f_log.write(f"# {len(self.structures_each_generation[-1])} structures in generation #0"+"\n")
    
    def evolve_single(self):
        time_0 = time.time()
        self.num_generation += 1
        with open(self.out_dir / f"{self.title}.out", "a", encoding="utf-8") as f_log:
            f_log.write("---------- Generation #"+str(self.num_generation)+" ----------"+"\n")
        self.structures_each_generation.append([])
        threads = self.threads
        parent_smis_pool = [[] for _ in range(threads*5)]
        for i,smi in enumerate(self.structures_each_generation[-2]):
            parent_smis_pool[i%(threads*5)].append(smi)
        with ProcessPoolExecutor(max_workers=threads) as ex:
            futures = [ex.submit(evolve_batch,
                                 smis=smis,
                                 operations=self.operations,
                                 is_EZ=self.is_EZ) for smis in parent_smis_pool]
            for future in as_completed(futures):
                with open(self.out_dir / f"{self.title}_nodes.csv", "a", newline="", encoding="utf-8") as f_nodes:
                    writer_nodes = csv.writer(f_nodes)
                    for operation, smis in future.result().items():
                        for parent_smi, child_smi in smis:
                            if not (child_smi in self.generated_structures):
                                self.generated_structures.add(child_smi)
                                self.structures_each_generation[-1].append(child_smi)
                                writer_nodes.writerow([child_smi, self.num_generation])
                if self.edge_save:
                    for operation, smis in future.result().items():
                        with open(self.out_dir / f"{self.title}_edges_{operation}.csv", "a", newline="", encoding="utf-8") as f_edges:
                            writer_edges = csv.writer(f_edges)
                            for parent_smi, child_smi in smis:
                                writer_edges.writerow([parent_smi,child_smi,self.num_generation])

        with open(self.out_dir / f"{self.title}.out", "a", encoding="utf-8") as f_log:
            f_log.write(f"Time for generation of {self.num_generation}: {time.time()-time_0}"+"\n")
            f_log.write(f"# {len(self.structures_each_generation[-1])} structures in generation #{self.num_generation}"+"\n")
              
def evolve_batch(smis,operations,is_EZ):
    gen_dict = defaultdict(set)
    for smi in smis:
        for operation, gen_smis in evol_single(smi=smi, operations=operations, is_EZ=is_EZ).items():
            for gen_smi in gen_smis:
                gen_dict[operation].add((smi,gen_smi))
    return gen_dict
   
if __name__ == '__main__':
    RDLogger.DisableLog('rdApp.*')
    # original_strucutres as list of SMILES
    enum = Enumeration(title="cv", initial_structures=("[H][H]",), operations=("connect","vinyl"), is_EZ=False, edge_save=True, threads=16, save_dir="1_enumeration")
    for _ in range(9999):
        enum.evolve_single()
        if len(enum.structures_each_generation) <= 1 or len(enum.structures_each_generation[-1])/len(enum.structures_each_generation[-2])*len(enum.structures_each_generation[-1]) > 10000000:
            break

    enum = Enumeration(title="cve", initial_structures=("[H][H]",), operations=("connect","vinyl","ethynyl"), is_EZ=False, edge_save=True, threads=16, save_dir="1_enumeration")
    for _ in range(9999):
        enum.evolve_single()
        if len(enum.structures_each_generation) <= 1 or len(enum.structures_each_generation[-1])/len(enum.structures_each_generation[-2])*len(enum.structures_each_generation[-1]) > 10000000:
            break

    enum = Enumeration(title="cv24p", initial_structures=("[H][H]",), operations=("connect","vinyl","annulate_pl2","annulate_pl4","phenyl"), is_EZ=False, edge_save=True, threads=16, save_dir="1_enumeration")
    for _ in range(9999):
        enum.evolve_single()
        if len(enum.structures_each_generation) <= 1 or len(enum.structures_each_generation[-1])/len(enum.structures_each_generation[-2])*len(enum.structures_each_generation[-1]) > 10000000:
            break
    
    enum = Enumeration(title="cve24p", initial_structures=("[H][H]",), operations=("connect","vinyl","ethynyl","annulate_pl2","annulate_pl4","phenyl"), is_EZ=False, edge_save=True, threads=16, save_dir="1_enumeration")
    for _ in range(9999):
        enum.evolve_single()
        if len(enum.structures_each_generation) <= 1 or len(enum.structures_each_generation[-1])/len(enum.structures_each_generation[-2])*len(enum.structures_each_generation[-1]) > 10000000:
            break

    enum = Enumeration(title="cv_EZ", initial_structures=("[H][H]",), operations=("connect","vinyl"), is_EZ=True, edge_save=False, threads=16, save_dir="1_enumeration")
    for _ in range(9999):
        enum.evolve_single()
        if len(enum.structures_each_generation) <= 1 or len(enum.structures_each_generation[-1])/len(enum.structures_each_generation[-2])*len(enum.structures_each_generation[-1]) > 10000000:
            break

    enum = Enumeration(title="cve_EZ", initial_structures=("[H][H]",), operations=("connect","vinyl","ethynyl"), is_EZ=True, edge_save=False, threads=16, save_dir="1_enumeration")
    for _ in range(9999):
        enum.evolve_single()
        if len(enum.structures_each_generation) <= 1 or len(enum.structures_each_generation[-1])/len(enum.structures_each_generation[-2])*len(enum.structures_each_generation[-1]) > 10000000:
            break

    enum = Enumeration(title="cv24p_EZ", initial_structures=("[H][H]",), operations=("connect","vinyl","annulate_pl2","annulate_pl4","phenyl"), is_EZ=True, edge_save=False, threads=16, save_dir="1_enumeration")
    for _ in range(9999):
        enum.evolve_single()
        if len(enum.structures_each_generation) <= 1 or len(enum.structures_each_generation[-1])/len(enum.structures_each_generation[-2])*len(enum.structures_each_generation[-1]) > 10000000:
            break
    
    enum = Enumeration(title="cve24p_EZ", initial_structures=("[H][H]",), operations=("connect","vinyl","ethynyl","annulate_pl2","annulate_pl4","phenyl"), is_EZ=True, edge_save=False, threads=16, save_dir="1_enumeration")
    for _ in range(9999):
        enum.evolve_single()
        if len(enum.structures_each_generation) <= 1 or len(enum.structures_each_generation[-1])/len(enum.structures_each_generation[-2])*len(enum.structures_each_generation[-1]) > 10000000:
            break





    


