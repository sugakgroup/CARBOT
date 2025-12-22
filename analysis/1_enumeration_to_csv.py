import csv
from pathlib import Path
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

proj_root = (Path.cwd().parent).resolve()
sys.path.insert(0, str(proj_root))

from rdkit import Chem
from rdkit.Chem import Descriptors

from src.rediscovery import *
from src.utils import *

class SkeletonMatching:
    def __init__(self,target_smi):
        self.target_smi_skeleton = target_smi.replace("c","C").replace("-","").replace("=","").replace("#","")
        self.target_mol_skeleton = Chem.MolFromSmiles(self.target_smi_skeleton)

    def score(self, smi):
        smi_skeleton = smi.replace("c","C").replace("-","").replace("=","").replace("#","")
        mol_skeleton = Chem.MolFromSmiles(smi_skeleton)
        return self.target_mol_skeleton.HasSubstructMatch(mol_skeleton)

def calculate_molecules(smi,calculate_items=()):
    data = dict()
    # MW
    if any([item in calculate_items for item in ["MW","clogP","RBC"]]):
        mol = Chem.MolFromSmiles(smi)
    if "MW" in calculate_items:
        data["MW"] = round(Descriptors.MolWt(mol),3)
    if "clogP" in calculate_items:
        data["clogP"] = round(Descriptors.MolLogP(mol),3)
    if "RBC" in calculate_items:
        data["RBC"] = round(Descriptors.NumRotatableBonds(mol),3)
    # connect, vinyl, ethynyl
    if any([item in calculate_items for item in ["num_connect","num_vinyl","num_ethynyl"]]):
        data["num_connect"], data["num_vinyl"], data["num_ethynyl"] = num_unit_method(smi)
    # benzene ring
    if "has_benzene" in calculate_items:
        rds = RediscoverySubstructure(smi)
        data["has_benzene"] = (rds.score("c1ccccc1") > 0.0)
    if any([item in calculate_items for item in ["has_ring_3","has_ring_4","has_ring_5","has_ring_6","has_ring_7","has_ring_8","is_bredt","is_small_fused"]]):
        rds_skeleton = SkeletonMatching(smi)
        # ring 3-8
        for ring_size in [i for i in range(3,9) if f"has_ring_{i}" in calculate_items]:
            ring_smi = "C1" + "".join(["C" for _ in range(ring_size-1)]) + "1"
            data[f"has_ring_{ring_size}"] = rds_skeleton.score(ring_smi)
        # Bredt
        if "is_bredt" in calculate_items:
            data["is_bredt"] = False
            for skeleton in ["C1(C2)CC2C1",
                            "C1(C2)CC2CC1",
                            "C12CCC(C2)CC1","C12CCCC(C2)C1",
                            "C12CCC(CC2)CC1","C12CCCC(CC2)C1","C1(C2)CCCCC2C1",
                            "C12CCCC(CC2)CC1","C12CCCC(CCC2)C1","C12CCCCC(CC2)C1","C1(CCCCC2)CC2C1",
                            "C12CCCC(CCC2)CC1","C12CCC(CCCC2)CC1","C12CCCCC(CCC2)C1","C12CCC(CCCCC2)C1","C1(CCCCCC2)CC2C1"]:
                if rds_skeleton.score(skeleton):
                    data["is_bredt"] = True
                    break
        # small fused
        if "is_small_fused" in calculate_items:
            data["is_small_fused"] = False
            for skeleton in ["C12CC1C2",
                            "C12CC1CC2",
                            "C12CC1CCC2","C12CCC1CC2",
                            "C12CC1CCCC2","C12CCC1CCC2"]:
                if rds_skeleton.score(skeleton):
                    data["is_small_fused"] = True
                    break
    return data

def worker(operations):
    with open(f"../data/1_enumeration/{operations}_nodes.csv", "r", encoding="utf-8") as f_in, open(f"./{operations}_analysis.csv", "w", newline="", encoding="utf-8") as f_out:
        reader = csv.reader(f_in)
        writer = csv.writer(f_out)
        header = ["smi","generation","MW","clogP","RBC","num_connect","num_vinyl","num_ethynyl","has_benzene","has_ring_3","has_ring_4","has_ring_5","has_ring_6","has_ring_7","has_ring_8","is_bredt","is_small_fused"]
        writer.writerow(header)
        for row in reader:
            if row[1] == "Generation":
                continue
            smi = row[0]
            if Chem.MolFromSmiles(smi) is None:
                continue
            generation = int(row[1])
            calc_data = calculate_molecules(smi,calculate_items=header[2:])
            writer.writerow([smi,generation,
                                calc_data.get("MW",""),
                                calc_data.get("clogP",""),
                                calc_data.get("RBC",""),
                                calc_data.get("num_connect",""),
                                calc_data.get("num_vinyl",""),
                                calc_data.get("num_ethynyl",""),
                                calc_data.get("has_benzene",""),
                                calc_data.get("has_ring_3",""),
                                calc_data.get("has_ring_4",""),
                                calc_data.get("has_ring_5",""),
                                calc_data.get("has_ring_6",""),
                                calc_data.get("has_ring_7",""),
                                calc_data.get("has_ring_8",""),
                                calc_data.get("is_bredt",""),
                                calc_data.get("is_small_fused","")])
    return operations


def main():
    with ProcessPoolExecutor(max_workers=4) as ex:
        futures = [ex.submit(worker,
                                operations=operations) for operations in ["cv","cve","cv24p","cve24p"]]
        for future in as_completed(futures):
            print(f"{future.result()} done!")

if __name__ == "__main__":
    main()