
from evolve import *
from utils import num_unit_method
import pandas as pd
from collections import defaultdict

def edge_extract():
    df_unit = pd.read_pickle(f'./data/tmt_unit_noEZ_node_5.pkl')
    df_comp = pd.read_pickle(f'./data/tmt_comp_noEZ_node_3.pkl')

    nodes = defaultdict(int)
    edges = set()
    generation = dict()

    nodes["[H][H]"] = 3

    generation["[H][H]"] = 0

    cnt = 0
    for j in df_unit.itertuples():
        if cnt%1000 == 0:
            print(cnt,j.generation)
        cnt += 1
        if j.generation == 5:
            break
        for gen_smi, info in evol_single(j.smi,evolmethods=("connect","vinyl","ethynyl"),is_EZ=False):
            edges.add((j.smi,gen_smi,info[1]))
            nodes[gen_smi] |= 1
            generation[gen_smi] = sum(num_unit_method(gen_smi))

    cnt = 0
    for j in df_comp.itertuples():
        if cnt%1000 == 0:
            print(cnt,j.generation)
        cnt += 1
        if j.generation == 3:
            break
        for gen_smi, info in evol_single(j.smi,evolmethods=("connect","vinyl","ethynyl","annulate_pl2","annulate_pl4","phenyl"),is_EZ=False):
            edges.add((j.smi,gen_smi,info[1]))
            nodes[gen_smi] |= 2
            generation[gen_smi] = sum(num_unit_method(gen_smi))

    pd.to_pickle(nodes,(f'data/analysis/nodes.pkl'))
    pd.to_pickle(edges,(f'data/analysis/edges.pkl'))
    pd.to_pickle(generation,(f'data/analysis/generations.pkl'))
    print(len(nodes.keys()))
    print(len(edges))

from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger
import numpy as np

def calc_fp_all():
    df_unit = pd.read_pickle(f'./data/tmt_unit_noEZ_node_9.pkl')
    df_comp = pd.read_pickle(f'./data/tmt_comp_noEZ_node_6.pkl')

    fps = dict()
    nodes = defaultdict(int)

    nodes["[H][H]"] = 3
    
    cnt = 0
    for j in df_unit.itertuples():
        if cnt%1000 == 0:
            print(cnt,j.generation)
        cnt += 1
        # if j.generation == 9:
        #     break
        if not j.smi in fps.keys():
            fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(j.smi), 2, 2048)
            arr = np.zeros((2048,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp, arr)
            fps[j.smi] = arr
        nodes[j.smi] |= 1
        cnt += 1

    cnt = 0
    for j in df_comp.itertuples():
        if cnt%1000 == 0:
            print(cnt,j.generation)
        cnt += 1
        # if j.generation == 6:
        #     break
        if not j.smi in fps.keys():
            fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(j.smi), 2, 2048)
            arr = np.zeros((2048,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp, arr)
            fps[j.smi] = arr
        nodes[j.smi] |= 2
        cnt += 1

    pd.to_pickle(nodes,(f'data/analysis/nodes.pkl'))
    pd.to_pickle(fps,(f'data/analysis/fps.pkl'))

from rdkit.Chem import Descriptors

def calc_properties():
    df_unit = pd.read_pickle(f'./data/tmt_unit_noEZ_node_9.pkl')
    df_comp = pd.read_pickle(f'./data/tmt_comp_noEZ_node_6.pkl')

    props = defaultdict(dict)
    nodes = defaultdict(int)

    nodes["[H][H]"] = 3


    desk_list = Descriptors.descList
    calcs = dict()
    for (i,j) in desk_list:
        if i in ("MolLogP","MolWt","NumRotatableBonds"):
            calcs[i] = j
    
    print(calcs)

    cnt = 0
    for j in df_unit.itertuples():
        if cnt%1000 == 0:
            print(cnt,j.generation)
        # if j.generation == 9:
        #     break
        if not j.smi in props.keys():
            mol = Chem.MolFromSmiles(j.smi)
            for desc,calc in calcs.items():
                props[j.smi][desc] = calc(mol)
        nodes[j.smi] |= 1
        cnt += 1

    cnt = 0
    for j in df_comp.itertuples():
        if cnt%1000 == 0:
            print(cnt,j.generation)
        # if j.generation == 6:
        #     break
        if not j.smi in props.keys():
            mol = Chem.MolFromSmiles(j.smi)
            for desc,calc in calcs.items():
                props[j.smi][desc] = calc(mol)
        nodes[j.smi] |= 2
        cnt += 1
    
    pd.to_pickle(nodes,(f'data/analysis/nodes.pkl'))
    pd.to_pickle(props,(f'data/analysis/props.pkl'))

from evolve import evol_single
import queue
from rediscovery import *

def retrosynthesis(smi,name):
    rds = Rediscovery(smi)

    edges = set()
    nodes = dict()
    nodes_pool = queue.Queue()
    nodes_pool.put("[H][H]")
    nodes["[H][H]"] = 1/(sum(num_unit_method(smi))+1)

    while not nodes_pool.empty():   
        smi_p = nodes_pool.get()
        for smi_c, info in evol_single(smi_p, evolmethods=("connect","vinyl","ethynyl") , is_EZ=False):
            score = rds.score(smi_c)
            if score > 0:
                edges.add((smi_p,smi_c,info[1]))
                if not smi_c in nodes.keys():
                    print(smi_c,score)
                    nodes[smi_c] = score
                    nodes_pool.put(smi_c)

    pd.to_pickle(edges,(f'data/analysis/retrosynthesis_edges_{name}_unit.pkl'))
    pd.to_pickle(nodes,(f'data/analysis/retrosynthesis_nodes_{name}.pkl'))
    print(len(edges))
    print(len(nodes))

    # edges = pd.read_pickle(f'../data/analysis/retrosynthesis_edges_{name}_unit.pkl')
    # nodes = pd.read_pickle(f'../data/analysis/retrosynthesis_nodes_{name}.pkl')

    for smi_p in nodes.keys():   
        for smi_c, info in evol_single(smi_p, evolmethods=("annulate_pl2","annulate_pl4","phenyl") , is_EZ=False):
            score = rds.score(smi_c)
            if rds.score(smi_c) > 0:
                print(smi_c,score)
                edges.add((smi_p,smi_c,info[1]))

    pd.to_pickle(edges,(f'data/analysis/retrosynthesis_{name}_comp.pkl'))
    print(len(edges))
    print(len(nodes))

RDLogger.DisableLog('rdApp.*')
retrosynthesis(smi="C12=CC=CC(C3=C(C4=CC=C5)C5=CC=C3)=C1C4=CC=C2",name="perylene")




