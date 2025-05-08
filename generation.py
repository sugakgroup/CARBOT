import glob
import pandas as pd
from rdkit import Chem, RDLogger
import shutil
import time

from carbot import *

def save_tmt_data(tmt,step,edge_type="unit",is_EZ=True):
    gen = tmt.num_generation
    pd.to_pickle(tmt,(f'data/tmt_{edge_type}_{"EZ" if is_EZ else "noEZ"}_dist{gen}_step{step}.pkl'))
    df = [{**{"smi": k},**v} for k,v in tmt.node_data.nodeinfo.items()]
    df = pd.json_normalize(df)
    df.to_csv(f'data/tmt_{edge_type}_{"EZ" if is_EZ else "noEZ"}_node_{gen}.csv', index=False, encoding='utf-8')
    df.to_pickle(f'data/tmt_{edge_type}_{"EZ" if is_EZ else "noEZ"}_node_{gen}.pkl')

    
from local.evolve import evolve

if __name__ == "__main__":
    RDLogger.DisableLog('rdApp.*')

    # Set Edge Type
    # edge_type = "unit"
    edge_type = "comp"
    # is_EZ = True
    is_EZ = False

    # Set Methods (Edge Types) for Transformation if you want
    evolmethods = tuple([])

    # Check Generation
    tmts = glob.glob(f'data/tmt_{edge_type}_{"EZ" if is_EZ else "noEZ"}_dist*_step1.pkl')
    gen = len(tmts)

    # Evolve for one distance
    if gen == 0:
        evolmethods = ("connect","vinyl","ethynyl") if edge_type == "unit" else (("connect","vinyl","ethynyl","annulate_pl2","annulate_pl4","phenyl") if edge_type == "comp" else evolmethods)
        tmt = Transmuter(original_structures=["[H][H]"],evolmethods=evolmethods,is_EZ=is_EZ)
    else:
        tmt = pd.read_pickle(f'data/tmt_{edge_type}_{"EZ" if is_EZ else "noEZ"}_dist{gen-1}_step1.pkl')
        tmt.evolve_single()
    save_tmt_data(tmt,1,edge_type=edge_type,is_EZ=is_EZ)

