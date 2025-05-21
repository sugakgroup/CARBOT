# CARBOT
CARBOT (Conjugated Aromatic Rule-Based Organic Transformation) is a mathematically verified rule-based molecular generator for π-conjugated  hydrocarbons.


## Requirement
- python: 3.11.9
- rdkit: 2023.9.1
- pandas 2.1.4

*Other versions may work, but have not been tested.


## Installation
```sh
git clone https://github.com/sugakgroup/CARBOT.git
cd CARBOT
```

## Usages
The simplest use of CARBOT is below:
```python
from carbot import Transmuter

original_structures = ["[H][H]"] # hydrogen molecule
evolmethods = ("connect","ethylene","acetylene") # only w/ unit edges, but w/o edge sequence
# if you want to use edge sequence, write this instead:
# evolmethods = ("connect","ethylene","acetylene","annulate_pl2","annulate_pl4","phenyl")
tmt = Transmuter(original_structures,evolmethods)
tmt.evolve_single()

```



## License
This package is distributed under the MIT License.


## How to cite CARBOT
K. Suga*, H. Takahashi, K. Terayama, M. Sumita, S. Saito, _ChemRxiv_ **2025**.
*Manuscript in preparation. Before the DOI is issued, please cite this GitHub repository: https://github.com/sugakgroup/CARBOT
