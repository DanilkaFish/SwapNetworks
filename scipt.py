from utils import Circuits, get_file_name
import json
import numpy as np


data = []
mol_name = "datah2h2copy/DH4" 

methods = Circuits.get_circs_names()
R = np.linspace(1,1.6,12)

for noise in ["D","X","Y","Z"]:
    for index, method in enumerate(methods[:]):
        file_name = get_file_name(mol_name, noise, method)
        with open(file_name, 'r') as rf:
            data = json.load(rf)
        for index in range(len(data)):
            data[index]["dist"] = R[index]
    with open(file_name, 'w') as file:
        json.dump(data, file, indent=4)   
                # probs.append(obj["prob"])