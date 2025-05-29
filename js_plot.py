from typing import List

import json
import numpy as np
from numpy import array

from utils import get_file_name, Circuits
from my_utils.tikz_graph import * 


def get_attrs(file_name, attrs: tuple):
    energy = []
    probs = []
    gate_count = ""
    datas = {i : [] for i in attrs}
    datas["gate_count"] = []
    print(file_name)
    with open(file_name, 'r') as rf:
        data = json.load(rf)
        for obj in data:
            datas["gate_count"].append(obj["gate_count"]["cx"])
            for label in attrs:
            # dic[obj["noise gates"]]
                datas[label].append(obj[label])
            # probs.append(obj["prob"])
    return datas


def curve_plot():
    jw_line = Line(None, None, Marks.triangle(), "blue", "JW")
    bk_line = Line(None, None, Marks.triangle(), "red", "BK")
    jw_opt_line = Line(None, None, Marks.square(), "blue", "JW opt")
    bk_opt_line = Line(None, None, Marks.square(), "red", "BK opt")
    swap2xn_line = Line(None, None, Marks.pentagon(), "green", "MSN")
    swapgens_line = Line(None, None, Marks.pentagon(), "black", "FSN short")
    swapgeny_line = Line(None, None, Marks.pentagon(), "gray", "FSN yor")
    ideal_line = Line(None, None, Marks.mro(), "black", "presice")
    ccsd_line = Line(None, None, Marks.square(), "black", "CCSD")
    legends = ["JW", "BK","JW opt", "BK opt", "SWAP 2xn", "SWAP short", "SWAP yor"]
    noises = ["D", "X", "Y", "Z"]
    mol_name = "datah2h2/DH4" 
    methods = Circuits.get_circs_names()
    axes = []
    for noise in noises:
        lines=[jw_line, bk_line, jw_opt_line, bk_opt_line, swap2xn_line, swapgens_line, swapgeny_line,  ccsd_line]
        for index, method in enumerate(methods[:]):
            datas = get_attrs(get_file_name(mol_name, noise, method), ("ref_ener", "energy", "dist"))
            lines[index].Y, lines[index].X = np.array(datas["energy"]), np.linspace(1,1.6,12)
            # lines[index].Y, lines[index].X = np.array(datas["energy"]), np.array(datas["dist"])
            lines[index].Y = lines[index].Y - np.array(datas["ref_ener"])
            # lines[index].legend = legends[index] + f" ({gc["cx"]})"
            # lines[index].legend = legends[index]
        datas = get_attrs(get_file_name(mol_name, noise, methods[-1]), ("ref_ener", "dist"))
        # lines[-2].Y, lines[-2].X = np.array(datas["ref_ener"]), np.array(datas["dist"])
        lines[-1].Y, lines[-1].X = np.array([-4.618176584268548, -4.5267440656257065, -4.443973799442764, 
                                             -4.372587187298123, -4.316057866013907, -4.267417963568857, 
                                             -4.228030673716043, -4.195361580192051, -4.166272967536891, 
                                             -4.139158229521049, -4.113257779037612, -4.088214875105436
                                             ]), np.linspace(1,1.6,12)
        lines[-1].Y = lines[-1].Y - np.array(datas["ref_ener"])
        
        axes.append(Axis(lines, "$\\lambda$", "$E - E_0$", AxisScale.usual()))
    gd = GraphData(axes)
    print(gd.generate_tikz())
    
    
if __name__ == "__main__":
    jw_line = Line(None, None, Marks.triangle(), "blue", "JW")
    bk_line = Line(None, None, Marks.triangle(), "red", "BK")
    jw_opt_line = Line(None, None, Marks.square(), "blue", "JW opt")
    bk_opt_line = Line(None, None, Marks.square(), "red", "BK opt")
    swap2xn_line = Line(None, None, Marks.pentagon(), "green", "MSN")
    swapgens_line = Line(None, None, Marks.pentagon(), "black", "FSN short")
    swapgeny_line = Line(None, None, Marks.pentagon(), "gray", "FSN yor")
    legends = ["JW", "BK","JW opt", "BK opt", "MSN", "FSN short", "FSN yor"]
    # legends = legends[:2] + legends[4:]
    mol_name = "datah2h2ExcSolProb/DH4" 
    mol_name = "datah2h2noise_level/DH4" 
    # -1.8716797649675656
    # mol_name, ref_en = ("datah2_4/H2_4", -1.8573730129353947)
    # mol_name, ref_en = ("data02/H2_8", -1.8716797649675656)
    noises = ["D", "Y", "sc"]
    noises = ["X","Z"]
    # noises = ["D", "Z"][0:1]
    methods = Circuits.get_circs_names()[-3:]
    # methods = methods[:2] + methods[4:]
    axes = []
    for noise in noises:
        lines = [jw_line, bk_line, jw_opt_line, bk_opt_line, swap2xn_line, swapgens_line, swapgeny_line][-3:]
        # lines=lines[:2] + lines[4:]
        for index, method in enumerate(methods[:]):
            # index = inde +  4
            datas = get_attrs(get_file_name(mol_name, noise, method), ("ref_ener", "energy", "prob"))
            lines[index].Y, lines[index].X = np.array(datas["energy"]), 1 - np.array(datas["prob"])
            lines[index].Y = -abs(lines[index].Y - np.array(datas["ref_ener"]))/np.array(datas["ref_ener"])
            if noise in {"D", "sc"}:
                lines[index].X = 1 - np.sqrt(1 - 3 * lines[index].X/4)
            else:
                lines[index].X = 1 - np.sqrt(1 - 4 * lines[index].X/4)
                
            lines[index].legend = legends[index + 4] + f" ({datas["gate_count"][0]})"
            # lines[index].legend = legends[index]
        axes.append(Axis(lines, "$\\lambda$", "$|E-E_0|$", "log" + AxisScale.loglog()))
    gd = GraphData(axes)
    print(gd.generate_tikz())
    # curve_plot()