from typing import List

import json
import numpy as np
from numpy import array

from .utils import get_file_name, Circuits
from my_utils.tikz_graph import * 


def get_en_prob(file_name):
    energy = []
    probs = []
    with open(file_name, 'r') as rf:
        data = json.load(rf)
        for obj in data:
            # dic[obj["noise gates"]]
            energy.append(obj["energy"])
            probs.append(obj["prob"])
    return array(energy), 1 - array(probs)

if __name__ == "__main__":
    jw_line = Line(None, None, Marks.triangle(), "blue", "ДВ")
    bk_line = Line(None, None, Marks.triangle(), "red", "БК")
    jw_opt_line = Line(None, None, Marks.triangle(), "blue", "ДВ опт")
    bk_opt_line = Line(None, None, Marks.triangle(), "red", "БК опт")
    swap2xn_line = Line(None, None, Marks.triangle(), "green", "SWAP 2xn")
    swapgens_line = Line(None, None, Marks.triangle(), "black", "SWAP short")
    swapgeny_line = Line(None, None, Marks.triangle(), "gray", "SWAP yor")

    mol_name, ref_en = "H2_8", -1.8716797649675605
    noises = ["D", "X", "Y", "Z"]
    methods = Circuits.get_circs_names()
    axes = []
    for noise in noises:
        lines=[jw_line, bk_line, jw_opt_line, bk_opt_line, swap2xn_line, swapgens_line, swapgeny_line]
        for index, method in enumerate(methods):
            lines[index].Y, lines[index].X = get_en_prob(get_file_name(mol_name, noise, method)) 
            lines[index].Y = -abs(lines[index].Y - ref_en)/ref_en
        axes.append(Axis(lines, "error parameter", "$\\frac{|E-E_0|}{|E_0|}$", AxisScale.logx()))
    gd = GraphData(axes)