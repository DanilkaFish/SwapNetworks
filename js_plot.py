from typing import List

import json
import numpy as np
from numpy import array

from Ternary_Tree.qiskit_interface.circuit_provider import get_file_name, Circuits
from my_utils.tikz_graph import * 
from scipy.stats import linregress
from copy import deepcopy as copy

def get_attrs(file_name, attrs: tuple):
    energy = []
    probs = []
    gate_count = ""
    datas = {i : [] for i in attrs}
    datas["gate_count"] = []
    with open(file_name, 'r') as rf:
        data = json.load(rf)
        for obj in data:
            for label in attrs:
                datas[label].append(obj[label])
    return datas


def energy_error_add_axes(num, lines, noise):
    for index, method in enumerate(methods[:]):
        # index = inde +  4
        datas = get_attrs(get_file_name(mol_name, noise, method), ("ref_ener", "energy", "prob", "gate_count"))
        if noise == "sc":
            lines[index].Y, lines[index].X = np.array(datas["energy"]), np.array(datas["prob"])
        else:
            if noise == "D":
                probs = np.array(datas["prob"])
                inf = 1 - np.sqrt(1 - 3/4*(1 - probs))
            else:
                probs = np.array(datas["prob"])
                inf = 1 - np.sqrt(1 - (1 - probs))
            lines[index].Y, lines[index].X = np.array(datas["energy"]), inf
        lines[index].Y = -abs(lines[index].Y - np.array(datas["ref_ener"]))/np.array(datas["ref_ener"])
        lines[index].legend = legends[index] + f" cx={datas["gate_count"][0]["cx"]}"
    axes.append(Axis(lines[:], "$\\lambda$", "$|E-E_0|$", AxisScale.usual(), xshift=num*8, title=noise))
    
def dn_add_axes(num, lines, noise, ne=2):
    for index, method in enumerate(methods[:]):
    # index = inde +  4
        datas = get_attrs(get_file_name(mol_name, noise, method), ("ref_ener", "energy", "prob", "gate_count", "addition_res:"))
        ar = np.array(datas["addition_res:"])
        n2 = ar[:,1] - 2*ne*ar[:,0] + ne*ne
        if noise == "sc":
            lines[index].Y, lines[index].X = n2, np.array([0.01, 0.1, 0.5, 1, 2, 3, 5])
        else:
            if noise == "D":
                probs = np.array(datas["prob"])
                inf = 1 - np.sqrt(1 - 3/4*(1 - probs))
            else:
                probs = np.array(datas["prob"])
                inf = 1 - np.sqrt(1 - (1 - probs))
            lines[index].Y, lines[index].X = n2, inf
        lines[index].legend = legends[index] + f" cx={datas["gate_count"][index]["cx"]}"
    axes.append(Axis(lines[:], "$1 -\\mathcal{{F}}$", f"$\\langle(n-{ne})^2\\rangle$", AxisScale.loglog(), xshift=num*8,yshift=-8, title=noise))

def round_with_uncertainty(value, uncertainty):
    # if uncertainty == 0 or np.isnan(uncertainty):
    #     return f"{value}", f"{uncertainty}"
    
    exponent = int(np.floor(np.log10(abs(uncertainty))))
    decimals = -exponent + 1 if exponent < 0 else 0
    if (uncertainty * (decimals << 1) < 30 ):
        # print(uncertainty * (1 << decimals) )
        # decimals -= 1
        pass
    decimals -= 2
    if decimals < 0:
        decimals = 0
    rounded_uncertainty = round(uncertainty, decimals)
    rounded_value = round(value, decimals)

    fmt = f"{{:.{decimals}f}}"
    return fmt.format(rounded_value) 
    # return "$" + fmt.format(rounded_value) + " \pm " + fmt.format(rounded_uncertainty) + "$"
s = '\n'.join([
    ">{\\centering\\arraybackslash}m{5em}", 
    ">{\\centering\\arraybackslash}m{4em}", 
    ">{\\centering\\arraybackslash}m{4em}", 
    ">{\\centering\\arraybackslash}m{4em}", 
    ">{\\centering\\arraybackslash}m{4em}", 
    ">{\\centering\\arraybackslash}m{6em}"
]
)
s = s 
def plot_table():
    table = Table("$E = \\beta p + E_0$", ["Method", *[noise for noise in noises], "\\# CX gates"], align=s, bordered=False)
    pre_table = []
    for index, method in enumerate(methods[:]):
        B = []
        for noise in noises:
        # index = inde +  4
            datas = get_attrs(get_file_name(mol_name, noise, method), ("ref_ener", "energy", "prob", "gate_count"))
            liney, linex = -abs(np.array(datas["energy"]) - np.array(datas["ref_ener"]))/np.array(datas["ref_ener"]), np.array(datas["prob"])
            liney = np.array(datas["energy"])
            if noise == "sc":
                pass
            else:
                if noise == "D":
                    linex = 1 - np.sqrt(1 - 3/4*(1 - linex))
                else:
                    linex = 1 - np.sqrt(1 - (1 - linex))
            liney = liney[5:]
            # linex = np.log(linex[3:])
            # liney = np.log(liney[3:])
            linex = linex[5:]

            # A,b = np.polyfit(linex, liney, 1)
            result = linregress(linex, liney)
            B.append(round_with_uncertainty(result.slope, result.stderr))
            # B.append(round_with_uncertainty(result.intercept, result.intercept_stderr))
        pre_table.append([legends[index], *B, datas["gate_count"][0]["cx"]])
    for j in range(1, len(pre_table[0])):
        coef = 2**18
        index = -1
        for i in range(len(pre_table)):
            if float(pre_table[i][j]) < coef:
                index = i
                coef = float(pre_table[i][j]) 
            pre_table[i][j] = f"${pre_table[i][j]}$"
        pre_table[index][j] = f"$\\mathbf{{{pre_table[index][j][1:-1]}}}$"
    for index, method in enumerate(methods[:]):
        table.add_row(*pre_table[index])

    print(table.generate_table())

if __name__ == "__main__":
    jw_line = Line(None, None, Marks.triangle(), "blue", "JW")
    bk_line = Line(None, None, Marks.triangle(), "red", "BK")
    jw_opt_line = Line(None, None, Marks.square(), "blue", "JWO")
    bk_opt_line = Line(None, None, Marks.square(), "red", "BKO")
    swap2xn_line = Line(None, None, Marks.pentagon(), "black", "MSN")
    swapgens_line = Line(None, None, Marks.pentagon(), "orange", "FSN S")
    swapgeny_line = Line(None, None, Marks.pentagon(), "green", "FSN Y")
    swap2xnshort_line = Line(None, None, Marks.pentagon(), "black, dashed", "MSN S")
    swap2xnyor_line = Line(None, None, Marks.pentagon(), "gray, dashed,line width=0.5pt", "MSN Y")
    legends = ["JW", "BK","JW opt", "BK opt", "FSN short", "FSN yor", "MSN short", "MSN yor"]
    mol_name = "datah2h2ExcSolProb/DH4" 
    mol_name = "datah2h2noise_level/DH4" 
    mol_name = "data/adapt_vqe_h2_8/"
    mol_name = "data/vqe_sc2/"
    mol_name = "data/LiH_8"
    # -1.8716797649675656
    # mol_name, ref_en = ("datah2_4/H2_4", -1.8573730129353947)
    # mol_name, ref_en = ("data02/H2_8", -1.8716797649675656)
    # noises = ["sc"]
    # noises = ["X","Z"]
    noises = ["D","X","Y","Z"]
    methods = Circuits.get_circs_names()[:4] + Circuits.get_circs_names()[5:] 
    # methods = methods[:2] + methods[4:]
    axes = []
    plot_table()
    # for num, noise in enumerate(noises):
        # lines = copy([jw_line, bk_line, jw_opt_line, bk_opt_line,  swapgens_line, swapgeny_line, swap2xnshort_line, swap2xnyor_line][:])
        # energy_error_add_axes()
    #     # dn_add_axes(4)

    # gd = GraphData(axes)
    # print(gd.generate_tikz())
