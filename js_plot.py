from typing import List

import json
import numpy as np
from numpy import array

from Ternary_Tree.qiskit_interface.circuit_provider import get_file_name, Circuits
from my_utils import * 
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
    
    # exponent = int(np.floor(np.log10(abs(uncertainty))))
    # decimals = -exponent + 1 if exponent < 0 else 0
    # if (uncertainty * (decimals << 1) < 30 ):
    #     # print(uncertainty * (1 << decimals) )
    #     # decimals -= 1
    #     pass
    decimals = 0
    if decimals < 0:
        decimals = 0
    rounded_uncertainty = round(uncertainty, decimals)
    rounded_value = round(value, decimals)

    fmt = f"{{:.{decimals}f}}"
    return fmt.format(rounded_value) 


def plot_table(mol_name, gates, ne=None):
    caption = "Noise susceptibility $\\chi$ for the H2, 8 qubits."
    iterate = noises
    col_names = ["Method", *iterate, "\\# CX gates"]
    align_displ = [5, *[4 for i in iterate], 5]
    table = Table(caption, 
                  col_names, 
                  align=align_displ,
                  bordered=False,
                  placement="hbtp!",
                  float_format=2)
    pre_table = []
    for index, method in enumerate(methods[:]):
        B = []
        for mol_name in mol_names:
            for noise in noises:
                datas = get_attrs(get_file_name(mol_name, noise, method), ("ref_ener", "energy", "prob", "gate_count", "addition_res:"))
                if ne is not None:
                    ar = np.array(datas["addition_res:"])
                    n2 = ar[:,1] - 2*ne*ar[:,0] + ne*ne
                    liney, linex = np.sqrt(n2), np.array(datas["prob"])
                else:
                    liney, linex = np.array(datas["energy"]), np.array(datas["prob"])

                if noise in {"ion", "sc"}:
                    pass
                else:
                    if noise == "D":
                        linex = 1 - np.sqrt(1 - 3/4*(1 - linex))
                    else:
                        linex = 1 - np.sqrt(1 - (1 - linex))
                liney = liney[:]
                linex = linex[:]
                # linex = np.log(linex[3:])
                # liney = np.log(liney[3:])

                # A,b = np.polyfit(linex, liney, 1)
                result = linregress(linex, liney)
                B.append(round_with_uncertainty(result.slope, result.stderr))
                # B.append(round_with_uncertainty(result.intercept, result.intercept_stderr))

        pre_table.append([legends[index], *B, datas["gate_count"][0]["cx"]])
        # pre_table.append([legends[index], *B, gates[index]])

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

    print(table.generate_latex())

def plot_cnot_depth(lines, qubits, depth, cnot, pauli_num):
    ax_depth = Axis(lines)
    dlines = copy(lines)
    ax_cnot = Axis(dlines)
    for i, line in enumerate(lines):
        lines[i].X = qubits
        lines[i].Y = depth[i]
        # ax_depth.add_line(lines[i])
        dlines[i].X = qubits
        dlines[i].Y = np.array(cnot[i])/np.array(pauli_num)
        # ax_cnot.add_line(dlines[i])
    # Axis.add_line()
    gd = GraphData([ax_depth, ax_cnot])
    print(gd.generate_tikz())
    
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
    legends = ["JW", "BK","JW GdBM", "BK GdBM",  "FSN a-t-a", 
            #    "FSN $2\\times N$", 
               "MSN $2\\times N$", "FSN $2\\times N$"][4:]
    mol_name = "datah2h2ExcSolProb/DH4" 
    mol_name = "datah2h2noise_level/DH4" 
    mol_name = "data/adapt_vqe_h2_8/"
    mol_names = ["data_last/LiH_10", "data_last/LiH_10", "data_last/LiH_10_2xn_"]
    # mol_name = "data_all-to-all/LiH_8"
    # gates = ["49", "35", "21", "15", "15", "21", "16"]
    # gates = ["417", "394", "151", "132", "134", "170", "96"]

    # mol_names = ["data_all-to-all/H2_4", "data_all-to-all/H2_8", "data_all-to-all/LiH_8"]
    # mol_names = ["data_last/H2_4", "data_last/H2_8", "data_last/LiH_8"]
    # -1.8716797649675656
    # mol_name, ref_en = ("datah2_4/H2_4", -1.8573730129353947)
    # mol_name, ref_en = ("data02/H2_8", -1.8716797649675656)
    noises = ["sc"]

    # noises = ["ion"]
    # noises = ["X","Z"]
    noises = ["D","X","Y","Z"]
    methods = Circuits.get_circs_names()[:4] + Circuits.get_circs_names()[5:6] + Circuits.get_circs_names()[5:6] + Circuits.get_circs_names()[4:5] 
    methods =  Circuits.get_circs_names()[4:6] + Circuits.get_circs_names()[5:6]
    axes = []
    # plot_table(mol_name=mol_name, gates=None, ne=None)
    # for num, noise in enumerate(noises):
        # lines = copy([jw_line, bk_line, jw_opt_line, bk_opt_line,  swapgens_line, swapgeny_line, swap2xnshort_line, swap2xnyor_line][:])
        # energy_error_add_axes()
        # dn_add_axes(num=1, lines=[jw_line, bk_line, jw_opt_line, bk_opt_line, swapgeny_line, swap2xnshort_line])
    cnot = [[413, 1340, 3087, 5910, 10065, 15808, 23395, 33082, 45125, 59780, 77303],
            [394, 1528, 2881, 6064, 9546, 13594, 17771, 25778, 34489, 43765, 53287],
            [127, 351, 737, 1237, 1894, 2865, 3898, 5451, 6605, 8906, 10882],
            [112, 372, 748, 1301, 2168, 2936, 4156, 5915, 7547, 10021, 12209],
            [96, 240, 448, 720, 1056, 1456, 1920, 2448, 3040, 3696, 4416],
            [170, 411, 756, 1205, 1758, 2415, 3176, 4041, 5010, 6083, 7260],
            [134, 321, 588, 935, 1362, 1869, 2456, 3123, 3870, 4697, 5604]
            ]
    depth = [[527, 1550, 3346, 6133, 10103, 15478, 22457, 31244, 42070, 55093, 70602],
             [497, 1840, 3257, 7118, 10723, 15592, 18397, 29528, 38484, 48184, 57396],
             [177, 437, 802, 1267, 1871, 2694, 3639, 4929, 5641, 7571, 8939],
             [148, 458, 817, 1348, 2118, 2791, 3968, 5355, 6373, 8249, 10319],
             [65, 97, 129, 161, 193, 225, 257, 289, 321, 353, 385],
             [141, 211, 281, 351, 421, 491, 561, 631, 701, 771, 841],
             [125, 187, 249, 311, 373, 435, 497, 559, 621, 683, 745]]
    pauli_num = [72.0, 180.0, 336.0, 540.0, 792.0, 1092.0, 1440.0, 1836.0, 2280.0, 2772.0, 3312.0]
    qubits = [2*i for i in range(4,26,2)]
    datas = get_attrs(get_file_name("data_last/LiH_10_", "ion", "swap 2xn"), ("ref_ener", "energy_array"))
    re, ea = datas["ref_ener"][0], datas["energy_array"][0]
    ion_line = Line(list(range(len(ea))), np.array(ea) - re, "black", "L_BFGS_B")
    print(GraphData([Axis([ion_line])]).generate_tikz())

    # plot_cnot_depth([jw_line, bk_line, jw_opt_line, bk_opt_line,swap2xnshort_line, swapgeny_line,  swapgens_line], qubits, depth, cnot, pauli_num)
    # gd = GraphData(axes)
    # print(gd.generate_tikz())