from my_utils import GraphData, Axis, Marks, AxisScale, Line
from numpy import array
jw2 = [413, 1340, 3087, 5910, 10065, 15808, 23395, 33082, 45125, 59780, 77303, 97950]
jw1 = [370, 963, 1832, 2977, 4398, 6095, 8068, 10317, 12842, 15643, 18720, 22073]
jwd = [527, 1550, 3346, 6133, 10103, 15478, 22457, 31244, 42070, 55093, 70602, 88710]
bk2 = [394, 1528, 2881, 6064, 9546, 13594, 17771, 25778, 34489, 43765, 53287, 64024]
bk1 = [298, 876, 1583, 2700, 4046, 5623, 7222, 9428, 11922, 14560, 17523, 20663]
bkd = [497, 1840, 3257, 7118, 10723, 15592, 18397, 29528, 38484, 48184, 57396, 70769]
jwr1 = [149, 392, 819, 1339, 2004, 2994, 3996, 5630, 6903, 9276, 11351, 14902]
jwr2 = [181, 416, 888, 1403, 2071, 3021, 4063, 5523, 6820, 8864, 10597, 13873]
jwrd = [204, 477, 882, 1350, 1980, 2811, 3750, 5097, 5885, 7876, 9318, 12194]
bkr2 = [126, 407, 804, 1406, 2269, 3109, 4413, 6111, 7957, 10523, 12695, 19096]
bkr1 = [148, 428, 784, 1458, 2110, 3079, 4194, 6062, 7234, 9845, 11714, 17083]
bkrd = [164, 498, 866, 1439, 2223, 2911, 4158, 5512, 6684, 8643, 10680, 15599]
s2xn2 = [100, 246, 456, 730, 1068, 1470, 1936, 2466, 3060, 3718, 4440, 5226]
s2xn1 = [174, 420, 772, 1230, 1794, 2464, 3240, 4122, 5110, 6204, 7404, 8710]
s2xnd = [67, 99, 131, 163, 195, 227, 259, 291, 323, 355, 387, 419]
sgen_s2 = [142, 342, 626, 994, 1446, 1982, 2602, 3306, 4094, 4966, 5922, 6962]
sgen_s1 = [240, 570, 1040, 1650, 2400, 3290, 4320, 5490, 6800, 8250, 9840, 11570]
sgen_sd = [99, 152, 205, 258, 311, 364, 417, 470, 523, 576, 629, 682]
sgen_yr2 = [136, 327, 598, 949, 1380, 1891, 2482, 3153, 3904, 4735, 5646, 6637]
sgen_yr1 = [260, 618, 1128, 1790, 2604, 3570, 4688, 5958, 7380, 8954, 10680, 12558]
sgen_yrd = [126, 192, 258, 324, 390, 456, 522, 588, 654, 720, 786, 852]
pauli_num = [60.0, 150.0, 280.0, 450.0, 660.0, 910.0, 1200.0, 1530.0, 1900.0, 2310.0, 2760.0, 3250.0]
n = len(jw2)
Xs = [4*i for i in range(n)]
jw_line = Line(None, None, Marks.triangle(), "blue", "ДВ")
bk_line = Line(None, None, Marks.triangle(), "red", "БК")
jw_opt_line = Line(None, None, Marks.square(), "blue", "ДВ опт")
bk_opt_line = Line(None, None, Marks.square(), "red", "БК опт")
swap2xn_line = Line(None, None, Marks.pentagon(), "green", "SWAP 2xn")
swapgens_line = Line(None, None, Marks.pentagon(), "black", "SWAP short")
swapgeny_line = Line(None, None, Marks.pentagon(), "gray", "SWAP yor")

if __name__ == "__main__":
    lines=[jw_line, bk_line, jw_opt_line, bk_opt_line, swap2xn_line, swapgens_line, swapgeny_line],
    depth = Axis(
        lines=lines,
        xlabel= "Число кубитов",
        ylabel= "Глубина схемы",
        axistype=AxisScale.usual(),
        )
    for line in lines:
        line.Y = line.Y/array(pauli_num)
    cx_num = Axis(
        lines=lines,
        xlabel= "Число кубитов",
        ylabel= "Глубина схемы",
        axistype=AxisScale.usual(),
        )
    # Xs = 
    gd = GraphData([cx_num])
    print(gd.generate_tikz())