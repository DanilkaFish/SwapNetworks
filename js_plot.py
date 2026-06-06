import argparse
import copy
import json
import math
from pathlib import Path

from utils import Axis, AxisScale, GraphData, Line, Marks, Table


JW = "jw"
BK = "bk"
JW_LEX = "jw_lex"
BK_LEX = "bk_lex"
SWAP_GEN_YOR = "swap gen yor"
SWAP_2XN = "swap 2xn"

ACTIVE_METHODS = (JW, BK, JW_LEX, BK_LEX, SWAP_GEN_YOR, SWAP_2XN)
ACTIVE_LEGENDS = ("JW", "BK", "JW GdBM", "BK GdBM", "FSN a-t-a", "MSN $2\\times N$")
ACTIVE_MOLECULES = ("data_last/H2_4", "data_last/H2_8", "data_last/LiH_8")
ACTIVE_MOLECULE_LABELS = ("H2 4q", "H2 8q", "LiH 8q")
ACTIVE_NOISES = ("ion",)

CNOT_DEPTH_LEGENDS = ("JW", "BK", "JWO", "BKO", "MSN S", "FSN Y", "FSN S")

CNOT_COUNTS = (
    (413, 1340, 3087, 5910, 10065, 15808, 23395, 33082, 45125, 59780, 77303),
    (394, 1528, 2881, 6064, 9546, 13594, 17771, 25778, 34489, 43765, 53287),
    (127, 351, 737, 1237, 1894, 2865, 3898, 5451, 6605, 8906, 10882),
    (112, 372, 748, 1301, 2168, 2936, 4156, 5915, 7547, 10021, 12209),
    (96, 240, 448, 720, 1056, 1456, 1920, 2448, 3040, 3696, 4416),
    (170, 411, 756, 1205, 1758, 2415, 3176, 4041, 5010, 6083, 7260),
    (134, 321, 588, 935, 1362, 1869, 2456, 3123, 3870, 4697, 5604),
)
CNOT_DEPTHS = (
    (527, 1550, 3346, 6133, 10103, 15478, 22457, 31244, 42070, 55093, 70602),
    (497, 1840, 3257, 7118, 10723, 15592, 18397, 29528, 38484, 48184, 57396),
    (177, 437, 802, 1267, 1871, 2694, 3639, 4929, 5641, 7571, 8939),
    (148, 458, 817, 1348, 2118, 2791, 3968, 5355, 6373, 8249, 10319),
    (65, 97, 129, 161, 193, 225, 257, 289, 321, 353, 385),
    (141, 211, 281, 351, 421, 491, 561, 631, 701, 771, 841),
    (125, 187, 249, 311, 373, 435, 497, 559, 621, 683, 745),
)
PAULI_COUNTS = (72.0, 180.0, 336.0, 540.0, 792.0, 1092.0, 1440.0, 1836.0, 2280.0, 2772.0, 3312.0)
QUBITS = tuple(2 * i for i in range(4, 26, 2))


def result_path(name: str, noise: str, method: str) -> Path:
    return Path(f"{name}_{noise}{method}.json")


def get_attrs(file_name: Path, attrs: tuple[str, ...]):
    data = {name: [] for name in attrs}
    with file_name.open() as result_file:
        for item in json.load(result_file):
            for attr in attrs:
                data[attr].append(item[attr])
    return data


def regression_slope(x_values, y_values):
    x_values = [float(value) for value in x_values]
    y_values = [float(value) for value in y_values]
    x_mean = sum(x_values) / len(x_values)
    y_mean = sum(y_values) / len(y_values)
    denominator = sum((value - x_mean) ** 2 for value in x_values)
    if denominator == 0:
        return 0
    return sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values)) / denominator


def noise_axis(noise: str, probs):
    if noise in {"ion", "sc"}:
        return [float(prob) * 1.5 for prob in probs]
    if noise == "D":
        return [0.75 * (1 - float(prob)) for prob in probs]
    return [1 - float(prob) for prob in probs]


def plot_axis(noise: str, probs):
    if noise == "sc":
        return [float(prob) for prob in probs]
    if noise == "D":
        return [1 - math.sqrt(1 - 0.75 * (1 - float(prob))) for prob in probs]
    return [1 - math.sqrt(float(prob)) for prob in probs]


def round_slope(value):
    return f"{value:.3f}"


def make_method_lines():
    styles = (
        (Marks.triangle(), "blue"),
        (Marks.triangle(), "red"),
        (Marks.square(), "blue"),
        (Marks.square(), "red"),
        (Marks.pentagon(), "green"),
        (Marks.pentagon(), "black"),
    )
    return [
        Line(None, None, mark, color, legend)
        for (mark, color), legend in zip(styles, ACTIVE_LEGENDS)
    ]


def plot_ion_table(ne=None):
    table = Table(
        "Noise susceptibility $\\chi$ for the active ion-noise data.",
        ["Method", *ACTIVE_MOLECULE_LABELS],
        align=[5, *[4 for _ in ACTIVE_MOLECULES]],
        bordered=False,
        placement="hbtp!",
    )
    rows = []

    for method_index, method in enumerate(ACTIVE_METHODS):
        slopes = []
        for molecule in ACTIVE_MOLECULES:
            for noise in ACTIVE_NOISES:
                attrs = ("ref_ener", "energy", "prob", "gate_count", "addition_res:")
                data = get_attrs(result_path(molecule, noise, method), attrs)
                if ne is None:
                    y_values = [float(value) for value in data["energy"]]
                else:
                    y_values = []
                    for additional in data["addition_res:"]:
                        variance = additional[1] - 2 * ne * additional[0] + ne * ne
                        y_values.append(math.sqrt(variance))

                x_values = noise_axis(noise, data["prob"])
                slopes.append(round_slope(regression_slope(x_values[:3], y_values[:3])))
        rows.append([ACTIVE_LEGENDS[method_index], *slopes])

    for column in range(1, len(rows[0])):
        best_index = min(range(len(rows)), key=lambda index: float(rows[index][column]))
        for row in rows:
            row[column] = f"${row[column]}$"
        rows[best_index][column] = f"$\\mathbf{{{rows[best_index][column][1:-1]}}}$"

    for row in rows:
        table.add_row(*row)
    print(table.generate_latex())


def plot_energy_error():
    axes = []
    for noise_index, noise in enumerate(ACTIVE_NOISES):
        lines = make_method_lines()
        for method_index, method in enumerate(ACTIVE_METHODS):
            data = get_attrs(result_path("data_last/H2_8", noise, method), ("ref_ener", "energy", "prob", "gate_count"))
            lines[method_index].X = plot_axis(noise, data["prob"])
            lines[method_index].Y = [
                -abs(float(energy) - float(ref)) / float(ref)
                for energy, ref in zip(data["energy"], data["ref_ener"])
            ]
            cx_count = data["gate_count"][0].get("cx", 0)
            lines[method_index].legend = f"{ACTIVE_LEGENDS[method_index]} cx={cx_count}"
        axes.append(Axis(lines, "$\\lambda$", "$|E-E_0|$", AxisScale.usual(), xshift=noise_index * 8, title=noise))
    print(GraphData(axes).generate_tikz())


def plot_number_variance(ne=2):
    axes = []
    for noise_index, noise in enumerate(ACTIVE_NOISES):
        lines = make_method_lines()
        for method_index, method in enumerate(ACTIVE_METHODS):
            data = get_attrs(result_path("data_last/H2_8", noise, method), ("prob", "gate_count", "addition_res:"))
            n2 = [
                additional[1] - 2 * ne * additional[0] + ne * ne
                for additional in data["addition_res:"]
            ]
            lines[method_index].X = plot_axis(noise, data["prob"])
            lines[method_index].Y = n2
            cx_count = data["gate_count"][0].get("cx", 0)
            lines[method_index].legend = f"{ACTIVE_LEGENDS[method_index]} cx={cx_count}"
        axes.append(
            Axis(
                lines,
                "$1 -\\mathcal{F}$",
                f"$\\langle(n-{ne})^2\\rangle$",
                AxisScale.loglog(),
                xshift=noise_index * 8,
                title=noise,
            )
        )
    print(GraphData(axes).generate_tikz())


def plot_cnot_depth():
    base_lines = [
        Line(None, None, Marks.triangle(), "blue", "JW"),
        Line(None, None, Marks.triangle(), "red", "BK"),
        Line(None, None, Marks.square(), "blue", "JWO"),
        Line(None, None, Marks.square(), "red", "BKO"),
        Line(None, None, Marks.pentagon(), "black, dashed", "MSN S"),
        Line(None, None, Marks.pentagon(), "green", "FSN Y"),
        Line(None, None, Marks.pentagon(), "orange", "FSN S"),
    ]
    depth_lines = copy.deepcopy(base_lines)
    cnot_lines = copy.deepcopy(base_lines)

    for index, line in enumerate(depth_lines):
        line.X = QUBITS
        line.Y = CNOT_DEPTHS[index]
    for index, line in enumerate(cnot_lines):
        line.X = QUBITS
        line.Y = [count / pauli_count for count, pauli_count in zip(CNOT_COUNTS[index], PAULI_COUNTS)]

    depth_axis = Axis(depth_lines, "$N$", "depth", AxisScale.usual(), title="Depth")
    cnot_axis = Axis(cnot_lines, "$N$", "CX per Pauli", AxisScale.usual(), yshift=-8, title="CX count")
    print(GraphData([depth_axis, cnot_axis]).generate_tikz())


def plot_optimizer_convergence():
    data = get_attrs(result_path("data_last/LiH_10_", "ion", SWAP_2XN), ("ref_ener", "energy_array"))
    ref_energy = float(data["ref_ener"][0])
    energy_array = [float(value) for value in data["energy_array"][0]]
    line = Line(range(len(energy_array)), [energy - ref_energy for energy in energy_array], color="black", legend="L_BFGS_B")
    print(GraphData([Axis([line], "step", "$E - E_0$", AxisScale.usual(), title="Optimizer convergence")]).generate_tikz())


def main(argv=None):
    parser = argparse.ArgumentParser(description="Print LaTeX/TikZ plots for experiment data.")
    parser.add_argument(
        "preset",
        choices=("ion-table", "energy-error", "number-variance", "cnot-depth", "optimizer-convergence"),
        help="Plot/table preset to print.",
    )
    args = parser.parse_args(argv)

    if args.preset == "ion-table":
        plot_ion_table()
    elif args.preset == "energy-error":
        plot_energy_error()
    elif args.preset == "number-variance":
        plot_number_variance()
    elif args.preset == "cnot-depth":
        plot_cnot_depth()
    elif args.preset == "optimizer-convergence":
        plot_optimizer_convergence()


if __name__ == "__main__":
    main()
