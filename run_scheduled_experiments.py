import argparse
import logging
from pathlib import Path

import numpy as np

import main as experiments
from Ternary_Tree.qiskit_interface.circuit_provider import Circuits, get_file_name


LOGGER = logging.getLogger(__name__)

MOLECULES = {
    "h2_4": ("H2_4", experiments.H2_4),
    "h2_8": ("H2_8", experiments.H2_8),
    "lih_8": ("LiH_8", experiments.LiH_8),
    "lih_10": ("LiH_10", experiments.LiH_10),
    "lih_12": ("LiH_12", experiments.LiH_12),
}

MOLECULE_PRESETS = {
    "all": ("h2_4", "h2_8", "lih_8", "lih_10", "lih_12"),
    "h2": ("h2_4", "h2_8"),
    "lih": ("lih_8", "lih_10", "lih_12"),
}

CIRCUITS = {
    "jw": Circuits.jw(),
    "bk": Circuits.bk(),
    "jw_lex": Circuits.jw_lex(),
    "bk_lex": Circuits.bk_lex(),
    "swap_gen_yor": Circuits.swap_gen_yor(),
    "swap_gen_yor_2xn": Circuits.swap_gen_yor_2xn(),
    "swap_2xn": Circuits.swap_2xn(),
    "swap_2xn_yor": Circuits.swap_2xn_yor(),
}

CIRCUIT_PRESETS = {
    "table": ("jw", "bk", "jw_lex", "bk_lex", "swap_gen_yor", "swap_2xn"),
    "scheduled": (
        "jw",
        "bk",
        "jw_lex",
        "bk_lex",
        "swap_gen_yor",
        "swap_gen_yor_2xn",
        "swap_2xn",
        "swap_2xn_yor",
    ),
    "optimized": ("jw_lex", "bk_lex", "swap_gen_yor", "swap_gen_yor_2xn", "swap_2xn", "swap_2xn_yor"),
}

NOISES = {
    "sc_scheduled": "sc_scheduled",
    "ion_scheduled": "ion_scheduled",
}

NOISE_PRESETS = {
    "scheduled": ("sc_scheduled", "ion_scheduled"),
}


def expand_names(values, items, presets, kind):
    expanded = []
    for value in values:
        key = value.lower().replace("-", "_").replace(" ", "_")
        if key in presets:
            expanded.extend(presets[key])
        elif key in items:
            expanded.append(key)
        else:
            valid = sorted(list(items) + list(presets))
            raise ValueError(f"Unknown {kind} {value!r}; expected one of {valid}")

    result = []
    seen = set()
    for key in expanded:
        if key not in seen:
            result.append(key)
            seen.add(key)
    return result


def parse_multipliers(values):
    if values is None:
        return experiments.SC_ION_MULTIPLIERS
    return np.array([float(value) for value in values])


def result_file(output_dir, molecule_key, noise_key, circuit_key):
    molecule_name, _ = MOLECULES[molecule_key]
    prefix = str(Path(output_dir) / molecule_name)
    return Path(get_file_name(prefix, NOISES[noise_key], CIRCUITS[circuit_key]))


def make_tasks(args):
    if args.start_at < 1:
        raise ValueError("--start-at must be at least 1")
    if args.limit is not None and args.limit < 1:
        raise ValueError("--limit must be at least 1")

    molecule_keys = expand_names(args.molecules, MOLECULES, MOLECULE_PRESETS, "molecule")
    circuit_keys = expand_names(args.circuits, CIRCUITS, CIRCUIT_PRESETS, "circuit")
    noise_keys = expand_names(args.noises, NOISES, NOISE_PRESETS, "noise")

    tasks = []
    for molecule_key in molecule_keys:
        for noise_key in noise_keys:
            for circuit_key in circuit_keys:
                tasks.append((molecule_key, noise_key, circuit_key))

    start = args.start_at - 1
    if args.limit is None:
        return tasks[start:]
    return tasks[start:start + args.limit]


def configure_experiment_globals(args):
    experiments.NOISY_INITIAL_POINT_MODE = args.initial_point_mode
    experiments.NOISELESS_VQE_RESTARTS = args.noiseless_restarts
    experiments.NOISY_VQE_RESTARTS = args.noisy_restarts
    experiments.CACHE_NOISELESS_INITIAL_POINTS = not args.no_noiseless_cache


def print_tasks(tasks, output_dir):
    for index, (molecule_key, noise_key, circuit_key) in enumerate(tasks, start=1):
        path = result_file(output_dir, molecule_key, noise_key, circuit_key)
        state = "exists" if path.exists() else "new"
        print(
            f"{index:03d}: {molecule_key} {NOISES[noise_key]} {CIRCUITS[circuit_key]} "
            f"-> {path} [{state}]"
        )


def run_tasks(args):
    tasks = make_tasks(args)
    configure_experiment_globals(args)
    probs = parse_multipliers(args.multipliers)

    if args.dry_run:
        print_tasks(tasks, args.output_dir)
        return 0

    current_key = None
    vqe_data = None
    failures = []

    for task_index, (molecule_key, noise_key, circuit_key) in enumerate(tasks, start=1):
        path = result_file(args.output_dir, molecule_key, noise_key, circuit_key)
        if args.skip_existing and path.exists():
            LOGGER.info("Skipping existing result: %s", path)
            continue

        molecule_name, molecule = MOLECULES[molecule_key]
        noise = NOISES[noise_key]
        circuit = CIRCUITS[circuit_key]
        prefix = str(Path(args.output_dir) / molecule_name)
        group_key = (molecule_key, noise_key)

        if current_key != group_key:
            LOGGER.info("Preparing %s with %s", molecule_name, noise)
            vqe_data = experiments.VQEData(
                prefix,
                molecule,
                experiments.OPTIMIZER,
                reps=args.reps,
                probs=probs,
                noise_type=noise,
                device=args.device,
            )
            current_key = group_key

        LOGGER.info("Task %s/%s: %s %s %s", task_index, len(tasks), molecule_name, noise, circuit)
        try:
            experiments.run_vqe(circuit, vqe_data, args.distance)
        except Exception as exc:
            failures.append((molecule_key, noise_key, circuit_key, exc))
            LOGGER.exception("Failed: %s %s %s", molecule_key, noise_key, circuit_key)
            if not args.keep_going:
                raise

    if failures:
        LOGGER.error("%s task(s) failed", len(failures))
        for molecule_key, noise_key, circuit_key, exc in failures:
            LOGGER.error("%s %s %s: %s", molecule_key, noise_key, circuit_key, exc)
        return 1
    return 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Recalculate 4-12 qubit molecule experiments with scheduled hardware noise."
    )
    parser.add_argument("--molecules", nargs="+", default=["all"], help="Molecules or preset: all, h2, lih.")
    parser.add_argument(
        "--circuits",
        nargs="+",
        default=["scheduled"],
        help="Circuits or preset: table, scheduled, optimized.",
    )
    parser.add_argument("--noises", nargs="+", default=["scheduled"], help="Noise modes: sc_scheduled, ion_scheduled.")
    parser.add_argument("--output-dir", default="data_scheduled")
    parser.add_argument("--device", default=experiments.DEVICE)
    parser.add_argument("--reps", type=int, default=experiments.REPS)
    parser.add_argument("--distance", type=float, default=experiments.DISTANCE)
    parser.add_argument("--multipliers", nargs="+", default=None)
    parser.add_argument("--initial-point-mode", choices=sorted(experiments.INITIAL_POINT_MODES), default="noiseless")
    parser.add_argument("--noiseless-restarts", type=int, default=experiments.NOISELESS_VQE_RESTARTS)
    parser.add_argument("--noisy-restarts", type=int, default=experiments.NOISY_VQE_RESTARTS)
    parser.add_argument("--no-noiseless-cache", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    parser.add_argument("--start-at", type=int, default=1)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))
    return run_tasks(args)


if __name__ == "__main__":
    raise SystemExit(main())
