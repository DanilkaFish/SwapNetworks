import argparse
import logging
from pathlib import Path

import numpy as np

import main as experiments
from Ternary_Tree.qiskit_interface.circuit_provider import Circuits, get_file_name


LOGGER = logging.getLogger(__name__)

PAULI_TASKS = (
    ("D", Circuits.jw()),
    ("D", Circuits.bk()),
    ("X", Circuits.jw()),
    ("X", Circuits.bk()),
    ("Y", Circuits.jw()),
    ("Y", Circuits.bk()),
    ("Z", Circuits.jw()),
    ("Z", Circuits.bk()),
)
ION_NETWORK_TASKS = (
    ("ion_scheduled", Circuits.swap_gen_yor()),
    ("ion_scheduled", Circuits.swap_gen_yor_2xn()),
    ("ion_scheduled", Circuits.swap_2xn()),
)
TASK_GROUPS = {
    "pauli-lih12-jw-bk": PAULI_TASKS,
    "ion-lih12-networks": ION_NETWORK_TASKS,
}


def parse_float_array(values, default):
    if values is None:
        return default
    return np.array([float(value) for value in values])


def configure_experiment_globals(args):
    experiments.NOISY_INITIAL_POINT_MODE = args.initial_point_mode
    experiments.NOISELESS_VQE_RESTARTS = args.noiseless_restarts
    experiments.NOISY_VQE_RESTARTS = args.noisy_restarts
    experiments.CACHE_NOISELESS_INITIAL_POINTS = not args.no_noiseless_cache


def output_prefix(args, noise):
    if noise.startswith("ion"):
        return str(Path(args.scheduled_output_dir) / "LiH_12")
    return str(Path(args.pauli_output_dir) / "LIH12")


def probabilities(args, noise):
    if noise.startswith("ion"):
        return parse_float_array(args.ion_multipliers, experiments.SC_ION_MULTIPLIERS)
    return parse_float_array(args.pauli_probs, experiments.PAULI_PROBS)


def result_path(args, noise, circuit):
    return Path(get_file_name(output_prefix(args, noise), noise, circuit))


def expanded_tasks(groups):
    tasks = []
    for group in groups:
        tasks.extend(TASK_GROUPS[group])
    return tasks


def print_tasks(args, tasks):
    for index, (noise, circuit) in enumerate(tasks, start=1):
        path = result_path(args, noise, circuit)
        state = "exists" if path.exists() else "new"
        print(f"{index:02d}: LiH_12 {noise} {circuit} -> {path} [{state}]")


def run_tasks(args):
    configure_experiment_globals(args)
    tasks = expanded_tasks(args.groups)

    if args.dry_run:
        print_tasks(args, tasks)
        return 0

    current_noise = None
    current_prefix = None
    vqe_data = None
    failures = []

    for task_index, (noise, circuit) in enumerate(tasks, start=1):
        path = result_path(args, noise, circuit)
        if path.exists() and not args.overwrite:
            LOGGER.info("Skipping existing result: %s", path)
            continue

        prefix = output_prefix(args, noise)
        if (noise, prefix) != (current_noise, current_prefix):
            LOGGER.info("Preparing LiH_12 with %s", noise)
            vqe_data = experiments.VQEData(
                prefix,
                experiments.LiH_12,
                experiments.OPTIMIZER,
                reps=args.reps,
                probs=probabilities(args, noise),
                noise_type=noise,
                device=args.device,
            )
            current_noise = noise
            current_prefix = prefix

        LOGGER.info("Task %s/%s: %s %s", task_index, len(tasks), noise, circuit)
        try:
            experiments.run_vqe(circuit, vqe_data, args.distance)
        except Exception as exc:
            failures.append((noise, circuit, exc))
            LOGGER.exception("Failed: %s %s", noise, circuit)
            if not args.keep_going:
                raise

    if failures:
        LOGGER.error("%s task(s) failed", len(failures))
        for noise, circuit, exc in failures:
            LOGGER.error("%s %s: %s", noise, circuit, exc)
        return 1
    return 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Run only the LiH 12q experiments needed for currently missing table entries."
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        choices=sorted(TASK_GROUPS),
        default=sorted(TASK_GROUPS),
        help="Task groups to run.",
    )
    parser.add_argument("--scheduled-output-dir", default="data_scheduled")
    parser.add_argument("--pauli-output-dir", default="data_revised")
    parser.add_argument("--device", default=experiments.DEVICE)
    parser.add_argument("--reps", type=int, default=experiments.REPS)
    parser.add_argument("--distance", type=float, default=experiments.DISTANCE)
    parser.add_argument("--ion-multipliers", nargs="+", default=None)
    parser.add_argument("--pauli-probs", nargs="+", default=None)
    parser.add_argument(
        "--initial-point-mode",
        choices=sorted(experiments.INITIAL_POINT_MODES),
        default=experiments.NOISY_INITIAL_POINT_MODE,
    )
    parser.add_argument("--noiseless-restarts", type=int, default=experiments.NOISELESS_VQE_RESTARTS)
    parser.add_argument("--noisy-restarts", type=int, default=experiments.NOISY_VQE_RESTARTS)
    parser.add_argument("--no-noiseless-cache", action="store_true")
    parser.add_argument("--overwrite", action="store_true", help="Recompute files that already exist.")
    parser.add_argument("--keep-going", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))
    return run_tasks(args)


if __name__ == "__main__":
    raise SystemExit(main())
