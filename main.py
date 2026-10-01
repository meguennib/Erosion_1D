"""
Entry point for the thesis reference solver (scientific freeze v1).
"""

from pathlib import Path

from src.config import Params
from src.simulation import run_simulation
from src.utils import setup_logger


def main():
    p = Params()
    out_dir = Path("output_thesis_v1")
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(out_dir / "run.log")

    result = run_simulation(p, logger)
    logger.info(
        "Reference run completed: "
        f"stop_reason={result['run_manifest']['stop_reason']}, "
        f"t_final/ter={result['run_manifest']['t_final_over_ter']:.5f}"
    )
    return result


if __name__ == "__main__":
    main()
