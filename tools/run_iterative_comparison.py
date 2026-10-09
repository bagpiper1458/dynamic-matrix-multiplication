"""Run sequential, unprofiled timings; save raw samples and source provenance."""
import argparse
import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--profile", choices=("historical", "packed"), default="historical")
    parser.add_argument("--max-exponent", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--iterations", type=int, default=20)
    args = parser.parse_args()
    packed = args.profile == "packed"
    executable = args.executable or ROOT / ("build/dm_compare_packed.exe" if packed else "build/dm_compare.exe")
    destination = ROOT / ("analysis/packed" if packed else "analysis/iterative")
    destination.mkdir(parents=True, exist_ok=True)
    metadata = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "processor": os.environ.get("PROCESSOR_IDENTIFIER", platform.processor()),
        "compile_flags": "-std=c++17 -O3 -Wall -Wextra -Wpedantic",
        "note": "Local run; fixed matrix per size; " + ("16 rotating vectors" if packed else "fixed vector") + "; no memory profiler during timing.",
        "source_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in (("cpp/dynamic_matrix.hpp", "cpp/iterative_matrix.hpp", "cpp/packed_matrix.hpp", "cpp/compare_packed.cpp")
                                       if packed else ("cpp/dynamic_matrix.hpp", "cpp/iterative_matrix.hpp", "cpp/compare.cpp"))},
        "runs": [],
    }
    for kind in ("toeplitz", "circulant", "unit-circulant"):
        command = [str(executable.resolve()), str(args.max_exponent), str(args.repeats), str(args.iterations), kind]
        run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=True)
        (destination / f"{kind}-raw.csv").write_text(run.stdout, encoding="utf-8")
        # Keep public artifacts independent of the local username/workspace path.
        try:
            recorded_executable = executable.resolve().relative_to(ROOT).as_posix()
        except ValueError:
            recorded_executable = executable.name
        metadata["runs"].append({"command": [recorded_executable, *command[1:]], "stderr": run.stderr})
        print(f"Saved {kind} raw samples", flush=True)
    metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
    (destination / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
