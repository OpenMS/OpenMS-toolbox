"""Standalone FASTA analysis runner for subprocess execution."""
import sys
import json
from pathlib import Path


def main():
    input_path, output_path = sys.argv[1], sys.argv[2]

    # Set up imports (project root = parent of utils/)
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from utils.fasta_analyzer import analyze_fasta

    fasta_text = Path(input_path).read_text(encoding="utf-8")
    results = analyze_fasta(fasta_text, "auto")

    # Drop fields not needed by UI (not JSON-serializable / large)
    results.pop("sequences", None)
    results.pop("per_sequence_stats", None)

    Path(output_path).write_text(json.dumps(results), encoding="utf-8")


if __name__ == "__main__":
    main()
