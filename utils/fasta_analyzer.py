"""
FASTA sequence analysis utilities using biotite.

Provides functions for parsing FASTA files, detecting sequence types,
and calculating statistics like sequence lengths and residue frequencies.
"""
from io import StringIO
from typing import List, Tuple, Dict, Any
from collections import Counter

import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta


def detect_sequence_type(sequence: str) -> str:
    """
    Auto-detect sequence type using biotite's alphabet validation.

    Args:
        sequence: The sequence string to analyze

    Returns:
        One of: "protein", "dna", "rna"
    """
    seq_upper = sequence.upper()

    # Try nucleotide first (biotite's approach)
    try:
        seq.NucleotideSequence(seq_upper)
        return "dna"
    except seq.AlphabetError:
        pass

    # Check for RNA (U instead of T)
    if "U" in seq_upper and "T" not in seq_upper:
        try:
            seq.NucleotideSequence(seq_upper.replace("U", "T"))
            return "rna"
        except seq.AlphabetError:
            pass

    # Fall back to protein
    return "protein"


def parse_fasta_biotite(fasta_text: str) -> List[Tuple[str, str, str]]:
    """
    Parse FASTA using biotite's FastaFile.

    Args:
        fasta_text: Raw FASTA text content

    Returns:
        List of tuples: (header, sequence, detected_seq_type)
    """
    fasta_file = fasta.FastaFile.read(StringIO(fasta_text))
    sequences = []

    for header, sequence_str in fasta_file.items():
        seq_type = detect_sequence_type(sequence_str)
        sequences.append((header, sequence_str, seq_type))

    return sequences


def calculate_sequence_lengths(sequences: List[Tuple[str, str, str]]) -> Dict[str, Any]:
    """
    Calculate length statistics for all sequences.

    Args:
        sequences: List of (header, sequence, seq_type) tuples

    Returns:
        Dictionary with length statistics
    """
    if not sequences:
        return {
            "lengths": [],
            "headers": [],
            "min": 0,
            "max": 0,
            "mean": 0,
            "median": 0,
            "total_residues": 0,
        }

    lengths = [len(sequence) for _, sequence, _ in sequences]
    headers = [header for header, _, _ in sequences]
    sorted_lengths = sorted(lengths)

    return {
        "lengths": lengths,
        "headers": headers,
        "min": min(lengths),
        "max": max(lengths),
        "mean": sum(lengths) / len(lengths),
        "median": sorted_lengths[len(sorted_lengths) // 2],
        "total_residues": sum(lengths),
    }


def calculate_residue_frequencies(
    sequences: List[Tuple[str, str, str]], seq_type: str
) -> Dict[str, Any]:
    """
    Calculate residue frequencies for all sequences.

    Args:
        sequences: List of (header, sequence, seq_type) tuples
        seq_type: Sequence type ("protein", "dna", "rna", or "auto")

    Returns:
        Dictionary with residue counts and percentages
    """
    # Get alphabets from biotite (use unambiguous alphabets to avoid overlap with amino acids)
    if seq_type == "protein":
        expected_residues = "".join(seq.ProteinSequence.alphabet.get_symbols())
    elif seq_type == "dna":
        expected_residues = "ACGTN"  # Unambiguous DNA + N for unknown
    elif seq_type == "rna":
        expected_residues = "ACGUN"  # Unambiguous RNA + N for unknown
    else:
        # Auto-detect based on majority of sequences
        type_counts = Counter(s_type for _, _, s_type in sequences)
        majority_type = type_counts.most_common(1)[0][0] if type_counts else "protein"
        return calculate_residue_frequencies(sequences, majority_type)

    # Count all residues
    all_residues = "".join(sequence.upper() for _, sequence, _ in sequences)
    residue_counts = Counter(all_residues)
    total = len(all_residues)

    # Build results with expected residues
    results = {
        "counts": {},
        "percentages": {},
        "total": total,
        "seq_type": seq_type,
    }

    for residue in expected_residues:
        count = residue_counts.get(residue, 0)
        results["counts"][residue] = count
        results["percentages"][residue] = (count / total * 100) if total > 0 else 0

    # Add any unexpected residues found
    for residue, count in residue_counts.items():
        if residue not in expected_residues and residue.isalpha():
            results["counts"][residue] = count
            results["percentages"][residue] = (count / total * 100) if total > 0 else 0

    return results


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC content for a nucleotide sequence.

    Args:
        sequence: Nucleotide sequence string

    Returns:
        GC content as a percentage (0-100)
    """
    seq_upper = sequence.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    total = sum(1 for c in seq_upper if c in "ACGTU")
    return (gc_count / total * 100) if total > 0 else 0


def get_per_sequence_stats(
    sequences: List[Tuple[str, str, str]]
) -> List[Dict[str, Any]]:
    """
    Calculate statistics for each individual sequence.

    Args:
        sequences: List of (header, sequence, seq_type) tuples

    Returns:
        List of dictionaries with per-sequence statistics
    """
    stats = []
    for header, sequence, seq_type in sequences:
        seq_stats = {
            "header": header,
            "length": len(sequence),
            "type": seq_type,
        }

        # Add GC content for nucleotide sequences
        if seq_type in ("dna", "rna"):
            seq_stats["gc_content"] = calculate_gc_content(sequence)

        stats.append(seq_stats)

    return stats


def analyze_fasta(
    fasta_text: str, seq_type_override: str = "auto"
) -> Dict[str, Any]:
    """
    Main analysis function returning complete results.

    Args:
        fasta_text: Raw FASTA text content
        seq_type_override: Override sequence type detection ("auto", "protein", "dna", "rna")

    Returns:
        Dictionary containing all analysis results
    """
    # Parse FASTA
    sequences = parse_fasta_biotite(fasta_text)

    if not sequences:
        return {
            "success": False,
            "error": "No valid sequences found in FASTA input",
            "sequences": [],
            "length_stats": {},
            "residue_frequencies": {},
            "per_sequence_stats": [],
        }

    # Override sequence types if specified
    if seq_type_override != "auto":
        sequences = [(h, s, seq_type_override) for h, s, _ in sequences]

    # Determine the overall sequence type for analysis
    type_counts = Counter(s_type for _, _, s_type in sequences)
    majority_type = type_counts.most_common(1)[0][0]

    # Calculate statistics
    length_stats = calculate_sequence_lengths(sequences)
    residue_frequencies = calculate_residue_frequencies(sequences, majority_type)
    per_sequence_stats = get_per_sequence_stats(sequences)

    return {
        "success": True,
        "error": None,
        "total_sequences": len(sequences),
        "detected_type": majority_type,
        "type_distribution": dict(type_counts),
        "sequences": sequences,
        "length_stats": length_stats,
        "residue_frequencies": residue_frequencies,
        "per_sequence_stats": per_sequence_stats,
    }
