"""
FASTA Analyzer Page

This module provides functionality for analyzing FASTA files to compute
sequence length histograms and residue frequency statistics. Supports
both protein and nucleotide (DNA/RNA) sequences.
"""

import gzip
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import sys
from pathlib import Path

# Add utils to path
sys.path.append(str(Path(__file__).parent.parent))

from src.common.common import page_setup
from utils.fasta_analyzer import analyze_fasta

params = page_setup()

def main():
    """Main function for the FASTA Analyzer page."""
    st.title("FASTA Analyzer")

    st.markdown("""
    **Analyze FASTA sequence files** to compute statistics about sequence lengths and residue composition.

    This tool helps you:
    - **Understand sequence datasets** by visualizing length distributions
    - **Analyze residue composition** with frequency statistics for amino acids or nucleotides
    - **Calculate GC content** for nucleotide sequences
    """)

    with st.expander("How FASTA Analysis Works"):
        st.markdown("""
        **Sequence Type Detection:**
        - Sequences are auto-detected as DNA, RNA, or protein based on their characters
        - Nucleotide sequences may contain standard bases, modified bases, and IUPAC ambiguity codes
        - Sequences with U (but no T) are classified as RNA; otherwise as DNA
        - Sequences containing characters outside the nucleotide alphabet are classified as protein

        **Supported Nucleotide Alphabet:**

        | Symbol | Meaning |
        |--------|---------|
        | A | Adenosine |
        | C | Cytidine |
        | G | Guanosine |
        | T | Ribosylthymine |
        | U | Uridine |
        | I | Inosine |
        | X | Xanthosine |
        | Q | Pseudouridine (\u03A8) |
        | R | Unspecified purine |
        | Y | Unspecified pyrimidine |
        | N | Unspecified/unknown |

        **Residue Frequencies:**
        - For proteins: counts of the 20 standard amino acids (ACDEFGHIKLMNPQRSTVWY)
        - For DNA: counts of A, C, G, T, I, N, X, Q, R, Y
        - For RNA: counts of A, C, G, U, I, N, X, Q, R, Y

        **Length Statistics:**
        - Histogram showing distribution of sequence lengths
        - Summary statistics: min, max, mean, median, 25th/75th percentile lengths

        **Gzipped Files:**
        - Gzipped FASTA files (`.gz`) are supported and automatically decompressed
        """)

    # Input section
    st.subheader("Input Parameters")

    # File uploader (outside form since file uploaders don't work well inside forms)
    uploaded_file = st.file_uploader(
        "Upload FASTA file",
        type=["fasta", "fa", "faa", "fna", "gz"],
        help="Upload a FASTA file (.fasta, .fa, .faa, .fna) or gzipped FASTA (.gz)",
    )

    # Analyze button
    submit = st.button("Analyze Sequences", type="primary")

    # Process submission
    if submit:
        if uploaded_file is None:
            st.error("Please upload a FASTA file to analyze.")
            return

        # Read file content (transparently handles gzipped files)
        try:
            raw_bytes = uploaded_file.read()
            if raw_bytes[:2] == b"\x1f\x8b":
                raw_bytes = gzip.decompress(raw_bytes)
            fasta_input = raw_bytes.decode("utf-8")
        except (UnicodeDecodeError, gzip.BadGzipFile):
            st.error("Could not read file. Please ensure it is a valid (optionally gzipped) FASTA file.")
            return

        if not fasta_input.strip():
            st.error("Please provide FASTA sequences to analyze.")
            return

        # Perform analysis
        with st.spinner("Analyzing sequences..."):
            try:
                results = analyze_fasta(fasta_input, "auto")
            except Exception as e:
                st.error(f"Error parsing FASTA: {str(e)}")
                return

        if not results["success"]:
            st.error(results["error"])
            return

        # Display results
        st.success(f"Successfully analyzed {results['total_sequences']} sequence(s)")

        # Build length DataFrame for metrics and histogram
        length_data = pd.DataFrame({
            "Header": results["length_stats"]["headers"],
            "Length": results["length_stats"]["lengths"],
        })
        q1 = length_data["Length"].quantile(0.25)
        q3 = length_data["Length"].quantile(0.75)

        # Summary metrics
        st.subheader("Summary Statistics")
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Total Sequences", results["total_sequences"])
        with col2:
            st.metric("Total Residues", f"{results['length_stats']['total_residues']:,}")
        with col3:
            st.metric("Avg Length", f"{results['length_stats']['mean']:.1f}")
        with col4:
            st.metric("Median Length", results["length_stats"]["median"])

        col5, col6, col7, col8 = st.columns(4)
        with col5:
            st.metric("Min Length", results["length_stats"]["min"])
        with col6:
            st.metric("25th Percentile", f"{q1:.1f}")
        with col7:
            st.metric("75th Percentile", f"{q3:.1f}")
        with col8:
            st.metric("Max Length", results["length_stats"]["max"])

        # Sequence length histogram
        st.subheader("Sequence Length Distribution")

        if len(length_data) > 1:
            # Calculate number of bins based on max bin width of 100
            length_range = length_data["Length"].max() - length_data["Length"].min()
            nbins = max(1, int(length_range / 100)) if length_range > 0 else 1
            fig_hist = px.histogram(
                length_data,
                x="Length",
                nbins=nbins,
                title="Distribution of Sequence Lengths",
                labels={"Length": "Sequence Length (residues)", "count": "Count"},
            )
            fig_hist.update_layout(
                showlegend=False,
                xaxis_title="Sequence Length (residues)",
                yaxis_title="Number of Sequences",
            )
            st.plotly_chart(fig_hist, use_container_width=True)
        else:
            st.info(f"Single sequence with length: {length_data['Length'].iloc[0]}")

        # Residue frequency analysis
        st.subheader("Residue Frequency Analysis")

        freq_data = results["residue_frequencies"]
        residue_df = pd.DataFrame({
            "Residue": list(freq_data["counts"].keys()),
            "Count": list(freq_data["counts"].values()),
            "Percentage": [f"{p:.2f}%" for p in freq_data["percentages"].values()],
        })
        residue_df = residue_df.sort_values("Count", ascending=False)
        nonzero_residues = residue_df[residue_df["Count"] > 0]

        # Bar chart
        fig_bar = px.bar(
            nonzero_residues,
            x="Residue",
            y="Count",
            title=f"Residue Frequencies ({freq_data['seq_type'].upper()})",
            labels={"Residue": "Residue", "Count": "Count"},
        )
        fig_bar.update_layout(xaxis_tickangle=0)
        st.plotly_chart(fig_bar, use_container_width=True)

        # Frequency table
        col_table, col_chart = st.columns([1, 1])

        with col_table:
            st.markdown("**Residue Counts**")
            st.dataframe(
                nonzero_residues,
                use_container_width=True,
                hide_index=True,
            )

        with col_chart:
            # Pie chart for residue composition
            fig_pie = px.pie(
                nonzero_residues,
                values="Count",
                names="Residue",
                title="Residue Composition",
            )
            st.plotly_chart(fig_pie, use_container_width=True)

        # Download section
        st.subheader("Download Results")

        col_dl1, col_dl2 = st.columns(2)

        # Prepare summary data
        summary_df = pd.DataFrame({
            "Metric": [
                "Total Sequences",
                "Total Residues",
                "Average Length",
                "Min Length",
                "25th Percentile",
                "Median Length",
                "75th Percentile",
                "Max Length",
            ],
            "Value": [
                results["total_sequences"],
                results["length_stats"]["total_residues"],
                f"{results['length_stats']['mean']:.2f}",
                results["length_stats"]["min"],
                f"{q1:.2f}",
                results["length_stats"]["median"],
                f"{q3:.2f}",
                results["length_stats"]["max"],
            ],
        })

        with col_dl1:
            # Combined results as TSV
            tsv_parts = [
                "# FASTA Analysis Summary",
                summary_df.to_csv(sep="\t", index=False),
                "\n# Residue Frequencies",
                residue_df.to_csv(sep="\t", index=False),
            ]
            tsv_data = "\n".join(tsv_parts)

            st.download_button(
                label="Download as TSV",
                data=tsv_data,
                file_name="fasta_analysis_results.tsv",
                mime="text/tab-separated-values",
            )

        with col_dl2:
            # CSV version
            csv_parts = [
                "# FASTA Analysis Summary",
                summary_df.to_csv(index=False),
                "\n# Residue Frequencies",
                residue_df.to_csv(index=False),
            ]
            csv_data = "\n".join(csv_parts)

            st.download_button(
                label="Download as CSV",
                data=csv_data,
                file_name="fasta_analysis_results.csv",
                mime="text/csv",
            )


main()
