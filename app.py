import streamlit as st
from pathlib import Path
import json
# For some reason the windows version only works if this is imported here
import pyopenms

if "settings" not in st.session_state:
        with open("settings.json", "r") as f:
            st.session_state.settings = json.load(f)

if __name__ == '__main__':
    pages = {
        str(st.session_state.settings["app-name"]) : [
            st.Page(Path("content", "quickstart.py"), title="Quickstart", icon="👋"),
        ],
        "Toolbox": [
            st.Page(Path("content", "digest.py"), title="In Silico Digest", icon="✂️"),
            st.Page(Path("content", "peptide_mz_calculator.py"), title="m/z Calculator", icon="⚖️"),
            st.Page(Path("content", "isotope_pattern_generator.py"), title="Isotopic Pattern Calculator", icon="📶"),
            st.Page(Path("content", "fragmentation.py"), title="Fragment Ion Generation", icon="💥"),
            st.Page(Path("content", "fasta_analyzer.py"), title="FASTA Analyzer", icon="🧬"),
        ],
    }

    pg = st.navigation(pages)
    pg.run()
