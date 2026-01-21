"""
Main page for the OpenMS Toolbox App.

This module sets up and displays the Streamlit app for the OpenMS Toolbox.

Usage:
Run this script to launch the OpenMS Toolbox App.

Note:
- If run in local mode, the CAPTCHA control is not applied.
- If not in local mode, CAPTCHA control is applied to verify the user.

Returns:
    None
"""

from pathlib import Path
import streamlit as st

from src.common.common import page_setup, v_space

page_setup(page="main")

st.markdown("# OpenMS Toolbox")
c1, c2 = st.columns(2)
c1.markdown(
    """
## Tools

A collection of web-based tools for mass spectrometry analysis powered by **pyOpenMS**.

- **In Silico Digest**: Generate peptide sequences from protein digestion
- **m/z Calculator**: Calculate mass-to-charge ratios for peptides
- **Isotopic Pattern Calculator**: Generate theoretical isotope distributions
- **Fragment Ion Generation**: Calculate fragment ions for peptide sequences
- **FASTA Analyzer**: Analyze FASTA files for sequence length distributions and residue frequencies
"""
)
v_space(1, c2)
c2.image("assets/openms_transparent_bg_logo.svg", width=300)
if Path("OpenMS-App.zip").exists():
    st.subheader(
        """
Download the latest version for Windows here by clicking the button below.
"""
    )
    with open("OpenMS-App.zip", "rb") as file:
        st.download_button(
            label="Download for Windows",
            data=file,
            file_name="OpenMS-App.zip",
            mime="archive/zip",
            type="primary",
        )
    st.markdown(
        """
Extract the zip file and run the installer (.msi) file to install the app. The app can then be launched using the corresponding desktop icon.
"""
    )