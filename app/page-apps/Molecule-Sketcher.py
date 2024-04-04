import streamlit as st
from streamlit_ketcher import st_ketcher
from streamlit_extras.add_vertical_space import add_vertical_space

# # Set page config
# st.set_page_config(layout="wide")
#
# # add logo
# st.sidebar.image('img/py50_logo_only.png', width=150)

# todo Go over background documentation
# https://blog.streamlit.io/introducing-a-chemical-molecule-component-for-your-streamlit-apps/

# Page text
st.markdown("# Molecule Sketcher")
st.markdown("### Draw Molecular Structures")
st.write(
    "Users can translate a SMILES string, draw molecules, or translate a 2D structure into a smiles string."
)
# add_vertical_space(1)

(
    col1,
    col2,
) = st.columns(2)
with col1:
    molecule = st.text_input(label="", placeholder="Input SMILES String")
    smile_code = st_ketcher(molecule, height=600)
    st.markdown(f"Structure SMILES string: {smile_code}")
with col2:
    if smile_code:
        st.header("Calculate Chemical Properties")

        # Include scripts for Rule of 5
