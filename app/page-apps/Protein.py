import streamlit as st
from stmol import showmol
import py3Dmol


st.title("Protein Visualizer")
query = st.text_input("Protein Visualizer")

# For quick testing. Remove when finished
query = '1A2C'
# 1A2C
# Structure of thrombin inhibited by AERUGINOSIN298-A from a BLUE-GREEN ALGA

# Todo organize selections
# Documentation: https://napoles-uach-stmol-home-pom051.streamlit.app/Documentation
# todo include options to upload a specific pdb file
# todo include options to visualize protein and give color options

bcolor = st.color_picker('Pick A Color','#ffffff')
style = st.selectbox('style',['cartoon','line','cross','stick','sphere'])
xyzview = py3Dmol.view(query=query)
xyzview.setStyle({style:{'color':'spectrum'}})
xyzview.setBackgroundColor(bcolor)
showmol(xyzview, height = 500,width=800)
