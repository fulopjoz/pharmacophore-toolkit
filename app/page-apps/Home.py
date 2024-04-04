import streamlit as st

# Adjust hyperlink colorscheme
links = """<style>
a:link , a:visited{
color: 3081D0;
background-color: transparent;
}

a:hover,  a:active {
color: forestgreen;
background-color: transparent;
}
"""

st.markdown("# Welcome to The Molecular Toolkit")

st.markdown("A Place to Practice Basic Cheminformatic Tools!")

st.markdown('Get started with links in the sidebar.')


# todo modify caption as a fooder (currently only able to use CSS to modify)
st.caption(
    "The Molecular Toolkit is created by [Tony E. Lin](https://github.com/tlint101/molecular-toolkit)"
)
