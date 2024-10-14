import streamlit as st
from src.common.common import page_setup
from src.UmetaFlowTOPPWorkflow import Workflow

# The rest of the page can, but does not have to be changed
params = page_setup()

wf = Workflow()

st.title(wf.name)

t = st.tabs(["⚙️ **Configure**", "🚀 **Run**", "📊 **Results**"])

with t[0]:
    wf.show_parameter_section()

with t[1]:
    wf.show_execution_section()
    
with t[2]:
    wf.show_results_section()

