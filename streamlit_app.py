import streamlit as st
import glob
import pandas as pd
from targeted_orthofinder import orthofinder
import time

apptitle = 'OrthoFinder'
st.set_page_config(page_title=apptitle, page_icon="🧬")
st.write('# Ortholog Finder')

inputs = st.container()
inputs.write('## Inputs')
inputs.write('### Gene Target')
inputs.write('If you want to upload a new gene here:. Otherwise, it will search for cas9 by default.')
genes1, genes2 = inputs.columns(2)
with genes1:
    upload_fasta = st.file_uploader("Upload fasta files here", 
                                type='fasta', accept_multiple_files=True)
with genes2:
    st.write('Default genes:')
    st.json({
        'default': [f.split('/')[-1] for f in glob.glob('./defaults/inputs/*.fasta')]
    })
inputs.write('### Organisms')
organisms = inputs.text_area("Comma-separated list of organisms you want to search in.",
                        value="Bacillus subtillus, Escherichia coli, Geobacillus stearothermophilus, Streptococcus pyogenes, Streptococcus thermophilus")

run_orthofinder = st.container()
run_orthofinder.write('### Running parameters')
default = run_orthofinder.checkbox("Use default data?", value=False)
clean = run_orthofinder.checkbox("Clean up the output files?")
e_value = run_orthofinder.text_input("E-value cutoff", value="1e-5")
email = run_orthofinder.text_input("Email", value="")
try:
    float(e_value)
except ValueError:
    st.error("E-value must be a number")

def run_orthofinder_func(clean, e_value, email):
    # with st.spinner('Running orthofinder...'):
    ## Save genes
    if len(upload_fasta) > 0:
        # Write the files to "./inputs"
        for f in upload_fasta:
            with open(f"./inputs/{f.name}", "wb") as fout:
                fout.write(f.read())
        st.session_state.genes = [f.name for f in upload_fasta]
    ## Save organisms
    with open(f"./inputs/organisms.txt", "w") as f:
        f.write(organisms)
    orthofinder(clean, e_value, email)
    
run_flow = run_orthofinder.button("Submit orthofinder query", on_click=run_orthofinder_func, args=[clean, float(e_value), email])


results = st.container()
results.write('## Results')
try:
    with open('./inputs/organisms.txt', 'r') as f:
        organisms = f.read()
        results.write('Organisms:  ' + organisms)
except FileNotFoundError:
    pass

# st.write('Genes:')
# st.write(st.session_state.genes)
# st.write('Organisms:')
# st.write(st.session_state.organisms)

# ## Open outputs/summary.csv and read as pandas dataframe
# if st.button("Open summary.csv"):
#     try:
#         summary = pd.read_csv(f"./outputs/summary.csv")
#         st.write(summary)
#     except:
#         print()