# Orthofinder

Orthofinder is a command line application which finds orthologs of select genes in select organisms. To use it, the user must identify their gene of interest, download the FASTA file from NCBI, and then write the list of organisms they want to query.

### Requirements

Orthofinder is written in Python and uses the BioPython package. It requires access to the internet.

To setup the application, it's recommended to use a virtual environment (venv or conda).
For example:

    conda create -n orthofinder python=3.8
    conda activate orthofinder
    conda install biopython

### Setup

Orthofinder takes two types of inputs: a FASTA file of a gene target, and a list of organisms that the user wants to search for orthologs in. 

    orthofinder
    ├── inputs
    │   ├── cas9.fasta
    │   └── organisms.txt
    └── outputs

Organisms.txt is a comma-separated list of organism names, such as the following:

    Bacillus subtillus, Escherichia coli, Geobacillus stearothermophilus, Streptococcus pyogenes, Streptococcus thermophilus

##### Flags
| Flag      | Alternate flag | Description      | Example | Required? |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| -c      | --clean       | Whether or not to clean up the folders once done.      | -c       | No, but recommended|
| -m      | --email       | The user's email. One of the API requests (getting the sequence of the gene target once found) is rate-limited, and they will email the user if there are any issues.      | -m mnemeth6@berkeley.edu       | No, but recommended|
| -e      | --e_value       | E-value for the BLASTp search.      | -e 1e-3       | No. The default is 0.01|


**Examples**

    python targeted_orthofinder.py -c -m mnemeth6@berkeley.edu -e 1e-3

    python targeted_orthofinder.py -m mnemeth6@berkeley.edu

### Notes

This runs quite slowly. NCBI throttles requests from API access, and recommends the user setting up their own cloud BLAST server for very large requests.