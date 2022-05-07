from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq, translate
import Bio.SeqIO as SeqIO
from Bio import Entrez

import argparse, csv, glob, os, re, sys
from itertools import product
from multiprocessing import Pool

# Get pairs of gene targets and organism targets
def get_blastp_queries():
    # Finds all files with .fasta extension in ./inputs
    files = glob.glob("./inputs/*.fasta")
    # Loads comma-separated organism names from ./inputs/organisms.txt
    organisms = [line.rstrip('\n').split(',') for line in open("./inputs/organisms.txt")][0]
    organisms = [o.strip() + "[orgn]" for o in organisms if o != ""]
    return list(product(files, organisms))

# blastp genes in the organism targets
def blastp(query):
    file, organism = query
    gene = file.split('/')[-1].split('.')[0]
    record = SeqIO.read(file, "fasta")
    # Translates the sequence
    translated = translate(record.seq)
    print(f"translated, sending {(gene, organism)} to blastp")
    try:
        results_handler = NCBIWWW.qblast("blastp","nr", translated,
                                            entrez_query=organism)
        with open(f"./outputs/blastp/blastp_results_{gene}_{organism}.xml", "w") as f:
            f.write(results_handler.read())
        print(f"blastp received, {(gene, organism)}")
    except:
        print(f"Error getting blastp results for {(gene, organism)}")

def get_tblastn_queries(e_value=1e-10):
    files = glob.glob("./outputs/blastp/*.xml")
    good_queries = []
    for f in files:
        with open(f, "r") as result_handle:
            blast_record = NCBIXML.read(result_handle)
            gene, organism = f.split('/')[-1].split('_')[2], f.split('/')[-1].split('_')[3].split('.')[0]
            found_alignment=False
            if len(blast_record.alignments) > 0:
                wrong_organisms = set()
                for alignment in blast_record.alignments:
                    # iterate through the alignments and find the one with the right organism name
                    regex = re.compile(r'\[(.*?)\]')
                    found_organism_name = regex.search(alignment.title).group(1) + "[orgn]"
                    same_organism = found_organism_name.strip() == organism.strip()
                    if not same_organism:
                        # If we already found a wrong organism, don't ask again
                        if found_organism_name in wrong_organisms:
                            continue

                        # If we haven't seen this organism yet, ask the user if it's the right one
                        ask_same_organism = input(f"Is the organism {found_organism_name} the same as {organism} (the query)? (y/n)")
                        if not (ask_same_organism.lower() == 'y'):
                            print(f"Okay, we won't use {found_organism_name}.")
                            wrong_organisms.add(found_organism_name)
                            continue
                    # If we get here, we found the right organism
                    print(f"Found the right organism: {found_organism_name}")
                    print(f"The alignment is: {alignment.title}")
                    print(f"The evalue is: {alignment.hsps[0].expect}")
                    print(f"The hit_id is: {alignment.hit_id}")
                    print('\n')
                    if alignment.hsps[0].expect < e_value:
                        good_queries.append((gene, organism, found_organism_name, alignment.accession))
                        found_alignment = True
                        break
            if not found_alignment:
                good_queries.append((gene, organism, None, None))
    return good_queries

# tblastn the found proteins in the organisms to get the original sequence
def tblastn(query):
    gene, organism, found_organism_name, accession = query
    print(f"found_organism_name: {found_organism_name}")
    if found_organism_name:
        print(f"sending {(gene, organism)} to tblastn")
        try:
            results_handler = NCBIWWW.qblast("tblastn","nt", accession,
                                                filter=None, entrez_query=found_organism_name)
            with open(f"./outputs/tblastn/tblastn_results_{gene}_{organism}.xml", "w") as f:
                f.write(results_handler.read())
            print(f"tblastn received, {(gene, organism)}")
        except:
            print(f"Error getting tblastn results for {(gene, organism)}")

# Get the original sequence from the tblastn results
def get_original_sequence(query):
    gene, organism, found_organism, accession, email = query
    cds_start, cds_end = None, None
    if found_organism:
        try:
            result_handle = open(f"./outputs/tblastn/tblastn_results_cas9_{organism}.xml", "r")
            blast_record = NCBIXML.read(result_handle)
            if len(blast_record.alignments) > 0:
                start = blast_record.alignments[0].hsps[0].sbjct_start
                padded_start = max(0, blast_record.alignments[0].hsps[0].sbjct_start - 200)
                cds_start = start - padded_start
                end = blast_record.alignments[0].hsps[0].sbjct_end
                padded_end = min(blast_record.alignments[0].hsps[0].sbjct_end + 200, blast_record.alignments[0].length)
                cds_end = padded_end - end
                gene_accession = blast_record.alignments[0].accession

                # Get the sequence from the accession
                handle = Entrez.efetch(email=email,db="nucleotide", id=gene_accession, rettype="gb", retmode="text", seq_start=padded_start, seq_stop=padded_end)
                with open(f"./outputs/genbank/genbank_results_{gene}_{organism}.gb", "w") as f:
                    f.write(handle.read())
                print(f"genbank received, {(gene, organism)}")
        except:
            print(f"No tblastn file found for {(gene, organism)}")
    return (gene, organism, found_organism, accession, cds_start, cds_end)

def get_csv_entries(query):
    gene, organism, found_organism, accession, cds_start, cds_end = query
    pre200, cds, post200 = None, None, None
    if found_organism:
        try:
            with open(f"./outputs/genbank/genbank_results_{gene}_{organism}.gb", "r") as f:
                record = SeqIO.read(f, "genbank")
                seq = record.seq
                pre200, cds, post200 = seq[:cds_start], seq[cds_start:-1*cds_end], seq[-1*cds_end:]
        except:
            print(f"No genbank file for {(gene, organism)}")
    return (gene, organism, found_organism, accession, cds_start, cds_end, pre200, cds, post200)

def validate_folder_structure():
    # Check that the inputs folder exists and has .fasta files and one organisms.txt file
    if not os.path.isdir("./inputs"):
        print("The inputs folder does not exist. Please create it and put the .fasta files and organisms file in it.")
        sys.exit(1)
    if not os.path.isfile("./inputs/organisms.txt"):
        print("The organisms.txt file does not exist. Please create it and put the organisms in it, separated by commas.")
        sys.exit(1)

    # Make sure .outputs, .outputs/blastp, .outputs/tblastn, and .outputs/genbank exist
    if not os.path.exists("./outputs"):
        os.mkdir("./outputs")
    if not os.path.exists("./outputs/blastp"):
        os.mkdir("./outputs/blastp")
    if not os.path.exists("./outputs/tblastn"):
        os.mkdir("./outputs/tblastn")
    if not os.path.exists("./outputs/genbank"):
        os.mkdir("./outputs/genbank")

if __name__ == "__main__":
    # Define flags
    parser = argparse.ArgumentParser(description='Run targeted orthofinder')
    parser.add_argument('-c', '--clean', action='store_true', help='Clean up files')
    parser.add_argument('-e', '--e_value', type=float, default=0.01, help='E-value threshold')
    parser.add_argument('-m', '--email', type=str, default="", help='Email address')
    args = parser.parse_args()
    e_value = args.e_value
    validate_folder_structure()

    # Start pipeline
    queries = get_blastp_queries()
    with Pool() as p:
        p.map(blastp, queries)

    queries = get_tblastn_queries(e_value)
    with Pool() as p:
        p.map(tblastn, queries)

    # Add args.email to the queries
    queries = [(gene, organism, found_organism_name, accession, args.email) for gene, organism, found_organism_name, accession in queries]
    with Pool() as p:
        queries = p.map(get_original_sequence, queries)

    with Pool() as p:
        csv_entries = p.map(get_csv_entries, queries)
        # Save array csv_entries to outputs/summary.csv
        with open(f"./outputs/summary.csv", "w") as f:
            writer = csv.writer(f)
            # Write each fasta file in inputs to the csv file
            genes = glob.glob("./inputs/*.fasta")
            for gene in genes:
                seq = SeqIO.read(gene, "fasta").seq
                gene = gene.split('/')[-1].split('.')[0]
                writer.writerow([gene + ":", seq])
            # Write headers gene, organism, found_organism, accession, cds_start, cds_end, pre200, cds, post200
            writer.writerow(["gene", "organism", "found_organism", "accession", "cds_start", "cds_end", "pre200", "cds", "post200"])
            writer.writerows(csv_entries)
        print("Done!")
        print("Summary saved to outputs/summary.csv")

    # Clean up files
    if args.clean:
        print("Cleaning up files")
        os.system("rm -rf ./outputs/blastp")
        os.system("rm -rf ./outputs/tblastn")
