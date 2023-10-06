# BIOPM helper scripts
Small helper scripts for data transformation and analysis 

List of scripts:
1. convertEdgeRcsv.py: read in miRNA counts <hairpin/mature>_counts.csv file from nf-core/smrnaseq, and converts it to `merged_gene_counts.tsv` as input run rnadeseq. Requires a*
2. convertMirge3Counts.py: read in miRNA counts mir.Counts.csv output file from mirge3 workflow,  and converts it to `merged_gene_counts.tsv` as input to run rnadeseq. Requires a*


Notes:
1. a*: an environment running python3 with packages `pandas` installed. 
