# BIOPM helper scripts
Small helper scripts for data transformation and analysis 

List of scripts:
1. convertEdgeRcsv.py: read in miRNA counts <hairpin/mature>_counts.csv file from nf-core/smrnaseq, and converts it to `merged_gene_counts.tsv` as input run rnadeseq. Requires a*
2. convertMirge3Counts.py: read in miRNA counts mir.Counts.csv output file from mirge3 workflow,  and converts it to `merged_gene_counts.tsv` as input to run rnadeseq. Requires a*
3. `filterSalmon.R`: to filter salmon files (quant.sf and quant.genes.sf) for genes, which were called DE using DESeq2 in qbic-pipelines/rnadeseq, but are extremely low expressed and cause wings in Volcano plots (see for example here: https://x.com/bioinformer/status/1570103140741152768).   
             The script has 4 positional arguments:  
             [1] The Salmon folder in the style of <star_salmon>/<QbiCBarcode>   
             [2] The final DE table from the results folder of qbic-pipelines/rnadeseq, usually named "final_DE_gene_list.tsv"  
             [3] The baseMean value to filter for (genes with baseMean < this value will be removed from [1]  
             [4] Metadata file from the results folder of qbic-pipelines/rnadeseq, usually named "metadata.tsv"  
    Typical command to run is  
    **Rscript filterSalmon.R star_salmon final_DE_gene_list.tsv 5 metadata.tsv 2>&1 | tee filterSalmon.log **
   
Notes:
1. a*: an environment running python3 with packages `pandas` installed. 
