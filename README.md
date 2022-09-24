# sci-RegulatoryClusteringModel Mac M1 version
[![PyPI](https://img.shields.io/pypi/v/scircm)](https://pypi.org/project/scircm/)

If you want to read more about how SiRCle works, please check out our paper: https://www.biorxiv.org/content/10.1101/2022.07.02.498058v1 

Note this is the version which doesn't contain the patient stratification using the VAE since tensorflow doesn't work on the M1 chip.

We will be looking into how we can make this work in the future.

## Install
Optionally create a new conda env.
```
conda create --name scircle python=3.8
conda activate scircle
```

``` 
pip install sircle
```

### Note on Windows
We have tested our code on Windows (10) and Mac (pro) computers, I'm not sure how it would go on a Windows 7 machine so if you have issues post an issue.

### Note on libraries/dependenices
If you find that things don't install please let us know! We've done our best to make it reproducible but different 
environments may mess things up so we're happy to help you debug, just post an issue on the github.

Note we expect python 3.8 so if things don't work first time, check you're running python 3.8 and then try again :) 

## Run
See the examples folder for a proper tutorial with data included that you can run!

*If you have any troubles running the tutorial on a windows machine, please let us know and we can help to sort out any issues. 

#### Quick version
```
from sircle import SciRCM
# FORMAT must be csv :) 
prot_file = f'path to the output from protein differential abundence file'
rna_file = f'path to the output from differential expression analysis file'
meth_file = f'path to the output from methylation DCpG analysis file'

# Note we assume your methylation CpGs map to a single gene, if they don't see the section below.
# logFC_rna = column name in your RNA file that has your RNA logFC (same for the protein and CpG)
# padj_rna = column name in your RNA file that has your padj value (same for protein and CpG)
# NOTE: these need to be unique from one another since we merge the datasets, if they aren't, you need
# to update your csv files.
# Lastly: ensembl_gene_id this is the gene ID column, All must use the same identifier, and this must be
# labelled the same in each file, if it isn't, update your column names before running.

rcm = SciRCM(meth_file, rna_file, prot_file, 
             "logFC_rna", "padj_rna", "CpG_Beta_diff", "padj_meth", "logFC_protein", "padj_protein",
             "ensembl_gene_id", sep=',',
             rna_padj_cutoff=0.05, 
             prot_padj_cutoff=0.05, 
             meth_padj_cutoff=0.05,
             rna_logfc_cutoff=1.0, 
             prot_logfc_cutoff=0.5, 
             meth_diff_cutoff=0.1, 
             output_dir='',
             non_coding_genes=['None'],
             output_filename='RCM_Output.csv',
             bg_type = '(P&M)|(P&R)|(M&R)'
         )
rcm.run()
df = rcm.get_df()
# That DF now has your rcm clustering results, how easy was that :D
```
#### Making your CpGs map to a single gene version
```
from scircm import filter_methylation_data_by_genes
meth_df = pd.read_csv(f'path to the output from methylation DCpG analysis file')
# Note: you need to pass it: 
# 1) the gene ID column, here it is 'ensembl_gene_id'
# 2) the padj column: here it is 'padj_meth'
# 3) the logFC or test statistic column: here it is 'CpG_Beta_diff'
filtered_meth_df = filter_methylation_data_by_genes(meth_df, 'ensembl_gene_id', 'padj_meth', 'CpG_Beta_diff')
```
Now you can run the first version :) 


#### R version
First install Rtools if you haven't done this yet. There are different versions (windows: https://cran.r-project.org/bin/windows/Rtools/, macOS: https://cran.r-project.org/bin/macosx/tools/)

If you don't have conda, you'll need to do the below, first make sure you have reticulate installed.

```
#install.packages('BiocManager')
#BiocManager::install('basilisk')
library(basilisk)
# Set this to be the path to the example data we downloaded
data_dir <- '../data_example/'

protFile <- paste0(data_dir, 'prot_DE_Stage IV_sircle.csv')
rnaFile <- paste0(data_dir, 'rna_DE_Stage IV_sircle_renamed-cols.csv')
methFile <- paste0(data_dir, 'filtered_cpg_DE_Stage IV_sircle.csv')
# Note if you use gene names here you would need to change this to be the column that has the gene ID in it
geneId <- 'ensembl_gene_id'

sircleFileName <- paste0(data_dir, "SircleR-RCM.csv")

# Use basilisk to create an environment we can use
bas_scircm <- BasiliskEnvironment(envname="simple_sircle",
                                   pkgname="sircle",
                                   packages=c("numpy==1.20"),
                                  pip=c("sircle")
)

#logFC_rna = column name in your RNA file that has your RNA logFC (same for the protein and CpG)
#padj_rna = column name in your RNA file that has your padj value (same for protein and CpG)
#NOTE: these need to be unique from one another since we merge the datasets, if they aren't, you need
#to update your csv files.
#Lastly: ensembl_gene_id this is the gene ID column, All must use the same identifier, and this must be
#labelled the same in each file, if it isn't, update your column names before running.
res <- basiliskRun(env=bas_scircm, fun=function(args) {
    rcm <- sircleRCM(rnaFile, methFile, protFile, geneId,  "logFC_rna", "padj_rna", "CpG_Beta_diff", "padj_meth", "logFC_protein", "padj_protein",
                 outputFileName = sircleFileName, 
                 envName="simple_sircle")
    # Do something with pandas
    return(rcm)
})


# Plot the sircle function
sirclePlot(sircleFileName, regLabels="Regulation_Grouping_2") 

# Note you need to have the entrez gene ID added to your csv file
# Run ORA on the groups
sircleORAHuman(sircleFileName, "entrezgene_id", "Regulation_Grouping_2")

```

## Regulatory clustering model 

The general table of how we define regulatory clusters.

| Methylation      | RNAseq    | Proteomics | Regulation driver_1          | Regulation driver_2     | Regulation_Grouping1 | Regulation_Grouping2 | Regulation_Grouping3 |
|------------------|-----------|------------|------------------------------|-------------------------|----------------------|----------------------|----------------------|
| Hypermethylation | DOWN      | DOWN       | Methylation increase (MDS)   | None                    | MDS                  | MDS                  | MDS                  |
| Hypermethylation | UP        | DOWN       | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
| Hypermethylation | UP        | UP         | mRNA increase (TPDE)         | None                    | TPDE                 | TPDE                 | TPDE                 |
| Hypermethylation | DOWN      | UP         | Methylation increase (MDS)   | Protein increase (TMDE) | MDS+TMDE             | TMDE                 | TMDE                 |
| Hypermethylation | No Change | UP         | mRNA increase (TPDE)         | Protein increase (TMDE) | TPDE+TMDE            | TMDE                 | TMDE                 |
| Hypermethylation | No Change | DOWN       | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TMDS                 | TMDS                 |
| Hypermethylation | UP        | No Change  | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
| Hypermethylation | DOWN      | No Change  | Methylation increase (MDS)   | Protein increase (TMDE) | MDS+TMDE             | MDS+TMDE             | TMDE                 |
| Hypermethylation | No Change | No Change  | Methylation increase (ncRNA) | None                    | MDS-ncRNA            | MDS_ncRNA            | MDS_ncRNA            |
| Hypomethylation  | DOWN      | DOWN       | mRNA decrease (TPDS)         | None                    | TPDS                 | TPDS                 | TPDS                 |
| Hypomethylation  | UP        | DOWN       | Methylation decrease (MDE)   | Protein decrease (TMDS) | MDE+TMDS             | TMDS                 | TMDS                 |
| Hypomethylation  | UP        | UP         | Methylation decrease (MDE)   | None                    | MDE                  | MDE                  | MDE                  |
| Hypomethylation  | DOWN      | UP         | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
| Hypomethylation  | No Change | UP         | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TMDE                 | TMDE                 |
| Hypomethylation  | No Change | DOWN       | mRNA decrease (TPDS)         | Protein decrease (TMDS) | TPDS+TMDS            | TMDS                 | TMDS                 |
| Hypomethylation  | UP        | No Change  | Methylation decrease (MDE)   | Protein decrease (TMDS) | MDE+TMDS             | MDE+TMDS             | TMDS                 |
| Hypomethylation  | DOWN      | No Change  | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
| Hypomethylation  | No Change | No Change  | Methylation decrease (ncRNA) | None                    | MDE+ncRNA            | MDE_ncRNA            | MDE_ncRNA            |
| No Change        | DOWN      | UP         | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
| No Change        | UP        | DOWN       | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
| No Change        | DOWN      | DOWN       | mRNA decrease (TPDS)         | None                    | TPDS                 | TPDS                 | TPDS                 |
| No Change        | UP        | UP         | mRNA increase (TPDE)         | None                    | TPDE                 | TPDE                 | TPDE                 |
| No Change        | No Change | UP         | Protein increase (TMDE)      | None                    | TMDE                 | TMDE                 | TMDE                 |
| No Change        | No Change | DOWN       | Protein decrease (TMDS)      | None                    | TMDS                 | TMDS                 | TMDS                 |
| No Change        | UP        | No Change  | mRNA increase (TPDE)         | Protein decrease (TMDS) | TPDE+TMDS            | TPDE+TMDS            | TMDS                 |
| No Change        | DOWN      | No Change  | mRNA decrease (TPDS)         | Protein increase (TMDE) | TPDS+TMDE            | TPDS+TMDE            | TMDE                 |
| No Change        | No Change | No Change  | NoChange                     | NoChange                | NoChange             | NoChange             | NoChange             |

Please post questions and issues related to sci-rcm on the `Issues <https://github.com/ArianeMora/scircm/issues>`_  section of the GitHub repository.

