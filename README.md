## TOOL NETWORK

This represents a bioinformatics tool that enables the construction of a bipartite graph from mutational data and the subsequent clustering of the constructed network. 
The tool then allows for the analysis of the identified clusters and their visualization through a dedicated dashboard.

## Installation and usage
### Conda
This script can be run in a Conda environment. To get started, create and activate the environment by running the following commands:

```bash
# Create the environment from the environment.yml file
conda env create -n tool-network --file=environment.yaml
```

Una volta creato attiviamolo con il comando:

```bash
conda activate tool-network
```

Avviamo l'istanza tool-network tramite il comando:

```bash
python dashboard.py
```

Si consiglia la creazione dello studio in ambiente "windows" o "linux", potendo lavorare in multiprocesso il tempo stimato è di 30 secondi, senza multiprocesso il tempo stimato è di 5 minuti

# Create study

* **Paths**: specify the location of input and output data and the name of the study. <br>
* **Mutation** : specify the setting for mutational analysis. <br> Inside *vaf* it is possibile to specify if you can apply a filter of **V**ariant **A**llele **F**requency (**VAF**) and in *vaf_score* the threshold of this filter. ⚠️Only values greater than or equal to the threshold will be considered. <br>
in *vaf_column* it is possibile to specify the column of VAF value into the mutational file; if it is not specify, the VAF value will be calculated.
* **Clinical_data**: specify the column to identify sample and patient name. <br>

* **Enrichment** : specify the parameters for enrichment analysis. <br> Inside *adjusted* it is possibile to specify if considering the p adjusted value for the enrichment analysis and inside *threshold* the cut-off of significance for p-value or adjusted p-value. <br> Inside *go*, *kegg* and *wiki* it is possibile to specify the database for the enrichment. <br>
## Input files
- **data mutational** (obbligatory)
- **data clinical sample** (optional)
- **data clinical patient** (optional)

It's necessary a data mutational file (maf, csv or txt extension) to start the analysis. 
Clinical informations are required to matching those to the clusters, but are not obbligatory
