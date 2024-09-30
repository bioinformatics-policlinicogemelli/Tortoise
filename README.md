## TOOL NETWORK

This represents a bioinformatics tool that enables the construction of a bipartite graph from mutational data and the subsequent clustering of the constructed network. 
The tool then allows for the analysis of the identified clusters and their visualization through a dedicated dashboard.

## Installation and usage
### Conda
This script can be run in a Conda environment. To get started, create and activate the environment by running the following commands:
```bash
# Create the environment from the environment.yml file
conda env create -f Tool_Network.yml
```
# Usage
Set correctly the configuration file *config.json*.

```
{
    "Paths":{
        "name_study":"lung",
        "data_mutational":"./lung/data_mutational_filtered.txt",
        "data_clinical_sample":"./lung/data_clinical_sample.txt",
        "data_clinical_patient": "./lung/data_clinical_patient.txt",
        "output_folder":"./lung/output"
    },

    "Mutation":{
        "column_mutation_name":"",
        "gene_of_interest":[],
        "column_gene_name":"Hugo_Symbol",
        "column_hgvsp_short":"HGVSp_Short",
        "column_chromosome":"Chromosome",
        "column_start":"Start_Position",
        "column_end":"End_Position",
        "vaf":true,
        "vaf_score":0.05,
        "vaf_column":""

    },

    "Clinical_data":{
        "column_sample_name":"SAMPLE_ID",
        "column_patient_name":"PATIENT_ID"
    },

    "Enrichment":{
        "adjusted":false,
        "threshold": 0.05,
        "go":true,
        "kegg":true,
        "wiki":true
    }
}

```
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
ì


