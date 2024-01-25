# Rhizosphere_assembly_Common_Bean
 Analysis for the rhizoplane and rhizosphere assembly project on common bean

## Github Repository for
# Elucidating the recruitment timing and development of the root microbiome over the common bean lifecycle
## Abby Sulesky-Grieb, A. Fina Bintarti, Keara Grady, Chad Niederhuth, and Ashley Shade 
<i>This work is not published but will soon be available on bioRxiv.</i>


### Data
The raw sequence data for this study are available in the NCBI SRA under bioproject [PRJNA1066866](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1066866/)

Root photos for this study can be found on Figshare at doi:

### To cite this work or code
Coming soon.


### Abstract
Coming soon.


### Contents

Code is split up into four directories: [Rhizo_assembly_seq_processing.Rmd](https://github.com/ShadeLab/Rhizosphere_assembly_Common_Bean/blob/main/Rhizo_assembly_seq_processing.Rmd), [Rhizo_assembly_decontam.Rmd](https://github.com/ShadeLab/Rhizosphere_assembly_Common_Bean/blob/main/Rhizo_assembly_decontam.Rmd), [Rhizo_assembly_analysis.Rmd](https://github.com/ShadeLab/Rhizosphere_assembly_Common_Bean/blob/main/Rhizo_assembly_analysis.Rmd), and [Rhizo_Assembly_health_analysis.Rmd](https://github.com/ShadeLab/Rhizosphere_assembly_Common_Bean/blob/main/Rhizo_Assembly_health_analysis.Rmd). Files necessary to run R code are located in [R_Analysis_Files](https://github.com/ShadeLab/Rhizosphere_assembly_Common_Bean/tree/main/R_Analysis_Files).

#### Sequence processing
Code used for sequence processing including read QC, ASV clustering, taxonomy assignment, and tree building can be found under  [Rhizo_assembly_seq_processing.Rmd](https://github.com/ShadeLab/Rhizosphere_assembly_Common_Bean/blob/main/Rhizo_assembly_seq_processing.Rmd). Scripts were run in QIIME2 using SLURM on the MSU HPCC using slurm batch files with suffix .sb. 

#### Sequence decontamination and rarefaction
Code for sequence decontamination by extraction group can be found in [Rhizo_assembly_decontam.Rmd](https://github.com/ShadeLab/Rhizosphere_assembly_Common_Bean/blob/main/Rhizo_assembly_decontam.Rmd). Output ASV table, taxonomy and metadata files from QIIME2 were used to create a Phyloseq object in R. Decontaminated and rarefied phyloseq object is labelled RA_phyloseq_rarefied.rds. Files are available in R_Analysis_Files folder. 

#### Microbiome Analysis
Formal analysis can be found under [Rhizo_assembly_analysis.Rmd](https://github.com/ShadeLab/Rhizosphere_assembly_Common_Bean/blob/main/Rhizo_assembly_analysis.Rmd). All analysis was run with R and code was run in Rmarkdown. Input files for Procrustes analysis can be found in R_Analysis_Files/Procrustes_analysis_files folder.

#### Plant Phenotypic trait analysis
Above and below-ground plant measurements are located in the R_Analysis_Files folder. Root measurements were performed in ImageJ. Analysis code can be found at [Rhizo_Assembly_health_analysis.Rmd](https://github.com/ShadeLab/Rhizosphere_assembly_Common_Bean/blob/main/Rhizo_Assembly_health_analysis.Rmd).

### Funding
This work was supported by the Michigan State University [Plant Resilience Institute](https://plantresilience.msu.edu). 

### More info
[ShadeLab](http://ashley17061.wixsite.com/shadelab/home)