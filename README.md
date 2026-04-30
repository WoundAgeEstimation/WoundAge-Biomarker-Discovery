# WoundAge-Biomarker-Discovery
The main goal of this project is to characterize gene expression changes in different phases of wound healing and identify biomarkers for wound age estimation (Manuscript in preparation). 

This repostitory hosts the code for all steps of the project, including exploratory gene epxression analysis, identification of differentially expressed genes (DEGs), functional analysis, selection of biomarkers and plotting of the final figures. The identification of DEGs was performed using TAC software (Thermo Fisher Scientific), and is not included in this repository.

## Structure of the repository
The repository is structured as follows:
- `data/`: Folder containing the data used in the analysis (Minus the expression files).
- `code/`: Folder containing all scripts used for the analysis.

- **01_data_exploration_analysis.R** – Code for the exploratory data analysis. Includes PCA analysis, hierarchical clustering and exploration of the covariates.  
- **02_DEG_analysis.R** – Code for the exploration of the DEGs. Run after analysis performed by TAC software.  
- **03_enrichment_analysis.R** – Code for the enrichment analysis of the DEGs.  
- **04_orsum_summary.sh** – Code for the summary of the ORSUM analysis.  
- **05_plot_enrichement_results.R** – Code for the plotting of the enrichment analysis.  
- **06_biomarker_analysis.R** – Code for the identification of the biomarkers.  

Paper figures were generated using these scripts and edited in Inkscape when required. 
