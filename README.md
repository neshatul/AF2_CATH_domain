# Large-Scale Comparative Analysis of Experimental and AlphaFold2 Predicted Structures

The CATH database is an invaluable resource for protein domain 3D models, providing structural coordinates for domains extracted from proteins available in the RCSB database. In this study, we utilized CATH domain structures and compared them with their corresponding AlphaFold2 (AF2) predicted models. The repository contains scripts for downloading CATH domain structures, AF2 predicted structures, and extracting the respective domains from AF2 models. Additionally, it includes scripts for structural calculations such as root mean square deviation (RMSD), surface area, and volume. Furthermore, it provides tools for analyzing the generated data and plotting the results.

## The analysis is conducted through a series of R notebooks, each designed to perform specific tasks:

Notebook_01_pipeline_02.Rmd: This notebook compiles the scripts used to generate the data for the study. Each section of the notebook performs a specific task, taking input from CATH or AF2 structures. The outputs are saved in the data directory of the repository. Users must download the CATH and AF2 structures before evaluating any code.

Notebook_02_pipeline_data_integration_02.Rmd: This notebook integrates the data generated in Notebook_01 and saves the compiled data in the data directory.

Notebook_03_pipeline_data_analysis_02.Rmd: This notebook generates tables and figures based on the compiled data produced in Notebook_02. The generated figures and tables are also saved in the data directory.

# Contribution
Michael Zimmermann:mtzimmermann@mcw.edu; Neshatul Haque:nehaque@mcw.edu;neshathaq@gmail.com
