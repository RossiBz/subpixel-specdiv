# Calculaiting the endmember diversity from remote sensing data

Here we provide the codes used to calculate the endmember diversity from hyperpsectral data in Rossi & Gholizadeh (2023). 

Rossi, C., and Gholizadeh, H. (2023). Uncovering the hidden diversity: Leveraging sub-pixel spectral diversity to estimate plant diversity from space. Remote Sensing of Environment https://doi.org/10.1016/j.rse.2023.113734

The R codes rely on the R package unmix https://github.com/RossiBz/unmix/ developed in the framework of the study.

* **simulated_com_experiment_1.R** provides the code to run experiment one on the simulated communities for the Vertex Component Analysis (VCA) and Purest Pixel Index (PPI) endmember extract techniques. 

* **simulated_com_experiment_1_sisal_mvcnmf.m** is the matlab code to estimate the endmember abundance when  SISAL or MVC-NMF algorithms are used for endmember extraction.
The codes for these algorithms are available from http://www.lx.it.pt/~bioucas/code.html and https://github.com/aicip/MVCNMF

* simulated_communities_field_red100_soil2_3speciesup03_SNR60.csv used in experiment one can be downloaded from https://drive.google.com/file/d/1Z01Hz-8s4KUdhb-YgLu8SyISLNCvZc_K/view?usp=sharing and includes the spectral signature (400-2400 nm), the soil_abundance and the plant diversity (Taxonomic and phylogenetic) of 15'300 simulated communities.
simulated_communities_field_red100_soil2_3speciesup03_SNR60.csv was generated using **generating_simulated_communities.R**. The spectral signatures of the most abundant species used for the simulation can be made available upon request.


* **DESIS_experiment_2.R** provides the code to run experiment two on the DESIS data divided in subimages for the PPI method.
The DESIS data from the Tallgrass prairie preserve can be made available on request. 


## Author

* Christian Rossi christian.rossi1990@gmail.com

## License

Licensed under the GNU General Public License, Version 3.0: https://www.gnu.org/licenses/gpl-3.0.html
