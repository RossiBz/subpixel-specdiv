# Calculaiting the endmember diversity 

Here we provide the codes used to calculate the endmember diversity in Rossi & Gholizadeh (2023). 

Rossi, C., and Gholizadeh, H. (2023). Uncovering the hidden diversity: Leveraging sub-pixel spectral diversity to estimate plant diversity from space. Submitted to Remote Sensing of Environment.

The codes rely on the R package unmix https://github.com/RossiBz/unmix/ developed during the study.

simulated_com_experiment_1.R provides the code to run experiment one on the simulated communities for the vca and ppi endmember extract techniques. 
For the endmember extraction techniques, SISAL and MVC-NMF, Matlab was used. The codes are available from http://www.lx.it.pt/~bioucas/code.html and https://github.com/aicip/MVCNMF

simulated_communities_field_red100_soil2_3speciesup03_SNR60.csv includes the spectral signature (400-2400 nm), the soil_abundance and the plant diversity (Taxonomic and phylogenetic) of 15'300 simulated communities to run experiment one.


DESIS_experiment_2.R provides the code to run experiment two on the DESIS data for the ppi method on subimages.
The DESIS data from the Tallgrass prairie preserve can be made available on request. 


## Author

* Christian Rossi christian.rossi1990@gmail.com

## License

Licensed under the GNU General Public License, Version 3.0: https://www.gnu.org/licenses/gpl-3.0.html
