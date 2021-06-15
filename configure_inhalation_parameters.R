# Reset environment
rm(list=ls())

# Load scripts
source("populate_model.R")
source("custom_deposition.R")

molecule_name <- "Ciprofloxacin"
particle_diameters_dm <- 2.6e-5
geomean_particle_radius_dm <- 2.6e-5/2
gsd_particle_radius_dm <- 1.8e-5/2

########## Empirical equations ##########
# Inhaled
pkml_file <- "inhaled_ciprofloxacin.pkml"
oral_bioavailability <- 0.7
lung_bioavailability <- 0.2509824

populate_model(pkml_file, molecule_name, particle_diameters_dm, geomean_particle_radius_dm, gsd_particle_radius_dm,
			   oral_bioavailability, lung_bioavailability, logScale = TRUE)

########## Stass et al (2017) deposition ##########
# Inhaled
pkml_file <- "inhaled_ciprofloxacin.pkml"
deposition_fractions <- matrix(c(0.394, rep(0.015875, 24)), nrow=25, ncol=1)
oral_bioavailability <- 0.7
# Lung and device bioavailabilities are already taken into account in deposition_fractions, so the values are left at 1

custom_deposition(pkml_file, molecule_name, particle_diameters_dm, geomean_particle_radius_dm, gsd_particle_radius_dm, 
				  deposition_fractions, oral_bioavailability, logScale = TRUE)