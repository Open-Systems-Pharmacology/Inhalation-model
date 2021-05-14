# Reset environment
rm(list=ls())

# Load scripts
source("populate_model.R")
source("custom_deposition.R")

molecule_name <- "Ciprofloxacin"
particle_diameters_dm <- 2.6e-5
mean_particle_radius_dm <- 2.6e-5/2
sd_particle_radius_dm <- 1.8e-5/2

########## Empirical equations ##########
# Inhaled
pkml_file <- "inhaled_ciprofloxacin.pkml"
oral_bioavailability <- 0.7
lung_bioavailability <- 0.222

populate_model(pkml_file, molecule_name, particle_diameters_dm, mean_particle_radius_dm, sd_particle_radius_dm,
               oral_bioavailability, lung_bioavailability)

########## Stass et al (2017) deposition ##########
# Inhaled
pkml_file <- "inhaled_ciprofloxacin.pkml"
deposition_fractions <- matrix(c(0.394, rep(0.015875, 24)), nrow=25, ncol=1)
oral_bioavailability <- 0.7
# lung bioavailability is already taken into account in deposition_fractions, so the value is left at 1
lung_bioavailability <- 1

custom_deposition(pkml_file, molecule_name, particle_diameters_dm, mean_particle_radius_dm, sd_particle_radius_dm, 
                  deposition_fractions, oral_bioavailability, lung_bioavailability)