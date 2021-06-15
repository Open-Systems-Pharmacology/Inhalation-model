#####
#
# Filename:     custom_deposition.R
# Author:       Moriah Pellowe
# Date created: March 12, 2021
#
# This script will populate a .pkml file with custom deposition fractions.
# The output was designed to match that of populate_model.R.
# Note that the bioavailabilities will further reduce the deposition fractions given.
#
#####

custom_deposition <- function(pkml_file, molecule_name, particle_diameters_dm, mean_particle_radius_dm, sd_particle_radius_dm,
                              deposition_fractions, 
                              oral_bioavailability = 1, lung_bioavailability = 1, device_bioavailability = 1,
                              logScale = FALSE) {
    
    # load in libraries
    library(ospsuite)
    
    # load in pkml file
    sim <- loadSimulation(pkml_file)
    
    # initialization
    numberOfBins <- length(particle_diameters_dm)
    
    # check deposition_fractions for proper number of fractions
    if ((nrow(deposition_fractions)!=25) | (ncol(deposition_fractions)!=numberOfBins)) {
        print("ERROR: deposition_fractions has the wrong dimensions for the number of bins in the pkml file.")
        return(0)
    }
    
    # read in distribution
    distribution_across_gens <- deposition_fractions
    
    # adjust the deposition fractions
    # note that the oral bioavailability is changed within the MoBi simulation so it is not accounted for here
    distribution_across_gens[1,] <- device_bioavailability*distribution_across_gens[1,]
    
    # absolute bioavailability after inhaled administration with oral charcoal 
    #   = fraction output by device * fraction of drug deposited in lung * lung bioavailability
    # i.e. F_inh,charcoal = F_device * df_lung * F_lung
    distribution_across_gens[2:dim(distribution_across_gens)[1],] <- 
        device_bioavailability*distribution_across_gens[2:dim(distribution_across_gens)[1],]*lung_bioavailability
    
    # set oral bioavailability
    paths <- "Organism|ExtrathoracicRegion|Oral bioavailability - F_oral"
    setParameterValuesByPath(paths, oral_bioavailability, sim)
    
    # set particle radii
    paths <- NULL
    for (bin in 1:numberOfBins) {
        temp <- paste("Applications|Administration Protocol - Lung|Formulation - Particle Dissolution - Polydisperse|Application_1|ParticleBin_",
                      toString(bin), "|Particle radius (at t=0)", sep="")
        paths <- c(paths, temp)
    }
    particle_radius_dm <- particle_diameters_dm/2
    setParameterValuesByPath(paths, particle_radius_dm, sim)
    
    ### calculate number_of_particles_factor - taken from deposition_interface.R
    
    # Calculate proportion of particles based on radius (in dm)     # NOTE: Boger did this in decimetres
    if (logScale) {
        # convert geometric mean and geometric standard deviation to meanlog and sdlog as used in dlnorm
        mu <- log(mean_particle_radius_dm)
        sdev <- sqrt(log(sd_particle_radius_dm)^2)
        proportion_particles <- dlnorm(particle_radius_dm, meanlog = mu, sdlog = sdev)  
        print("Note: Since logScale==TRUE, the mean is interpreted as the geometric mean and the sd as the geometric standard deviation.")
    } else {
        proportion_particles <- dnorm(particle_radius_dm, mean=mean_particle_radius_dm, sd = sd_particle_radius_dm)    
    }
    
    # calculate mass of drug in each bin
    pdf_particles <- matrix(0, nrow=nrow(distribution_across_gens), ncol=ncol(distribution_across_gens))
    # deposited_particles <- matrix(0, nrow=nrow(total_deposition), ncol=ncol(total_deposition))
    for (gen in 1:nrow(distribution_across_gens)) {
        pdf_particles[gen,] <- distribution_across_gens[gen,]*proportion_particles
    }
    
    ## MoBi - see email from Juri Solodenko (AW: Question about dissolution in MoBi) from May 1, 2020
    # calculate total volume of all particles, NOTE: volume will be in L
    total_volume_ <- 0
    add_up_r3 <- 0
    for (gen in 1:nrow(pdf_particles)) {
        for (radius in 1:ncol(pdf_particles)) {
            add_up_r3 <- add_up_r3 + pdf_particles[gen,radius]*(particle_radius_dm[radius]^3)
        }
    }
    total_volume <- (4/3)*pi*add_up_r3
    # calculate number of particles factor 
    number_of_particles_factor <- matrix(0, ncol=ncol(pdf_particles))
    number_of_particles_factor <- colSums(pdf_particles)/total_volume
    
    # set Number_Of_Particles_Factor
    paths <- NULL
    for (bin in 1:numberOfBins) {
        temp <- paste("Applications|Administration Protocol - Lung|Formulation - Particle Dissolution - Polydisperse|Application_1|ParticleBin_",
                      toString(bin), "|Number_Of_Particles_Factor", sep="")
        paths <- c(paths, temp)
    }
    setParameterValuesByPath(paths, number_of_particles_factor, sim)
    
    # set generations with slices
    paths <- NULL
    values <- NULL
    gens_w_slices <- data.frame("Generation" = c(1:2,5:16), "Slices" = c(2,2,3,4,6,7,9,12,14,17,19,21,22,24))
    for (bin in 1:numberOfBins) {
        for (index in 1:nrow(gens_w_slices)) {
            gen <- gens_w_slices$Generation[index]
            numberOfSlices <- gens_w_slices$Slices[index]
            for (slice in 1:numberOfSlices) {
                # add zero to generation or slice if single digit value
                genZero <- ifelse(gen < 10, toString(0), "")
                sliceZero <- ifelse(slice < 10, toString(0), "")
                temp <- paste("Applications|Administration Protocol - Lung|Formulation - Particle Dissolution - Polydisperse|Application_1|ParticleBin_",
                              toString(bin), "|Generation ", genZero, toString(gen), " - Slice ", sliceZero, toString(slice), sep="")
                paths <- c(paths, temp)
                # note that the generation is off by one because 1 corresponds to ET and i+1 corresponds to generation i for row > 1
                values <- c(values, distribution_across_gens[gen+1,bin]/numberOfSlices)
            }
        }
    }
    setParameterValuesByPath(paths, values, sim)
    
    # set ET region
    paths <- NULL
    values <- NULL
    for (bin in 1:numberOfBins){
        temp <- paste("Applications|Administration Protocol - Lung|Formulation - Particle Dissolution - Polydisperse|Application_1|ParticleBin_",
                      toString(bin), "|Number of particles fraction - Extrathoracic", sep="")
        paths <- c(paths, temp)
    }
    setParameterValuesByPath(paths, distribution_across_gens[1,], sim)
    
    # set generation without slices
    paths <- NULL
    values <- NULL
    gens_wo_slices <- c(3:4, 17:24)
    for (bin in 1:numberOfBins) {
        for (gen in gens_wo_slices) {
            # add zero to generation or slice if single digit value
            genZero <- ifelse(gen < 10, toString(0), "")
            temp <- paste("Applications|Administration Protocol - Lung|Formulation - Particle Dissolution - Polydisperse|Application_1|ParticleBin_",
                          toString(bin), "|Number of particles fraction - Generation ", genZero, toString(gen), sep="")
            paths <- c(paths, temp)
            
            # note that the generation is off by one because 1 corresponds to ET and i+1 corresponds to generation i for row > 1
            
            values <- c(values, distribution_across_gens[gen+1,bin])
        }
    }
    setParameterValuesByPath(paths, values, sim)
    par <- getParameter(paths[1], sim)
    par$value
    
    saveSimulation(sim, paste("custom_", pkml_file, sep=""))
    
    deposition_output <- list("number_of_particles_factor" = number_of_particles_factor,
                              "distribution_across_gens" = distribution_across_gens)
    
    return(deposition_output)
}