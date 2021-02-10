####################
#
# File name :       populate_model.R
# Author :          Moriah Pellowe
# Date created :    September 10, 2020
# Updated :         February 9, 2020
#
# This script will:
#   - calculate the deposition fractions in R
#   - populate the corresponding parameters into a MoBi pkml file
#   - save the simulation
#
# NOTE:
#   - the number of bins (or particle sizes) and the molecule name should match what is in the MoBi pkml file
#
# TIPS:
#   - getSimulationTree() to get proper file path
#   - paste() to break down long lines of code
#
####################

populate_model <- function(pkml_file, molecule_name, particle_sizes_dm, mean_particle_sizes_dm, sd_particle_sizes_dm, logScale = FALSE,
                           breathing_frequency = 15, fraction_inspiratory = 0.5, breath_hold_time_sec = 0,
                           delay_volume_mL = 0, tidal_volume_mL = 1000, bolus_volume_mL = 1000) {
    
    # Load in libraries and scripts
    library(ospsuite)
    #source("deposition_interface_v2.R")
    
    # load in lung model as pkml
    sim <- loadSimulation(pkml_file)
    
    ### PARAMETERS: Define arguments for deposition function ###
    # particle_sizes_dm <- c(2.1e-5)      #c(1.35e-5, 2.1e-5, 2.85e-5)    #c(4e-5, 1e-4, 1.6e-4)
    # mean_particle_sizes_dm <- 2.1e-5    #1e-4
    # sd_particle_sizes_dm <- 0.75e-5     #3e-5
    # logScale <- FALSE
    # molecule_name <- "Salbutamol"
    
    # initialization
    numberOfBins <- length(particle_sizes_dm)
    # read in drug density from simulation
    density_kg_dm3 <- getParameter(paste(molecule_name, "|Density (drug)", sep=""), sim)
    drug_density_kg_m3 <- density_kg_dm3$value*1000
    
    # calculate deposition fractions
    deposition_output <- deposition_interface(particle_sizes_dm, mean_particle_sizes_dm, sd_particle_sizes_dm, drug_density_kg_m3, logScale,
                                              breathing_frequency, fraction_inspiratory, breath_hold_time_sec, delay_volume_mL, tidal_volume_mL, bolus_volume_mL)
    # reduce amount deposited by fraction deposited in lung
    # should it also include ET region? or just lung? depends on assumptions and what F_inh represents
    #deposition_output$distribution_across_gens[1:25,] <- deposition_output$distribution_across_gens[2:25,]*fraction_deposited_in_lung
    
    # set particle radii
    paths <- NULL
    for (bin in 1:numberOfBins) {
        temp <- paste("Applications|Administration Protocol - Lung|Formulation - Particle Dissolution - Polydisperse|Application_1|ParticleBin_",
              toString(bin), "|Particle radius (at t=0)", sep="")
        paths <- c(paths, temp)
    }
    setParameterValuesByPath(paths, particle_sizes_dm, sim)
    
    # set Number_Of_Particles_Factor
    paths <- NULL
    for (bin in 1:numberOfBins) {
        temp <- paste("Applications|Administration Protocol - Lung|Formulation - Particle Dissolution - Polydisperse|Application_1|ParticleBin_",
                      toString(bin), "|Number_Of_Particles_Factor", sep="")
        paths <- c(paths, temp)
    }
    setParameterValuesByPath(paths, deposition_output$number_of_particles_factor, sim)
    
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
                values <- c(values, deposition_output$distribution_across_gens[gen+1,bin]/numberOfSlices)
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
    setParameterValuesByPath(paths, deposition_output$distribution_across_gens[1,], sim)
    
    # set generations without slices
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
            values <- c(values, deposition_output$distribution_across_gens[gen+1,bin])
        }
    }
    setParameterValuesByPath(paths, values, sim)
    par <- getParameter(paths[1], sim)
    par$value

    saveSimulation(sim, paste("populated_", pkml_file, sep=""))
    
    return(deposition_output)
}




#####
#
# Filename:     deposition_interface.R
# Author:       Moriah Pellowe
# Date created: March 9, 2020
#
# v2:           September 30, 2020
#   - if statement added so that beta is in proper format for case of 1 particle bin
#
# This script will calculate the number of particles per L of drug volume, the pdf_particles over the sizes and 24 generations,
# as well as the distribution over the generations for each particle size.
#
#####

deposition_interface <- function(particle_sizes_dm, mean_particle_size_dm, sd_particle_size_dm, drug_density_kg_m3, log_flag=FALSE,
                                 breathing_frequency, fraction_inspiratory, breath_hold_time_sec,
                                 delay_volume_mL, tidal_volume_mL, bolus_volume_mL){
    library(pracma)
    
    #particle_sizes_dm <- c(10^(-9:-5))
    #mean_particle_size <- 1.5e-6
    #sd_particle_size <- 6e-7
    
    ## Parameters
    breath_f_br <- breathing_frequency       # breathing frequency, [1/min]
    breath_fr_in <- fraction_inspiratory     # fraction of breath as inspiratory
    breath_t_b <- breath_hold_time_sec       # breath-hold time [s]
    
    numGens <- 24
    numSizes <- length(particle_sizes_dm)
    
    breath_V_D <- delay_volume_mL                   # Delay volume [mL] 
    breath_V_T <- tidal_volume_mL                   # Tidal volume [mL]  
    breath_V_B <- bolus_volume_mL                   # Bolus volume [mL] 
    
    # functional residual capacity is hard-coded since the Weibel structure is scaled to this value of 3000 mL
    breath_FRC <- 3000                              # Functional residual capacity, FRC [mL]
    
    data_V_lung_tot = 3000
    data_alv_frac_value = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                            0.002, 0.007, 0.02, 0.07, 0.139, 0.282, 0.48)
    
    
    
    ## Initialization
    N <- c(2^(0:(numGens-1)))
    
    stk_in <- matrix(0, nrow=numGens, ncol=numSizes)
    stk_exp <- matrix(0, nrow=numGens, ncol=numSizes)
    IMP_in <- matrix(0, nrow=numGens, ncol=numSizes)
    IMP_exp <- matrix(0, nrow=numGens, ncol=numSizes)
    IMP_b <- matrix(0, nrow=numGens, ncol=numSizes)
    
    eps_in <- matrix(0, nrow=numGens, ncol=numSizes)
    eps_exp <- matrix(0, nrow=numGens, ncol=numSizes)
    SED_in <- matrix(0, nrow=numGens, ncol=numSizes)
    SED_exp <- matrix(0, nrow=numGens, ncol=numSizes)
    SED_b <- matrix(0, nrow=numGens, ncol=numSizes)
    
    re_in <- matrix(0, nrow=numGens, ncol=numSizes)
    re_exp <- matrix(0, nrow=numGens, ncol=numSizes)
    delta_in <- matrix(0, nrow=numGens, ncol=numSizes)
    delta_exp <- matrix(0, nrow=numGens, ncol=numSizes)
    DIF_in <- matrix(0, nrow=numGens, ncol=numSizes)
    DIF_exp <- matrix(0, nrow=numGens, ncol=numSizes)
    DIF_b <- matrix(0, nrow=numGens, ncol=numSizes)
    
    particle_sizes_m <- particle_sizes_dm/10
    particle_radius_dm <- particle_sizes_dm/2
    
    
    
    ## Lung morphometry [dm]
    D <- c(0.1539, 0.1043, 0.071, 0.0479, 0.0385, 0.0299, 0.0239, 0.0197, 0.0159,
           0.0132, 0.0111, 0.0093, 0.0081, 0.007, 0.0063, 0.0056, 0.0051, 0.0046,
           0.0043, 0.004, 0.0038, 0.0037, 0.0035, 0.0035)
    L <- c(1.026, 0.407, 0.1624, 0.065, 0.1086, 0.0915, 0.0769, 0.065, 0.0547,
           0.0462, 0.0393, 0.0333, 0.0282, 0.0231, 0.0197, 0.0171, 0.0141, 0.0121,
           0.01, 0.0085, 0.0071, 0.006, 0.005, 0.0043)
    
    V <- (L*10)*N*pi*(D*10/2)^2         # [cm3]
    
    # Scale alveolar region
    V_add  <- data_V_lung_tot - sum(V)
    V      <- V + c(data_alv_frac_value)*V_add
    
    phi <- rep((pi/4),length(D))
    
    
    
    ## BolusScaling.m
    TLC <- sum(V)
    scale_init <- breath_FRC/TLC
    Vi_init    <- V*scale_init
    
    cum_V <- cumsum(Vi_init)
    cum_Vi_init <- sum(Vi_init)*rep(1, length(Vi_init))- c(0, cum_V[1:(length(cum_V)-1)])
    
    scale_f <- rep(1, length(Vi_init))
    for (i in 1:length(Vi_init)) {
        scale_f[i] <- (sum(Vi_init) + breath_V_D)/cum_Vi_init[i]
    }
    
    scale_t <- rep(1, length(Vi_init))
    alpha <- rep(0, length(Vi_init))
    for (i in 1:length(Vi_init)) {
        # Special treatment for the last generation
        if (i == length(Vi_init)) {
            alpha[i] <- sum(Vi_init[i:length(Vi_init)])
            scale_t[i] <- (sum(Vi_init) + breath_V_D + breath_V_B)/alpha[i]
        } else {
            alpha[i] <- sum(Vi_init[(i+1):length(Vi_init)])
            scale_t[i] <- (sum(Vi_init) + breath_V_D + breath_V_B)/alpha[i]
        }
    } 
    
    checking <- rep(0, length(Vi_init)-1)
    for (i in (1:(length(Vi_init)-1))){
        checking[i] <- (breath_V_T < ((cum_V[i] + breath_V_D + breath_V_B)/(1 - cum_V[i]/sum(Vi_init)))) # This formula accounts for scaling
    }
    # Check if bolus is not washed out in mouth (before trachea)
    if (breath_V_T < (breath_V_D + breath_V_B)){
        i_wash <- length(Vi_init)
    } else {
        i_wash <- which.max(checking)
    }
    scale_t[i_wash:length(scale_t)] <- (sum(Vi_init) + breath_V_T)/sum(Vi_init)    
    
    f_ave <- (scale_f + scale_t)/2
    
    L <- L*f_ave^(1/3)
    D <- D*f_ave^(1/3)
    D_m <- D/10
    L_m <- L/10
    
    V_scaled <- f_ave*Vi_init
    V_cum_scaled <- cumsum(V_scaled)
    
    checking2 <- breath_V_T - breath_V_D < (V_cum_scaled[1:(length(V_cum_scaled)-1)])
    # If checking2 is false for every generation
    if (sum(checking2)==0) {
        imax <- length(Vi_init)
    } else {
        imax <- which.max(checking2)
    }
    V_scaled[imax] <- (breath_V_T - breath_V_D) - V_cum_scaled[imax-1]
    
    if (imax < length(Vi_init)){
        V_scaled[(imax+1):length(V_scaled)] <- Vi_init[(imax+1):length(Vi_init)]
    }
    V <- V_scaled
    V_cum_scaled_updated <- cumsum(V_scaled)
    
    
    
    ## Constants - extrathoracic deposition
    lamda <- 0.066e-6   # 0.066 um, mean free path of air molecule [m]
    
    k <- 1.38064852E-23    # Boltzmann constant [m^2·kg/(s^2·K)]
    T <- 310.65            # 37.5 degree Celsius in Kelvin
    ne <- 1.9224364E-05    # viscosity of air at 37.5 dgr C [kg/m/s] (1.846*10^-5 kg/m/s at 300K)
    
    ## Constants - inertial impaction
    po <- drug_density_kg_m3             # Unit particle density, 1 g/cm3 = 1000 kg/m3
    
    ## Constants - gravitational sedimentation
    g  <- 9.81               # gravitational acceleration [m/s^2]
    
    ## Constants - diffusion
    pa = 1.1372;             # density of air 37.5 degr C [kg/m3]; http://www.gribble.org/cycling/air_density.html
    
    
    
    ## Formulas - extrathoracic deposition
    breath_t_in <- breath_fr_in*(1/breath_f_br)        # min
    breath_t_exp <- (1-breath_fr_in)*(1/breath_f_br)   # min
    
    # Formulas - inertial impaction
    Cd <- 1 + (lamda/particle_sizes_m)*(2.514 + 0.8*exp(-0.55*(particle_sizes_m/lamda)))
    
    # Flow rate (inspiratory and expiratory)
    Q_in <- (breath_V_T/breath_t_in)/1000     # L/min
    Q_exp <- (breath_V_T/breath_t_exp)/1000   # L/min
    
    # Brownian diffusion coefficient
    Dmol <- ((k*T*Cd)/(3*pi*ne*particle_sizes_m)) # (m^2)/s
    Dmol_cm2_s <- Dmol*(100^2)
    
    ## Formulas - inertial impaction
    Q_in_gen_i <- (Q_in*1000/60)/N      # [cm3/s]
    Q_exp_gen_i <- (Q_exp*1000/60)/N   # [cm3/s]
    
    A_gen_i <- pi*(D*10/2)^2  # [cm2] A=r^2*pi
    
    v_in_gen_i <- (Q_in_gen_i/A_gen_i)/100     # [m/s]
    v_exp_gen_i <- (Q_exp_gen_i/A_gen_i)/100   # [m/s]
    
    theta <- L/(4*D)
    
    ## Formulas - Gravitational sedimentation
    vg   <- (po*((particle_sizes_m^2)*g*Cd)/(18*ne))                     # Gravitational settling velocity of a particle
    
    t_i_in    <- V/N/Q_in_gen_i # [s]
    t_i_exp    <- V/N/Q_exp_gen_i # [s]
    
    
    
    ## Extrathoracic deposition
    oral_in  <- 1 - exp(-0.000278*Q_in *(particle_sizes_m*1e6)^2 - 20.4*(Dmol_cm2_s)^0.66*Q_in ^(-0.31))
    oral_exp <- 1 - exp(-0.000278*Q_exp*(particle_sizes_m*1e6)^2 - 20.4*(Dmol_cm2_s)^0.66*Q_exp^(-0.31))
    
    # print(oral_in)
    # print(oral_exp)
    
    
    
    for(j in 1:numSizes){
        for (i in 1:numGens) {
            r_i   <- D_m[i]/2  # radius of tube
            r_i_p <- particle_sizes_m[j]/2   # radius of particle
            
            
            
            ## Inertial impaction
            stk_in[i,j] <- po*(particle_sizes_m[j]^2)*v_in_gen_i[i]*Cd[j]/(9*ne*D_m[i])     # With cunningham, Zhang et al. 1997
            stk_exp[i,j] <- po*(particle_sizes_m[j]^2)*v_exp_gen_i[i]*Cd[j]/(9*ne*D_m[i])   # With cunningham, Zhang et al. 1997
            
            IMP_in[i,j] <- 0.768*theta[i]*stk_in[i,j]                       
            IMP_exp[i,j] <- 0.768*theta[i]*stk_exp[i,j]
            
            
            
            ## Gravitational sedimentation
            eps_in[i,j]   <- 3*vg[j]*t_i_in[i]*cos(phi[i])/(4*D_m[i]) # upright position
            eps_exp[i,j]   <- 3*vg[j]*t_i_exp[i]*cos(phi[i])/(4*D_m[i]) # upright position
            
            # to handle case when eps > 1, since it is an argument of asin(x)
            if (eps_in[i,j] > 1) {
                SED_in[i,j] <- 1    # SEE FEB 10,11 NOTES
            } else {
                SED_in[i,j] <- 
                    2/pi*(2*eps_in[i,j]*(1-eps_in[i,j]^(2/3))^(1/2) - (eps_in[i,j]^(1/3))*(1-eps_in[i,j]^(2/3))^(1/2) + asin(eps_in[i,j]^(1/3)))
            }
            
            if (eps_exp[i,j] > 1) {
                SED_exp[i,j] <- 1   # SEE FEB 10, 11 NOTES
            } else {
                SED_exp[i,j] <- 
                    2/pi*(2*eps_exp[i,j]*(1-eps_exp[i,j]^(2/3))^(1/2) - (eps_exp[i,j]^(1/3))*(1-eps_exp[i,j]^(2/3))^(1/2) + asin(eps_exp[i,j]^(1/3)))                
            }
            
            # Breath hold
            SED_b[i,j] <- 1 - exp( (-4*g*Cd[j]*(r_i_p^2)*breath_t_b*cos(phi[i]))/(9*pi*ne*r_i) )  
            
            
            
            ## Diffusion
            delta_in[i,j] <- Dmol[j]*L_m[i]/(v_in_gen_i[i]*(D_m[i]^2))
            delta_exp[i,j] <- Dmol[j]*L_m[i]/(v_exp_gen_i[i]*(D_m[i]^2))
            
            re_in[i,j]  <- pa*D_m[i]*v_in_gen_i[i]/ne
            re_exp[i,j]  <- pa*D_m[i]*v_exp_gen_i[i]/ne
            
            if (re_in[i,j] > 2000) {   # turbulent flow
                DIF_in[i,j] <- 4*(delta_in[i,j]^0.5) * (1 - 0.444*(delta_in[i,j]^0.5)) # the eq. continues with a (...) in Yu and Diu 1982
            } else {                # laminar flow
                DIF_in[i,j] <- 1-0.819*exp(-14.63*delta_in[i,j])-0.0976*exp(-89.22*delta_in[i,j])-0.0325*exp(-228*delta_in[i,j]) - 0.0509*exp(-125.9*delta_in[i,j]^(2/3))
            }
            if (re_exp[i,j] > 2000) {   # turbulent flow
                DIF_exp[i,j] <- 4*(delta_exp[i,j]^0.5) * (1 - 0.444*(delta_exp[i,j]^0.5)) # the eq. continues with a (...) in Yu and Diu 1982
            } else {                # laminar flow
                DIF_exp[i,j] <- 1-0.819*exp(-14.63*delta_exp[i,j])-0.0976*exp(-89.22*delta_exp[i,j])-0.0325*exp(-228*delta_exp[i,j]) - 0.0509*exp(-125.9*delta_exp[i,j]^(2/3))
            }
            
            # Breath-hold
            DIF_b[i,j] <- 1 - exp(-5.784*k*T*Cd[j]*breath_t_b/(6*pi*ne*r_i_p*(r_i^2)))
        }
    }
    
    # print(IMP_in)
    # print(IMP_exp)
    
    # Remove any complex numbers from asin(x), only necessary when eps > 1
    #SED_in <- Re(SED_in)
    #SED_exp <- Re(SED_exp)
    
    # print(SED_in)
    # print(SED_exp)
    # print(SED_b)
    
    # print(DIF_in)
    # print(DIF_exp)
    # print(DIF_b)
    
    
    
    # Add row at top of each matrix to represent mouth (before trachea)
    IMP_in <- rbind(oral_in, IMP_in)
    IMP_exp <- rbind(oral_exp, IMP_exp)
    IMP_b <- rbind(rep(0,length(oral_in)), IMP_b)
    SED_in <- rbind(rep(0,length(oral_in)), SED_in)
    SED_exp <- rbind(rep(0,length(oral_in)), SED_exp)
    SED_b <- rbind(rep(0,length(oral_in)), SED_b)
    DIF_in <- rbind(rep(0,length(oral_in)), DIF_in)
    DIF_exp <- rbind(rep(0,length(oral_in)), DIF_exp)
    DIF_b <- rbind(rep(0,length(oral_in)), DIF_b)
    V_cum_scaled_updated <- c(0, V_cum_scaled_updated)
    V_scaled <- c(0, V_scaled)
    
    # Compensate for extra row 
    i_wash <- i_wash+1
    imax <- imax+1
    
    
    ## ApplyFractions.m
    numGens_dep <- imax
    
    f <- matrix(0, nrow=(numGens_dep+1), ncol=numSizes)
    DEP_in <- matrix(0, nrow=numGens_dep, ncol=numSizes)
    DEP_exp <- matrix(0, nrow=numGens_dep, ncol=numSizes)
    DEP_b <- matrix(0, nrow=numGens_dep, ncol=numSizes)
    
    P_in <- 1-(1-IMP_in)*(1-SED_in)*(1-DIF_in)
    P_exp <- 1-(1-IMP_exp)*(1-SED_exp)*(1-DIF_exp)
    P_b <- 1-(1-IMP_b)*(1-SED_b)*(1-DIF_b)
    
    for (j in 1:numSizes){
        f[,j] <- cumprod(c(1, (1-P_in[1:numGens_dep,j])))
    }
    
    # BolusScaling.m
    frac_bolus <- rep(0, numGens+1)
    frac_bolus[1:(i_wash-1)] <- rep(1, (i_wash-1))
    
    VF_cum <- rep(0, numGens+1)
    for (i in i_wash:imax) {
        VF_cum[i] = (breath_V_T - breath_V_D - V_cum_scaled_updated[i-1])/breath_V_B
    }
    VF_cum[VF_cum>1] <- 1
    frac_bolus[i_wash:imax] <- VF_cum[i_wash:imax]
    
    for (j in 1:numSizes){
        DEP_in[,j] <- f[1:(dim(f)[1]-1),j]*P_in[1:imax,j]*frac_bolus[1:numGens_dep]
        DEP_b[,j] <- f[1:(dim(f)[1]-1),j]*(1-P_in[1:imax,j])*P_b[1:imax,j]*V_scaled[1:numGens_dep]
    }
    
    # print(DEP_in)
    # print(DEP_b)
    
    # frac_pause -> BolusScaling.m
    frac_pause <- rep(0, length(frac_bolus))
    for (i in 2:imax) {
        if (i < imax) {
            frac_pause[i] <- frac_bolus[i] - frac_bolus[i+1]
        } else if (i == imax) {
            frac_pause[i] <- 1 - sum(frac_pause)   
        }
    }
    
    frac_pause_matrix <- matrix(rep(frac_pause[1:imax], times=numSizes), ncol=numSizes)
    x <- matrix(0, nrow=numGens_dep, ncol=numSizes)
    
    if (numSizes==1) {
        beta <- array(f[3:dim(f)[1],]*frac_pause_matrix[2:dim(frac_pause_matrix)[1],], dim = c(dim(f)[1]-2,1))
    } else {
        beta <- f[3:dim(f)[1],]*frac_pause_matrix[2:dim(frac_pause_matrix)[1],]    
    }
    
    for (i in (numGens_dep-1):1) {
        x[i,]   <- (1-P_exp[i+1,])*x[i+1,] + (1-P_b[i+1,])*beta[i,]     #*(1-P_b(i+1,:))
    }
    
    DEP_exp <- x*P_exp[1:imax,]
    
    # print(DEP_exp)
    
    
    
    # add up the deposition fractions
    total_deposition <- DEP_in + DEP_exp + DEP_b
    # add rows for generations with no deposition
    if (dim(total_deposition)[1] < (numGens+1)) {
        total_deposition <- rbind(total_deposition, matrix(0, nrow=(numGens+1-dim(total_deposition)[1]),
                                                           ncol = dim(total_deposition)[2]))
    }
    
    # Calculate proportion of particles based on radius (in dm)     # NOTE: Boger did this in decimetres
    if (log_flag) {
        proportion_particles <- dlnorm(particle_radius_dm, meanlog=mean_particle_size_dm, sdlog = sd_particle_size_dm)  
        print("Note: a log-normal distribution has not yet been tested.")
    } else {
        proportion_particles <- dnorm(particle_radius_dm, mean=mean_particle_size_dm, sd = sd_particle_size_dm)    
    }
    
    # calculate mass of drug in each bin
    pdf_particles <- matrix(0, nrow=nrow(total_deposition), ncol=ncol(total_deposition))
    # deposited_particles <- matrix(0, nrow=nrow(total_deposition), ncol=ncol(total_deposition))
    for (gen in 1:nrow(total_deposition)) {
        pdf_particles[gen,] <- total_deposition[gen,]*proportion_particles
    }
    
    ## MoBi
    # calculate total volume of all particles, NOTE: volume will be in L
    total_volume_ <- 0
    add_up_r3 <- 0
    for (gen in 1:nrow(pdf_particles)) {
        for (radius in 1:ncol(pdf_particles)) {
            add_up_r3 <- add_up_r3 + pdf_particles[gen,radius]*(particle_sizes_dm[radius]^3)
        }
    }
    total_volume <- (4/3)*pi*add_up_r3
    # calculate number of particles factor 
    number_of_particles_factor <- matrix(0, ncol=ncol(total_deposition))
    number_of_particles_factor <- colSums(pdf_particles)/total_volume
    
    distribution_across_gens <- matrix(0, nrow=nrow(total_deposition), ncol=ncol(total_deposition))
    for (column in 1:ncol(pdf_particles)) {
        distribution_across_gens[,column] <- pdf_particles[,column]/sum(pdf_particles[,column])
    }
    
    
    ## Boger
    # normalizing_factor <- 0
    # for (gen in 1:nrow(pdf_particles)) {
    #     normalizing_factor <- normalizing_factor + trapz(particle_radius_dm, pdf_particles[gen,])
    # }
    # pdf_particles <- (pdf_particles / normalizing_factor) * dose
    # 
    # # calculate number of particles in each bin
    # mass_of_one <- (particle_radius_dm^3)*4*pi/3 * drug_density_g
    # for (gen in 1:nrow(pdf_particles)) {
    #     deposited_particles[gen,] <- pdf_particles[gen,] / mass_of_one    
    # }
    
    output <- list("number_of_particles_factor" = number_of_particles_factor,
                   #"pdf_particles" = pdf_particles,
                   "distribution_across_gens" = distribution_across_gens)
    
    return(output)
}