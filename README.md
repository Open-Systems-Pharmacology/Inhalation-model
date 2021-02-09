# Introduction:

This contribution contains the MoBi implementation of a mechanistic physiologically-based pharmacokinetic inhalation model that brings together lung anatomy, particle deposition, particle dissolution, and mucociliary clearance from previous work and connects it to a two-compartment model.

In this model, the lung is separated into the extrathoracic region and 24 lung generations (or sections) that each consist of concentric layers: epithelial lining fluid, epithelium, and subepithelium. The model allows for absorption of particles in the lung as well as oral absorption of particles in the extrathoracic region.

Using MoBi and the ospsuite R package, the user enters parameters describing the physicochemical properties of the molecule (e.g. molecular weight, solubility, etc), the micro-constants pertaining to the two-compartment model (e.g. elimination rate, volume of central compartmnet, rate of transfer from central compartment to peripheral compartment, etc.) ,  bioavailability (oral and inhalation), particle size distribution (mean, SD), and dose.

Simulation of the model returns:

- amount or concentration of drug at each compartment and at each generation of the lung over the time of the simulation,
- amount or concentration of drug at each compartment of the two-compartment model as well as the extrathoracic region.

## This repository is organized as follows:

| File name                            | Description                                                                                                                                                                                                                                                                                                                     |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Inhalation_model_user_guide.pdf | Includes:<br>1. Introduction to the model and model details<br>2. A description of the model implementation in Mobi<br>3. Step-by-step tutorial of a basic inhalation simulation                                                                             |
| inhalation_model_two_compt_1_bin.mbp3           | MoBi project template containing the inhalation model connected to a two-compartment model.                                                                                                                                                                                                                                 |
| ParticleBin_1.pkml                   | File containing a single particle bin. This file can be used to add additional particle bins in the case of a polydisperse formulation.                                                                                                                                                                                                                                                       |
| populate_model.R              | File containing R script that is used to calculate deposition fractions. |
