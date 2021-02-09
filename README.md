# Introduction:

This contribution contains the MoBi implementation of a mechanistic physiologically-based pharmacokinetic inhalation model that brings together lung anatomy, particle deposition, particle dissolution, and mucociliary clearance from previous work and connects it to a two-compartment model.

In this model, the lung is separated into the extrathoracic region and 24 lung generations (or sections) that each consist of concentric layers: epithelial lining fluid, epithelium, and subepithelium. The model allows for absorption of particles in the lung as well as oral absorption of particles in the extrathoracic region.

Using MoBi and the ospsuite R package, the user enters parameters describing the physicochemical properties of the molecule, the micro-constants pertaining to the two-compartment model, oral bioavailability, particle size distribution, and dose.

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

## References
[1] [Boger, E., & Fridén, M. (2019). Physiologically based pharmacokinetic/pharmacodynamic modeling accurately predicts the better bronchodilatory effect of inhaled versus oral salbutamol dosage forms. Journal of aerosol medicine and pulmonary drug delivery, 32(1), 1-12.](https://www.liebertpub.com/doi/full/10.1089/jamp.2017.1436)

[2] [Boger, E., & Wigström, O. (2018). A partial differential equation approach to inhalation physiologically based pharmacokinetic modeling. CPT: pharmacometrics & systems pharmacology, 7(10), 638-646.](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12344)

[3] [Weibel, E. R., Cournand, A. F., & Richards, D. W. (1963). Morphometry of the human lung (Vol. 1). Berlin: Springer.](https://link.springer.com/book/10.1007%2F978-3-642-87553-3)

[4] [Willmann, S., Thelen, K., Becker, C., Dressman, J. B., & Lippert, J. (2010). Mechanism-based prediction of particle size-dependent dissolution and absorption: cilostazol pharmacokinetics in dogs. European journal of pharmaceutics and biopharmaceutics, 76(1), 83-94.](https://www.sciencedirect.com/science/article/abs/pii/S0939641110001517)

[5] [Yu, C. P., & Diu, C. K. (1982). A comparative study of aerosol deposition in different lung models. American Industrial Hygiene Association Journal, 43(1), 54-65.](https://oeh.tandfonline.com/doi/abs/10.1080/15298668291410891)
