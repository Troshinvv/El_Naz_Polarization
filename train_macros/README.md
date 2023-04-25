# Global Polarization Wagon

Wagon for the analysis of hyperon global polarization. Utilizes Centrality and Event Plane wagons. Consists of two separate parts: 
- MpdGlobalPolarizationMC: for testing the simulation and the performance of polarization extraction method on MCTracks
- MpdGlobalPolarizationRECO: for the polarization analysis of fully reconstructed Lambda-hyperons
## Overview (MpdGlobalPolarizationMC)

Uses MCTracks branch of the simulated dataset to obtain model distributions of global polarization of hyperons ($`P_{y}`$), which is a component of polarization vector **$`P`$**, saved for all hyperons from PHSD model. Obtains angular distributions of daughter particles (e.g. proton for Lambda hyperon), which can be fitted to obtain average value of global polarization.

Input config file (pGlobalPolMC.txt) contains following parameters, which can be modified:

- NITER: (default value = 20) number of angular bins $`\Delta (\phi) = \Psi_{RP} - \phi_{p}`$ bins, the angle difference between RP(or EP) and azimuthal angle of proton in rest frame of $`\Lambda`$.
- NITER_CENT: (default value = 4) number of centralty bins for the analysis. Defined options are 4 (0-10%,10-20%,20-50%,50-100% centrality intervals), 7 (10%-centrality bins up to 70% centrality), 10(10%-centrality bins up to 100% centrality).
- cent_cut: (default value = 70) cut-off value for centrality (rejects events if centrality > cent_cut).
- cent_cut_choice: (default value = 0) choice of centrality cut (0 -no centrality cut, 1 -cut on centrality using cent_cut value).
- particle_choice: (default value = 3122) particle choice for the analysis (pdg number of hyperon). Currently defined for Lambda (pdg = 3122) and anti-Lambda (pdg = -3122) hyperons.

Before the wagon is added into mpdroot, one can start it as:

- root -b -q MpdAnalysisTask2.cxx MpdGlobalPolarizationMC.cxx RunAnalysesMC.C\(\"pGlobalPolMC\"\)

Output file of the wagon contains the following histograms:

- hCentrality: Distribution of centrality for accepted events.
- NCentr: Amount of events in each centrality bin anaylyzed.
- Resolution_EP1_true: $`\cos (\Psi_{EP} - \Psi_{RP}))`$ for calculating true EP resolution.
- Resolution_EP1_exp: $`\cos (\Psi_{EP}^{N} - \Psi_{EP}^{S}))`$ for calculating reconstructed EP resolution.
- Lpolar_y: For each bin of centrality analyzed, the distribution of model $`P_{y}`$ for full hyperons (primary + secondary).
- Lpolar_y_prim: For each bin of centrality analyzed, the distribution of model $`P_{y}`$ for primary hyperons, the mean value of which represents average global polarization.
- PstarRP_hist: $`\Delta (\phi) = \Psi_{RP} - \phi_{p} `$ (w.r.t. RP angle)  distribution of daughter particles (for full hyperons), which can be used to obtain average polarization from fitting.
- PstarRP_hist_prim: $`\Delta (\phi) = \Psi_{RP} - \phi_{p} `$ (w.r.t. RP angle)  distribution of daughter particles (for primary hyperons), which can be used to obtain average polarization from fitting.
- PstarEP_hist: $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$ (w.r.t. EP angle)  distribution of daughter particles (for full hyperons), which can be used to obtain average polarization from fitting.
- PstarEP_hist_prim: $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$ (w.r.t. EP angle)  distribution of daughter particles (for primary hyperons), which can be used to obtain average polarization from fitting.

