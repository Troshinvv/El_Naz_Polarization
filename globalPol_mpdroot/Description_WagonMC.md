# Overview (MpdGlobalPolarizationMC)

Uses MCTracks branch of the simulated dataset to obtain model distributions of global polarization of hyperons ($`P_{y}`$), which is a component of polarization vector **$`P`$**, saved for all hyperons from PHSD model. Obtains angular distributions of daughter particles (e.g. proton for $`\Lambda`$ hyperon), which can be fitted to obtain average value of global polarization. Currently set up to work for either $`\Lambda`$ (pdg = 3122) or $`\bar\Lambda`$ (pdg = -3122) hyperons.

Note, that the globalPolarization wagon is made to work with the dataset, where hyperon polarization has been included in the model output, and processed within the detector simulation. Currently it will work only with Request 30 from the official MPD productions.

## Structure of input config file (pGlobalPolMC.txt)

- NITER: (default value = 20) number of angular bins $`\Delta (\phi) = \Psi_{RP} - \phi_{p}`$, the angle difference between RP(or EP) and azimuthal angle of proton in rest frame of $`\Lambda`$.
- NITER_CENT: (default value = 4) number of centralty bins for the analysis. Defined options are 4 (0-10%,10-20%,20-50%,50-100% centrality intervals), 7 (10%-centrality bins up to 70% centrality), 10(10%-centrality bins up to 100% centrality).
- cent_cut: (default value = 70) cut-off value for centrality, events with centrality > cent_cut are rejected.
- cent_cut_choice: (default value = 0) choice of centrality cut (0 - no centrality cut, 1 - cut on centrality using cent_cut value).
- particle_choice: (default value = 3122) particle choice for the analysis (pdg number of hyperon). Currently defined for $`\Lambda`$ (pdg = 3122) and $`\bar\Lambda`$ (pdg = -3122) hyperons.

## Usage of the wagon

One can start the wagon using main analysis macro *RunAnalysesMC.C* in *globalPolarization/macros/*.
```bash
root -b -q RunAnalysesMC.C\(\"pGlobalPolMC\"\)
```
The macro will use the provided list of input files *list.txt* and the corresponding parameter files for the three used wagons (evCentrality, evPlane and globalPolarization).

An example of the script is provided (see *globalPolarization/macros/*) to start the wagon on the nica-cluster: *sNICAclusterGlobalPolMC.sge*. Make sure that the paths are changed to correct ones. When all the output files are ready, combine them together using *hadd*:
```bash
hadd -k -f -j 20 Output_globalPolMC.root Anal_bin*.root
```
## Structure of Output file of the wagon

- hCentrality: Distribution of centrality for accepted events.
- hNevCentr: Number of events in each centrality bin (necessary for calculating EP resolution).
- hResolution_EP1_true: $`\cos (\Psi_{EP} - \Psi_{RP}))`$ for calculation of true 1st-order EP resolution.
- hResolution_EP1_reco: $`\cos (\Psi_{EP}^{N} - \Psi_{EP}^{S}))`$ for calculation of reconstructed 1st-order EP resolution.
- hPolarY_Full: For each bin of centrality analyzed, the distribution of model $`P_{y}`$ for full hyperons (primary + secondary).
- hPolarY_Prim: For each bin of centrality analyzed, the distribution of model $`P_{y}`$ for primary hyperons, the mean value of which represents average global polarization.
- hDeltaPhiRP_Full: $`\Delta (\phi) = \Psi_{RP} - \phi_{p} `$ (w.r.t. RP angle)  distribution of daughter particles (for full hyperons), which can be used to obtain average polarization from fitting.
- hDeltaPhiRP_Prim: $`\Delta (\phi) = \Psi_{RP} - \phi_{p} `$ (w.r.t. RP angle)  distribution of daughter particles (for primary hyperons), which can be used to obtain average polarization from fitting.
- hDeltaPhiEP_Full: $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$ (w.r.t. EP angle)  distribution of daughter particles (for full hyperons), which can be used to obtain average polarization from fitting.
- hDeltaPhiEP_Prim: $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$ (w.r.t. EP angle)  distribution of daughter particles (for primary hyperons), which can be used to obtain average polarization from fitting.

## Utilizing the output

Run the macro *mAnalyzeMC.C* on the obtained output file from the wagon to get the final results of the MC analysis of global polarization. 
Make sure that the parameters (number of centrality bins and the particle pdg) correspond to the ones from the *pGlobalPolMC.txt* file used to run the wagon.
```bash
root -l mAnalyzeMC.C'("Output_Wagon.root","Output_Final.root",4,3122)'
```
Provided example macro *mPlotMC.C* can be used to plot the results.
