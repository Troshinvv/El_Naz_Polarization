# Overview (MpdGlobalPolarizationRECO)

Using PID reconstructs $`\Lambda`$ hyperons via their weak decay channel (by collecting proton and pion pairs, then performing the selection procedure based on the topological parameters). 
Obtains angular distributions of daughter particles for the reconstructed hyperons (e.g. proton for $`\Lambda`$ hyperon), which can be fitted to obtain average value of global polarization. 
The model distributions of global polarization of hyperons ($`P_{y}`$), which is a component of polarization vector **$`P`$**, is obtained for the associated MC tracks and can be used for comparison.

## Structure of input config file (pGlobalPolRECO.txt) 
<details><summary>expand</summary>

**Parameters for event and track selection (the same as in evCentrality wagon)**

- mZvtxCut: (default value = 130) cut on vertex z coordinate

**General parameters for analysis**

- NITER: (default value = 20) number of angular bins $`\Delta (\phi) = \Psi_{RP} - \phi_{p}`$ bins, the angle difference between RP(or EP) and azimuthal angle of proton in rest frame of $`\Lambda`$
- NITER_CENT: (default value = 4) number of centralty bins for the analysis. Defined options are 4 (0-10%,10-20%,20-50%,50-100% centrality intervals), 7 (10%-centrality bins up to 70% centrality), 10(10%-centrality bins up to 100% centrality)
- NITER_ETA: (default value = 6) number of eta bins for the analysis of reconstructed polarization in eta bins 
- NITER_PT: (default value = 5) number of $`p_{T}`$ bins for the analysis of reconstructed polarization in $`p_{T}`$ bins, currently done only for 5 bins
- cent_cut: (default value = 70) cut-off value for centrality (rejects events if centrality > cent_cut)
- cent_cut_choice: (default value = 0) choice of centrality cut (0 - no centrality cut, 1 - cut on centrality using cent_cut value)
- particle_choice: (default value = 3122) particle choice for the analysis (pdg number of hyperon). Currently defined for $`\Lambda`$ (pdg = 3122) and $`\bar\Lambda`$ (pdg = -3122) hyperons
- nMix: (default value = 5) number of events to mix (used for the event mixing method in signal reconstruction of $`\Lambda`$), if nMix=0 event mixing is not used

**PID parameters for analysis**

- MCFile: path to the MC file with geometry (used for refit)
- sigM: (default value = 3.0) sigma for M
- sigE: (default value = 3.0) sigma for E
- energy: (default value = 9.2) energy
- coef: (default value = 1.0) coefficient
- generator: (default value = NSIG) name of generator/PID method 
- tracking: (default value = CFHM) type of tracking

**Topology selection parameters for $`\omega_{2}`$ selection**

- NITER_Selections: (default value = 30) number of values of $`\omega_{2}`$ for selection, number of invariant mass histograms created with corresponding $`\omega_{2}`$ value
- omega_start: (default value = 1.0) starting value for $`\omega_{2}`$ parameter
- omega_step: (default value = 0.1) step value for $`\omega_{2}`$ parameter

**Topology selection parameters for $`\chi`$ selection**

- chi_pi_start: (default value = 7.6) starting value for $`\chi_{\pi}`$ parameter - DCA of pion to the primary vertex, taken in $`\chi^{2}`$-space
- chi_p_start: (default value = 4.2) starting value for $`\chi_{p}`$ parameter - DCA of proton to the primary vertex, taken in $`\chi^{2}`$-space
- chi_V0_start: (default value = 5.6) starting value for $`\chi_{V_{0}}`$ parameter - two-track separation, taken in $`\chi^{2}`$-space
- lambda_path_start: (default value = 1.6) starting value for $`path_{\Lambda}`$ parameter - decay length of $`\Lambda`$ 
- lambda_angle_start: (default value = 0.06) starting value for $`angle_{\Lambda}`$ parameter - pointing angle
- chi_pi_step: (default value = 0.2) step value for $`\chi_{\pi}`$ parameter
- chi_p_step: (default value = 0.2) step value for $`\chi_{p}`$ parameter
- chi_V0_step: (default value = 0.2) step value for $`\chi_{V_{0}}`$ parameter
- lambda_path_step: (default value = 0.2) step value for $`path_{\Lambda}`$ parameter
- lambda_angle_step: (default value = 0.02) step value for $`angle_{\Lambda}`$ parameter

**Analysis parameters**
- selections_values: (e.g. Omega2Selection_values_MB.txt or ChiSelection_values_MB.txt) File with topology selection values
</details>

## Usage of the wagon

The main macro to start the wagon is *RunAnalysesRECO.C*, which has following arguments:
- output: name of the output file without extension (e.g. pGlobalPolRECO)
- analysis_choice: which iteration to perform (either "analysis" and "selection")
- selection_choice: which selection to use (defined for "omega2" and "chi")

For example:
```bash
root -b -q RunAnalysesRECO.C\(\"pGlobalPolRECO\",\"selection\",\"omega2\"\)
```
The macro will use the provided list of input files *list.txt* and the corresponding parameter files for the three used wagons (evCentrality, evPlane and globalPolarization).

First of all, one needs to find the optimal selection parameters for the topology selection of reconstructed hyperons. 
This can be done by choosing "selection" as *analysis_choice* within the main macro. 
Currently the types of selection defined are "omega2" and "chi", which can be specified as *selection_choice* in the main macro. 

In the case of "selection" + "omega2" (one-dimensional scan using parameter $`\omega_{2}`$), the wagon creates invariant mass distributions for each of the corresponding parameters $`\omega_{2}`$ (number of used values is defined in the input file as NITER_Selections, as well as the starting value and the corresponding step value). 
The obtained distributions can then be processed using the macro *mOmega2_Selection_SidebandsFit_MB.C* to find the optimal value of $`\omega_{2}`$ (corresponding to maximum significance) using sidebands method, and save it to the file. 
If in the config file the parameter *nMix > 0*, the obtained distributions can be also processed using event mixing technique, in which case macro *mOmega2_Selection_Mixing_MB.C* should be used to find the optimal value of $`\omega_{2}`$.

In the case of "selection" + "chi" (multi-dimensional scan using topological parameters of the $`\Lambda`$ decay, taken in $`\chi^{2}`$ space), the wagon creates a tree, containing $`\Lambda`$ candidates and all the necessary information about them. 
Once the tree is created, one can perform a multidimensional scan using:

1. macro *Run_SelectHistos_ChiSelection.cc* to collect invariant mass distributions for each set of five parameters scanned: $`\chi_{\pi}`$, $`\chi_{p}`$, $`\chi_{V_{0}}`$, $`path_{\Lambda}`$, $`angle`$. 
To run the macro on the nica-cluster, an example script is provided: *sNICAclusterChiSelection.sge*.
Make sure that the paths are changed to correct ones. 
2. macro *mChi_Selection_SidebandsFit_MB.C* to fit the obtained distributions using sidebands method in order to find the optimal values of selection parameters, and save them to the file. 
macro *mPlotBestChi_Selection_SidebandsFit_MB.C* can be used to plot the results of the fit for the optimal set of selection parameters. 
3. macro *mChi_Selection_Mixing_MB.C* to fit the obtained distributions using event mixing in order to find the optimal values of selection parameters, and save them to the file. 
macro *mPlotBestChi_Selection_Mixing_MB.C* can be used to plot the results of the fit for the optimal set of selection parameters.

The obtained selections file can then be used in the "analysis" (+ "omega2" or + "chi") mode of the wagon, to obtain final distributions.

An example of the script is provided (see *globalPolarization/macros/*) to start the wagon on the nica-cluster: *sNICAclusterGlobalPolRECO.sge*. Make sure that the paths are changed to correct ones. When all the output files are ready, combine them together using *hadd* (do not attempt this for the trees obtained using the mode "selection" + "chi", use the provided script to analyze them):
```bash
hadd -k -f -j 20 Output_globalPolRECO.root Anal_bin*.root
```
## Structure of Output file of the wagon (after "analysis" mode)
<details><summary>For QA</summary>

- hEvents: Number of events (same as in evCentrality wagon)
- hVertex: Vertex distribution (same as in evCentrality wagon)
- hCentrality: Distribution of centrality for accepted events
- hMassL: Full MB invariant mass
- hMassLsig: Full MB invariant mass (signal)
- hMassLbkg: FUll MB invariant mass (background)
- hPIDflag: PID flags
- hLambFlag: Flags for $`\Lambda`$
- hXiFlag: Flags for $`\Xi`$
</details>
<details><summary>Main</summary>

- hNevCentr: Number of events in each centrality bin (necessary for calculating EP resolution).
- hResolution_EP1_true: $`\cos (\Psi_{EP} - \Psi_{RP}))`$ for calculation of true 1st-order EP resolution.
- hResolution_EP1_reco: $`\cos (\Psi_{EP}^{N} - \Psi_{EP}^{S}))`$ for calculation of reconstructed 1st-order EP resolution.
- hPolarY_Full: For each bin of centrality analyzed, the distribution of model $`P_{y}`$ for full hyperons (primary + secondary), obtained from associated MC tracks.
- hPolarY_Prim: For each bin of centrality analyzed, the distribution of model $`P_{y}`$ for primary hyperons, the mean value of which represents average global polarization, obtained from associated MC tracks.
- hDeltaPhiRP_Full: $`\Delta (\phi) = \Psi_{RP} - \phi_{p} `$ (w.r.t. RP angle)  distribution of daughter particles (for full hyperons), which can be used to obtain average polarization from fitting, obtained from associated MC tracks.
- hDeltaPhiRP_Prim: $`\Delta (\phi) = \Psi_{RP} - \phi_{p} `$ (w.r.t. RP angle)  distribution of daughter particles (for primary hyperons), which can be used to obtain average polarization from fitting, obtained from associated MC tracks.
- hDeltaPhiEP_Full: $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$ (w.r.t. EP angle)  distribution of daughter particles (for full hyperons), which can be used to obtain average polarization from fitting, obtained from associated MC tracks.
- hDeltaPhiEP_Prim: $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$ (w.r.t. EP angle)  distribution of daughter particles (for primary hyperons), which can be used to obtain average polarization from fitting, obtained from associated MC tracks.
- hPolvsPt: $`P_{y}`$ vs $`p_{T}`$ of $`\Lambda`$
- hPolvsEta: $`P_{y}`$ vs $`\eta`$ of $`\Lambda`$
- hm0: Invariant mass in bins of $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$, for each centrality interval
- hm0_mixed: Invariant mass in bins of $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$, for each centrality interval (mixed background)
- hm0_ptbin: Invariant mass in bins of $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$, for each $`p_{T}`$ bin for $20 - 50 \%$ centrality interval
- hm0_ptbin_mixed: Invariant mass in bins of $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$, for each $`p_{T}`$ bin for $20 - 50 \%$ centrality interval (mixed background)
- hm0_etabin: Invariant mass in bins of $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$, for each $`\eta`$ bin for $20 - 50 \%$ centrality interval
- hm0_etabin_mixed: Invariant mass in bins of $`\Delta (\phi) = \Psi_{EP} - \phi_{p} `$, for each $`\eta`$ bin for $20 - 50 \%$ centrality interval (mixed background)

</details>

## Utilizing the output

Run the macro *mAnalyzeRECO.C* on the obtained output file from the wagon to get the final results of the MC analysis of global polarization. 