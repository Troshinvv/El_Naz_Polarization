# MCtests 

[[_TOC_]]

## Overview

This is an overview of possible MC tests, which one can perform on the analyzed dataset.
First part concerns the usage of "Single particle" tests - creation of related histograms and plotting of the results for one type of particle ($`\Lambda`$ or $`\bar\Lambda`$).
Second part concerns the usage for "Full" version of the tests: all possible results, that one can obtain using MC tracks for different particle types, with the corresponding plotting script.

MCTracks of the simulated dataset can be used to obtain distributions of global polarization of hyperons (-$`P_{y}`$), which is a component of polarization vector **$`P`$**, saved for all hyperons from PHSD model. 
Also, the average value of global polarization can be obtained from the fitting of angular distribution of daughter particles (e.g. proton for Lambda hyperon).

## Procedure -- Single Particle

The main script for the MC tests is *MCTest_Single_Particle.cc* with corresponding script for nice-cluster *script_MCTest_Single_Particle.sge*. 
One can start it on the cluster with (from the directory you wish to collect the files in):
- qsub /.../script_MCTest_Single_Particle.sge

Alternatively, you can start it from your working directory (either to check, whether code works, or to avoid using the cluster nodes):
- make MCTest_Single_Particle && ./MCTest_Single_Particle -NITER_CENT 4 -NITER 20 -dca_choice 1 -cut_pt 0.15 -cut_eta 1.0 -dca_cut 1.0 -_N_Hits 16 -angle_choice RP -particle_choice Lambda -CENT_CUT 70. -cent_cut_choice 0 -inname $INFILE_PATH -outname $OUTFILE -centname $CENTFILE

Here the configuration parameters can be modified, according to the purpose of the study:
- NITER_CENT: the amount of centrality bins analyzed (by default, we analysed 4 bins: 0-10%,10-20%,20-50%,50-100%), can be either 4, 7, or 10.
- NITER: the amount of $`\Delta (\phi) = \Psi_{NP} - \phi_{p}`$ bins, the angle difference between RP(EP) and azimuthal angle of proton in rest frame of $`\Lambda`$. Value of 20 was used as default.
- dca_choice: whether the centrality file used was created without (0) or with (1) DCA cut.
- cut_pt, cut_eta, dca_cut, _N_Hits - values for cuts used to obtain Multiplicity distribution. Should correspond to the ones used to obtain centrality file.
- angle_choice: whether we want to calculate the angle difference w.r.t. reaction plane (RP) or event plane (EP).
- particle_choice: analysis for Lambda or ALambda.
- cent_cut_choice: whether we want to cut after a particular value of centrlity (1) or not (0).
- CENT_CUT: the value for centrality cut, after which we disregard the events (if chose to cut).
- inname: path to the input file/files.
- outname: path to the output file/files.
- centname: path to the used centrality file FINAL.root (see [Centrality](../Centrality/) for details about how to obtain it).

Use hadd to combine the files into one, e.g.:
- hadd -k -f -j 100 MCTest_Lambda.root *.root

If you want to collect a particular amount of files (e.g. 500 files for analyzing 1M events):
- hadd -k -f -j 20 MCTest_Lambda_1M.root Anal_bin{0..499}.root

The output contains histograms, corresponding to model distributions of polarization vector, angular distributions of daughter baryons, event plane resolution (if the EP choice was taken). 

To plot the basic set of pictures, one can use the script *Plot_single_particle.C*:
- root -l Plot_single_particle.C'()'

Make sure, that you put the same choices of parameters as for the previous script.
