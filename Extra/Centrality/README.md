# Centrality calibration

[[_TOC_]]

## Overview

For this analysis the MC-Glauber based centrality framework is utilized in order to calibrate centrality (developed by P. Parfenov et al, NRNU MEPhI for the MPD collaboration). Here I will make a short overview of the procedure, including a few modifications that were introduced, please see the following link for the detailed description of the framework: [https://github.com/FlowNICA/CentralityFramework](https://github.com/FlowNICA/CentralityFramework)

This framework calibrates centrality of an event dependent on the multiplicity within TPC. For the analysis described here the multiplicity is calculated using MpdGlobalTracks with following criteria:
- 500k events
- $`| \eta |`$ < 0.5
- $`| p_{T} |`$ > 0.15 GeV
- $` N_{hits} `$ > 16
- 10%-centrality bins
- |DCA| < 0.5 cm (optional)

Please note, that these criteria might change, which needs to be taken into account.

## Step-by-step procedure

The instructions given here are a short version of the ones provided within the documentation for the Centrality Framework.

**1.** First of all, you need to install and set up the Glauber Monte Carlo package. 
On the cluster this can be done using "wget":
- wget https://tglaubermc.hepforge.org/downloads/?f=TGlauberMC-3.2.tar.gz
- tar -xf *.tar.gz
- rm *.tar.gz

The main ROOT script is *runglauber_v3.2.C*, you will need to manually add the two colliding systems (Au+Au and Bi+Bi) in this script under the line number 1172:

- else if (TString(name) == "Au3")
    {fN = 197; fR = 6.5541; fA = 0.523; fW = 0; fF = 1; fZ =79;}
- else if (TString(name) == "Bi")
    {fN = 209; fR = 6.75; fA = 0.468; fW = 0; fF = 1; fZ =83;}

As soon as it's ready, you can run the Glauber Monte Carlo, e.g.:
- root -l
- .L runglauber_v3.2.C
- runAndSaveNtuple(5000000, "Au3", "Au3", 29.7)
- .q

Here in the order of appearance: Nev = 5000000 is the number of events for Glauber Monte Carlo (needs to be 10:1 w.r.t. the data), sysA = "Au3" and sysB = "Au3" are the colliding systems in question, signn = 29.7 represents cross-section parameter, which is dependend on the energy of the collision. In this case the example produces the root file for the analysis of Au+Au collisions at $` \sqrt{s_{NN}} = 7.7 \mathrm{GeV} `$ (e.g. *GlauberFile.root* with *TNtuple_nt_Au3_Au3* in it).

The script's running time is rather long, it is recommended to run it in a screen session.

**2.** Secondly, one needs to install the Centrality Framework project:
- source /opt/fairsoft/mpd/new/bin/thisroot.sh 
- git clone https://github.com/FlowNICA/CentralityFramework.git
- cd CentralityFramework/Framework/centrality-master/
- mkdir build/ && cd build/
- cmake ..
- make

In order to start the calibration procedure, several steps must be followed.

**2.1.** Preparing the data: 

For the analysis described here, the script MpdDstReader_modified.C was used (using mpddev version of mpdroot), which is a modified version of the default script (CentralityFramework/Readers/MpdDstReader.C). 
The paths for input and output files need to be set within the script (line 46). 
This script uses the aforementioned selection criteria to calculate the multiplicity within TPC (the main difference w.r.t. the default script is that the absolute value of $` p_{T}`$ is used rather then just $` p_{T}`$, which is important in case of GlobalTracks). The output file (e.g. *MultHisto.root* with *hRefMultSTAR* in it) consists of Reco multuplicity histograms wihout/with the DCA cut, multuplicity using MCTracks, $` p_{T}`$ distributions for Reco and MC for comparison.

The script's running time is rather long, it is recommended to run it in a screen session.

**2.2.** Configuring the framework: 

Prior to starting the calibration procedure, one needs to configure the scripts stored in *CentralityFramework/scripts/template* directory. Change the following lines in *config.txt*:
- line 1: path to the MC-Glauber file 
- line 2: name of the TNtuple from the file (e.g. *TNtuple_nt_Au3_Au3*)
- line 3: path to the multiplicity file 
- line 4: name of the multiplicity histogram in the file (e.g. *hRefMultSTAR*) 

Within this analysis the default parametrization (line 14) was used, which can be changed as well.

Next one needs to configure *parameter.list* script. Here the columns represent the values of $` f_{min}`$, $` f_{max}`$, $` k_{min}`$, $` k_{max}`$, $` Mult_{min}`$, $` Mult_{max}`$, accordingly. Only the last two columns need to be changed, dependend on the range of the multiplicity distribution in question, e.g. if one wishes to have the range [10,270] instead of default [20,360]:
- sed -i "s/20:360/10:270/" parameter.list

The last script to be configure is the main *start.sh*. Please take care of the following lines:
- lines 4,12,13: temporal directory of your choosing
- line 25: full path to the framework core (/path/CentralityFramework/Framework/centrality-master/)
- line 32: short name for this centrality determination run (*COMMIT*)

Note, that you are supposed to start your jobs from this directory (*/path/CentralityFramework/scripts/template*). When everything is ready, start the procedure:
- qsub start.sh

The procedure will run for 2-3 hours, upon completion the results will be stored in the output directory *OUT* within the same template directory as the main script *start.sh*.

**2.3.** Processing the results:
 
 When the output is ready, one has to merge the results:
 - cd /.../CentralityFramework/scripts/template/OUT/COMMIT/jobid/file/root/
 - hadd -k -f -j 20 fit_merged.root fit/*.root

 The next step is to use the *Chi2.C* script with the obtained merged file:
 - cd /.../CentralityFramework/Framework/
 - root -l -b -q Chi2.C'("/.../CentralityFramework/scripts/template/OUT/COMMIT/jobid/file/root/fit_merged.root")'

 This macro will result in a line of the following form:
 - f = 0.34+/-0.122 mu = 0.221908+/-0.197214 k = 41+/-8.243 chi2 = 0.991703+/-0.0813411

 One has to find in the file *parameter.list* the line, corresponding to the obtained parameters *f* and *k*. The number of the line will then correspond to the file in the */.../CentralityFramework/scripts/template/OUT/COMMIT/jobid/file/root/glauber_qa/* directory. For example, the values shown above correspond to the line 65 (f = 0.34 is within the range [0.31, 0.35] and k = 41 is within the range [41,50]). This means that the desired glauber file will be */.../CentralityFramework/scripts/template/OUT/COMMIT/jobid/file/root/glauber_qa/glauber_qa_jobid_65.root*.

 The next macro, *HistoCut.C*, divides the multiplicity distributions (from the data and fitted ones) into centrality classes. Please set the paths to the multiplicity file and to the *glauber_qa* file. In the case that the name of the multiplicity histogram is different, modify it as well in line 5.
 - root -l -b -q HistoCut.C'(10)'

 Resulting file *HistoCutResult.root* contains the multiplicity distributions for each 10%-centrality cut. Finally, the parameters can be mapped from MC-Glauber with the centrality classes using macro *CentralityClasses.C*. You will also need to modify the lines 2 and 3 with the paths to the *glauber_qa* file and the obtained from previous step *HistoCutResult.root*. Also, either comment out or change the hard-coded paths in lines 152, 153, 165, 166, 178, 179, 187, 188, 196, 197, 205, 206, 217, 218, 229, 230, 241, 242, 262, 263.
 - root -l -b -q CentralityClasses.C'(10)'

 This produces the final result of the centrality calibration: **FINAL.root**, which can then be used in the analysis.

## Utilizing the results

There are a few ways to use the results. You can see in the main scripts of the analysis, how the centrality calibration is utilized.
If you wish to print the calibration in an easily readable way, you can use the script *printFinal.C*:
- root -l -b -q printFinal.C'("/.../FINAL.root")'

To save results in a .tex file:

- root -l -b -q printFinal.C'("/.../FINAL.root","./example.tex")'

If you wish to plot the obtained results, you can use the script *Plot_centrality_full.C* provided here. 
Lines 18,19,20 need to be changed to correspond to the full paths of the files, obtained during the steps of the centrality calibration procedure (*glauber_qa*, *FINAL.root*, *HistoCutResult.root*). 
To obtain the file, that is referred to in line 21, use the script *Impact_par_from_tpc_centrality.cc*.
As input for this script, one can either use the files obtained in the first step of the main analysis procedure, or the script *Collect_b_from_dataset.cc*, provided here. 
Change the paths in the script *script_collect_b_from_dataset.sge* for your directories (macro, output and input files), take care that the cuts in line 46 correspond to the ones used for the centrality calibration. Then start it on the cluster:
- qsub script_collect_b_from_dataset.sge 

When the files are ready (e.g. in the directory *output*), one can calculate the required impact parameter distributions:
- make Impact_par_from_tpc_centrality && ./Impact_par_from_tpc_centrality -dca 0 -inname output/*.root -centname /.../FINAL.root -outname B_average_model.root

Here the paths for input files, final centrality file and output file need to provided. parameter dca = 0 (1) corresponds to the choice of calibration without (with) DCA cut.

When all is ready, use the following script to plot the results:
- root -l Plot_centrality_full.C'(10)'
