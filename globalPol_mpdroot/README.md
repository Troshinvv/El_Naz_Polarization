# Global Polarization Wagon

Wagon for the analysis of hyperon global polarization. Utilizes Centrality (evCentrality) and Event Plane (evPlane) wagons. Consists of two separate parts: 
- MpdGlobalPolarizationMC: for testing the simulation and the performance of polarization extraction method on MCTracks
- MpdGlobalPolarizationRECO: for the polarization analysis of fully reconstructed Lambda ($`\Lambda`$) hyperons

The wagon is intended to work with datasets containing polarization (e.g. Request 30).
The reconstruction of hyperons is largely based on the code by A. Zinchenko (Hyperons wagon), with modification to include the necessary information for polarization reconstruction.
## [Overview (MpdGlobalPolarizationMC)](Description_WagonMC.md)

Uses MCTracks branch of the simulated dataset to obtain model distributions of global polarization of hyperons ($`P_{y}`$), which is a component of polarization vector **$`\vec{P}`$**, saved for all hyperons from PHSD model. Obtains angular distributions of daughter particles (e.g. proton for $`\Lambda`$ hyperon), which can be fitted to obtain average value of global polarization. Currently works for either $`\Lambda`$ (pdg = 3122) or $`\bar\Lambda`$ (pdg = -3122) hyperons.

## [Overview (MpdGlobalPolarizationRECO)](Description_WagonRECO.md)

Using PID reconstructs $`\Lambda`$ hyperons via their weak decay channel (collecting proton and pion candidates, then performing the selection procedure based on the topological parameters). Obtains angular distributions of daughter particles for the reconstructed hyperons (e.g. proton for $`\Lambda`$ hyperon), which can be fitted to obtain average value of global polarization. The model distributions of global polarization of hyperons ($`P_{y}`$), which is a component of polarization vector **$`\vec{P}`$**, is obtained for the associated MC tracks and can be used for comparison.

