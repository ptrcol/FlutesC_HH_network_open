# Code for contact/epidemic analysis of Flutes-C data

The code is divided into several folders, depending on the part of the analysis. All sensitive information (data, intermediate information, files generated starting from the data, etc) is removed via the .gitignore and is available upon reasonable request.

# Pipeline to execute
1) run "network_obj_generation/generate_network_objects_full_data.R" (up to the generation of the three time-point networks line ~ 230).
2) run "epi_fitting/generate_epidemic_info.R".
3) run the remaining part of "network_obj_generation/generate_network_objects_full_data.R" (from line ~ 230, to generate only networks with participant "in common" between time-points).
4) run "network_fitting/network_modelling_full_data_recursive.R" to fit (recursively) candidate network model to the contact data and obtain the best fit model.
5) run "network_fitting/simulate_networks.R" to simulate physical contact networks that are consistent with the data.
6) run "epi_fitting/fitting_main_analysis.R" to get best-fit parameters for the agent based model.
7) run "epi_simulation/results_main_analysis.R" to simulate the impact of physical distancing strategy (agent based model with best-fit parameters coming from 6 and network of pysical contacts coming from 5).


