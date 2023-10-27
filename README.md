# agent-based_multi-scale_model
The agent-based multi-scale model for investigating T-cell inheritance, written and developed by Emil Andersson*.

*Computational Biology and Biological Physics group, Center of Environmental and Climate Science, Lund University.

DOI: https://doi.org/10.1101/2023.10.18.562905 

The model simulates a developing T-cell colony starting from a single ETP cell.

The simulations are initialised from run_agents.py where the type of colony (WT or knockdown) and the numbers of colonies can be set. Here, other transcriptional and epigenetic networks could be specified to simulate completely different systems. Helper files for the simulations are: Colony.py, Cell.py, Colourmap.py, Gillespie_methyl_collab.py, OdeSolver.py, timer.py.

Deterministic_analysis.py is an independent program and can be used to run the transcriptional level deterministically. 

After the colonies have been simulated, LCA_classification.py can be run to classify all cells into the LCA-classification system. plot_LCA_statistics.py plots statistics for the LCA categories.

plot_expression.py is used to plot the expression levels resulting from the stochastic simulations for chosen cell lineages.
