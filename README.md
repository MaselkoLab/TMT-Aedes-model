# 1. System requirements
GAMA 1.9 or above (code was tested using v1.9.2 on Windows 11). Can be installed on Windows, macOS, Linux, or Docker. Code and installation instructions available at [Installation | GAMA Platform (gama-platform.org)](https://gama-platform.org/wiki/Installation).
# 2. Installation guide
1. Download and install GAMA 1.9 (or above)
2. Download the model file at [MaselkoLab/TMT-Aedes-model: GAMA 1.9 model of Ae. aegypti control with TMT, SIT, and fsRIDL (github.com)](https://github.com/MaselkoLab/TMT-Aedes-model)
3. Create a new project within GAMA, then copy the model file into the `models` folder within the project.
4. Typical install time = 5 minutes
# 3. Demo
The simulations are run in batches of 10 for each treatment, with each batch taking 5 ~ 10 minutes to complete on a relatively modern computer. The simulations will produce four .csv files: 
1. TMT_PR50.csv
	1. The time to 50% female population reduction (relative to starting population size) for each simulation
2. TMT_PR95.csv
	1. The time to 95% female population reduction (relative to starting population size) for each simulation
3. TMT_step_output.csv
	1. At each time point after the initial release, all simulations will write key parameter values to this file
	2. *The first time point may result in several header rows being written to the file. Before further analysis, I've found it easiest to remove these extraneous headers manually .*
1. TMT_female_output.csv
	1. After the initial release, adult females will write to this file when they die with key parameters 
# 4. Instructions for use
1. Adjust any experimental variables within the `global` section of the code
2. Run the batch simulations
	1. Takes 5 ~ 10 minutes per batch
	2. By default, simulations are run in batches of 10 for each treatment
	3. Simulations will run for a maximum of the release delay (i.e. days given to the WT population to reach equilibrium prior to introduction of transgenic males) + 120 days
