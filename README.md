[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11439090.svg)](https://doi.org/10.5281/zenodo.11439090)

# 1. System requirements
GAMA 1.9 or above (code was tested using v1.9.2 on Windows 11). Can be installed on Windows, macOS, Linux, or Docker. Code and installation instructions available at [Installation | GAMA Platform (gama-platform.org)](https://gama-platform.org/wiki/Installation).
# 2. Installation guide
1. Download and install GAMA 1.9 (or above)
2. Download the model file at [MaselkoLab/TMT-Aedes-model: GAMA 1.9 model of Ae. aegypti control with TMT, SIT, and fsRIDL (github.com)](https://github.com/MaselkoLab/TMT-Aedes-model)
3. Create a new project within GAMA, then copy the model file into the `models` folder within the project.
4. Typical install time = 5 minutes
# 3. Demo
For visualisation and demonstration purposes, simulations can be run singly (i.e. `WT`, `TMT`, `fsRIDL`, `SIT` simulations). Batch simulations of all parameter and technology combinations (i.e. the `bulk_batch` simulation) are used to collect data for analysis. The simulations are run in batches of 10 for each combination, with each batch taking ~ 5 minutes to complete on a relatively modern computer. With 3 x 3 parameters and 4 technologies, 144 simulations will be run in total. We performed these simulations on a Ryzen 5 7640U with 32 GB of memory, which took 8 ~ 12 hours to complete.
## Data
The simulations will produce four .csv files: 
1. TMT_PR50.csv
	1. The time to 50% female population reduction (relative to starting population size) for each simulation
2. TMT_PR95.csv
	1. The time to 95% female population reduction (relative to starting population size) for each simulation
3. TMT_step_output.csv
	1. At each time point after the initial release, all simulations will write key parameter values to this file
	2. *The first time point in the first batch may result in multiple header rows being written to the file, which must be removed prior to analysis.*
4. TMT_female_output.csv
	1. After the initial release, adult females will write to this file when they die with key parameters

 We have found it easiest to work with this data by copying it to an SQLite database file. Example data and analysis are available from the following repository:
 
 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11465043.svg)](https://doi.org/10.5281/zenodo.11465043)
# 4. Instructions for use
1. Adjust any experimental variables within the `global` section of the code
2. Adjust any parameter variables within the `bulk_batch` experiment section
3. Run the batch simulation
