# SubNetX

The data and scripts contained in this repository allow the user to generate novel 
pathways predicitons for any of the compounds available in the network.

## Installation

The installation can be completed in less than 10 minutes, including installation of 
dependencies and fetching the data from the git repository.

### Requirements

- python 3
- rdkit environment
- networkx

rdkit is required for balance calculation and visualisation.

Please, install rdkit and use rdkit environment for running SubNetX.
First install anaconda: https://docs.anaconda.com/anaconda/install/index.html

Then install rdkit as described here: https://www.rdkit.org/docs/Install.html.

`$ conda create -c conda-forge -n my-rdkit-env rdkit`
`$ conda activate my-rdkit-env`

Since networkx package is not part of the default rdkit environment, install it to the environment as follows when the environment is activated:

`$ conda install networkx`

### Download repository

`$ git clone https://github.com/EPFL-LCSB/SubNetX.git`

If you are installing on macOS, make sure you have Homebrew installed, otherwise you might get "git: 'lfs' is not a git command." error.
Once you installed Homebrew, run

`$ brew install git-lfs`
`$ run git-lfs install`

### Note

Data files are stored using git large file storage (lfs). The make file will install git lfs 
automatically. However, if lfs was not installed previously, the repository has to be 
updated after installation:

`$ git-lfs pull`

This is needed to retrieve the data files from the repository after installation.

The default data used in SubNetX is ARBRE repository data (https://doi.org/10.1016/j.ymben.2022.03.013)

# Test
Ex.  ajmalicine, can be repeated for any other compound
- copy a test project folder from `SubNetX/1_subnetwork_extraction/tutorials/any_mode_of_tutorial/ajmalicine` to `SubNetX/1_subnetwork_extraction/projects/ajmalicine`

Run the code as following

`$ cd SubNetX/1_subnetwork_extraction/code`

`$ python3 Main.py ajmalicine`
 
# Usage

- create a folder with the name of your project in the "projects" directory (e.g. `.../SubNetX/1_subnetwork_extraction/projects/your_compound`)
- copy the parameters file from the `.../1_subnetwork_extraction/defaults` folder to your project folder
- follow the instructions for the parameters.txt adjustment specified in the parameters file (for more details consult the manuscript)

## Network expansion

- Set target compound as LCSB ID of compound (e.g., 1467874237, can be found in data/ARBRE/compounds.csv, cUID column)
- The minimal amount of adjustments for your search is substituting the target ID by your target ID.


You will get the following output:

- output_optimization_input: folder that should be passed to the optimisation stage of the algorithm
- stats: number of compounds and network statistics at the different stages of subnetwork extraction
- auxilary_output: detailed overview of boundaries at each stage and initial pathways
- figures: view of the extracted subnetwork

The .gdf files ready for visualisation in Gephi software are available at arbre/output/{projectname}/visualization/gephifiles

The generation is represented by the color and is labeled in the edges part of the .gdf file in "color VARCHAR" column.

You can install Gephi from https://gephi.org.

# Finishing work with SubNetX

Deactivate your rdkit environment as follows:

`$ conda deactivate`
