Lib-INVENT Dataset Tutorials: A Guide to Data Preprocessing
======================================================================================================
Two tutorials are provided here to illustrate the data preparation process for purging a molecular dataset 
and slicing it according to reaction definitions to obtain scaffold and decorations. The tutorials are in the
form of jupyter notebooks with a combination or markup and code cells. Upon appropriate replacement of paths
in the code cells, it is possible to run the notebooks on a local machine to assemble the JSON configuration files 
necessary to run the Lib-INVENT Dataset code.

The first tutorial, data_preparation_demo, illustrates the process of obtaining statistics about the chemical properties
of the dataset and filtering based on these.

The second tutorial, reaction_based_slicing, illustrates the assembly of a JSON file necessary to perform slicing based 
on reaction SMIRKS. The SMIRKS used for slicing of the data for the Lib-INVENT decorator training are provided in the 
data folder. 

### Running the notebooks (command line)
1. Navigate to the Lib-INVENT Dataset directory
`cd <path/to/Lib-INVENT-Dataset>`
2. Install the `lib_invent_data` environment:
`conda env -f environment.yml`
3. Activate the environment:
`conda activate lib_invent_data`
4. Execute `jupyter`:
`jupyter notebook`
5. Copy the link to a browser`