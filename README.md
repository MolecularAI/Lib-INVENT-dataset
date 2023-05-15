**Please note: this repository is no longer being maintained.**

Data preparation and handling for the Lib-INVENT project
========================================================================================================================
This repository holds the code used for preprocessing and data analysis in the Lib-INVENT project 
([Lib-INVENT: Reaction Based Generative Scaffold Decoration for in silico Library Design](https://chemrxiv.org/articles/preprint/Lib-INVENT_Reaction_Based_Generative_Scaffold_Decoration_for_in_silico_Library_Design/14473980 "ChemRxiv preprint")). 


 
The code is organised
in two main folders:
- data_preparation contains scripts necessary for the preparation of the data before the model is trained. The functionality
includes reaction based slicing, dataset profiling and filtering based on user defined criteria.

- results_analysis works with the output of the Lib-INVENT model, scaffold memory csv files, as well as tensorboard logs
associated to the trainings to extract relevant statistics in a csv format.

Note the computation of statistics in the preprocessing stage is relatively time demanding for large datasets despite
the use of PySpark.

Usage
------------------------------------------------------------------------------------------------------------------------
The sole entry point for all runs is the input.py file.

The script needs to be passed a JSON configuration file specifying both the type of run to be implemented (this determines
the action to be taken) and specific parameters relating to the given running mode. The format of the JSON configuration
is as follows:
{
"run_type": ...,
"parameters": {
    ...
    }
}
For example inputs, see the folder example_configurations.

PySpark is used to execute computations on large datasets. If the script uses PySpark, run the code by invoking spark:

>> spark-submit --driver-memory=32g --conf spark.driver.maxResultSize=16g input.py example_configurations/smiles_randomiser.json

If PySpark is not used (e.g. in results analysis), it suffices to call python as:

>> python input.py example_configurations/scaffold_memory_analyser.json

Tutorials
------------------------------------------------------------------------------------------------------------------------
Two tutorials walking through the process of data preparation for Lib-INVENT experiments is provided in the folder
tutorial. We encourage the interested user to walk through the jupyter notebooks contained in this folder to better understand the
general workflow and usability of the project. 

The notebooks give guidance on assembling the JSON input files passed to the input.py script.
Using the JSONs with the arguments presented in the notebooks, it is possible to completely reproduce the data preprocessing 
executed during the preparation of the data used for the training of the Lib-INVENT model.

The filtered [ChEMBL data](https://zenodo.org/record/6627127/files/chembl_train.gz), as well as the [sliced training dataset](https://zenodo.org/record/6627127/files/purged_chembl_sliced.smi.gz?download=1), obtained by this procedure, are provided at Zenodo. A download script is provided as part of the [Lib-INVENT repository](https://github.com/MolecularAI/Lib-INVENT "Lib-INVENT GitHub").
