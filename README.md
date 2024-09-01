# README for Reduced Graph Lead Optimisation Tool

This repository contains the implementation and dataset for our paper ""

## Code 

### Installation
The code can be installed directly GitHub with:

```
$ pip install git+https://github.com/SheffieldChemoinformatics/reduced_graph_visualisation.git
```


### Creating a Conda Environment
Please install a conda environment with all the requirements in visualisation_conda_env.yml  </br>

### Instantiating The Server
 
Start the visualisation server by using the script </br>

```
$ ./start_server.sh
```

or  </br>

```
$ python lead_optimisation_visualisation.py runserver
``` 

The visualisation can then be run on a laptop or computer at the following address: </br>

http://127.0.0.1:5000/

---

## How To Use The Reduced Graph Visualisation

For instructions on how to use the Reduced Graph visualisation can be found in docs/How_to_use_the_RG_core_tool.docs </br>

### To Run Your Own Dataset Through The Reduced Graph Visualisation Tool

To run a new dataset - select 'New Dataset' and select your chosen file, it must be in the format SMILES ID pIC50.  </br>
This may take a while to run depending on the size of the dataset.  </br>

</br>
The workflow can be done independently and then the corresponding output files added to the visualisation  </br>

**Workflow:**
</br>
&nbsp;&nbsp;&nbsp;Step 1: Run python/reduced_graph_code/reduced_graph.py  </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Example code:  </br>

```
$ python reduced_graph.py -i input_smiles.smi -o output_rg.txt
```

</br>
&nbsp;&nbsp;&nbsp;Step 2: Run python/MCS/mcs_similarity_matrix.py  </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Example code:  </br>

```
$ python mcs_similarity_matrix.py -i output_rg.txt -s output_rg.sdf -o rg_mcs.txt
```

</br>
&nbsp;&nbsp;&nbsp;Step 3: Run python/reduced_graph_core_extraction/finding_cores_from_whole_dataset.py  </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Example code:  </br>

```
$ python finding_cores_from_whole_dataset.py -i output_rg.txt -m rg_mcs.txt -o output_rg_core_extraction.txt
```

</br>
&nbsp;&nbsp;&nbsp;Step 4: Run python/generating_files_for_visualisation/generating_file_for_visualisation_coordinates.py  (The -o must be in the format <dataset>_coordinates.txt)</br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Example code:  </br>
  
```
$ python generating_file_for_visualisation_coordinates.py -i output_rg.txt -s output_rg.sdf -o testdataset_coordinates.txt
```

</br>
&nbsp;&nbsp;&nbsp;Step 5: Run python/generating_files_for_visualisation/creating_visualisation_file.py  (The -o must be in the format <dataset>_node_information.txt) </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Example code:  </br>

```
$ python creating_visualisation_file.py-i output_rg.txt -s output_rg.sdf -a testdataset_activities.txt -c output_rg_core_extraction.txt -o testdataset_node_information.txt
```

</br>
&nbsp;&nbsp;&nbsp;Step 6: Run python/generating_files_for_visualisation/creating_core_breakdown_analysis_file.py  (The -o must be in the format <dataset>_core_analysis.txt)  </br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Example code:  </br>

```
$ python creating_core_breakdown_analysis_file.py -i testdataset_node_information.txt -o testdataset_core_analysis.txt
```

</br>

</br>
The output files from Step 4, 5 and 6 need to be added to the datasets folder and add the dataset name too datasets/datasets.json (for this example it would be testdataset) and restart the server to be able to see the dataset within the Reduced Graph Lead Optimisation Tool


# Citing
If you use the Reduced Graph Visualisation in your analysis, please cite our paper.

# Contacts
Please contact me at jessiestacey@msn.com for any questions or comments.
