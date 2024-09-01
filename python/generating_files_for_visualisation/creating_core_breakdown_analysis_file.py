# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:06:46 2019

@author: Jessica Stacey
"""
import argparse
import ast
import re
import pandas as pd
from rdkit import Chem

def input_args():
    parser = argparse.ArgumentParser(description='Arguments for the Generation of the Core Visualisation File')
    parser.add_argument('-i',
                        dest='input',
                        default=r'..\example_input_files\output_node_information.txt',
                        help='Name of the input file')
    parser.add_argument('-o',
                        dest='output',
                        default=r'..\example_input_files\output_core_analysis.txt',
                        type=str,
                        help='Name of the output file containing')
    parser.add_argument('-t',
                        dest='separator',
                        default='\t',
                        help='File separator in the input file')
    args = parser.parse_args()

    return args


def cleaning_function(smarts):
    ### remove numbered wild atoms to just wild atoms
    smarts_cleaned = re.sub(r'\[\d+\*\]', '[*]', smarts)
    ### canonicalise the smarts
    mol = Chem.MolFromSmarts(smarts_cleaned)
    canonical_smarts = Chem.MolToSmiles(mol)
    
    return canonical_smarts

if __name__ == '__main__':
    args = input_args()

    detailed_dataframe = pd.read_csv(args.input, sep=args.separator, header=0)
    
    detailed_dataframe['core_numbered'] = detailed_dataframe['core_numbered'].apply(ast.literal_eval)
    detailed_dataframe['core_info'] = detailed_dataframe['core_info'].apply(ast.literal_eval)
    detailed_dataframe['core_edges'] = detailed_dataframe['core_edges'].apply(ast.literal_eval)
    
    core_nodes_table_info = []
    whole_dataframe = pd.DataFrame(columns=['name', 'group', 'core', 'number_of_molecules', 'functional_groups', 'functional_groups_list', 'functional_groups_breakdown'])
    for core in set(detailed_dataframe['core'].tolist()):
        ## filter the df so that all the core data in
        core_info_data_filtered = detailed_dataframe[detailed_dataframe.core == core]
        
        node_dict=[]
        for index, (idx, row) in enumerate(core_info_data_filtered.iterrows()):
            if index == 0:
                for key, value in row['core_info'].items():
                    value = cleaning_function(value)
                    node_dict.append({'name': row['core_numbered'][key],
                                     'group': ''.join([i for i in row['core_numbered'][key] if not i.isdigit()]),
                                     'core': core,
                                     'number_of_molecules': 1,
                                     'functional_groups': 1,
                                     'functional_groups_list': [value],
                                     'functional_groups_breakdown': [{'SMILES': value,
                                                                     'Image': value,
                                                                     'Occurrences': 1,
                                                                     'core':core,
                                                                     'Activities': [row['pIC50']]}]})
            else:
                for key, value in row['core_info'].items():
                    value = cleaning_function(value)
                    if value in node_dict[key]['functional_groups_list']:
                        fg_idx = node_dict[key]['functional_groups_list'].index(value)
                        node_dict[key]['number_of_molecules'] += 1
                        node_dict[key]['functional_groups_breakdown'][fg_idx]['Occurrences'] += 1
                        node_dict[key]['functional_groups_breakdown'][fg_idx]['Activities'].append(row['pIC50'])
                    else:
                        node_dict[key]['number_of_molecules'] += 1
                        node_dict[key]['functional_groups'] += 1
                        node_dict[key]['functional_groups_list'].append(value)
                        node_dict[key]['functional_groups_breakdown'].append({'SMILES': value,
                                                                              'Image': value,
                                                                              'Occurrences': 1,
                                                                              'core': core,
                                                                              'Activities': [row['pIC50']]})
    
        node_df = pd.DataFrame(node_dict)
        whole_dataframe = whole_dataframe.append(node_df, sort=False)

    whole_dataframe.to_csv(args.output, sep='\t', index=False)
        