##############################################################################
# Finds the maximum common substructure
# Finds the disconnected and connected for the chemical and reduced graph
#
#
# Jess Stacey
##############################################################################

import argparse
import time
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFMCS

def input_args():
    ''' Function: input_args
    This function establishes the arguments for the python code
    Input: None
    Output: args - argparser
    '''
    parser = argparse.ArgumentParser(description='Arguments for the Reduced Graph Code')
    parser.add_argument('-i',
                        dest='input',
                        default=r'..\example_input_files\output_rg_file.txt',
                        type=str,
                        help='Name of the input file')
    parser.add_argument('-s',
                        dest='output_smiles',
                        default=r'..\example_input_files\output_rg_file.sdf',
                        type=str,
                        help='Name of the output of the smiles mcs file')
    parser.add_argument('-r',
                        dest='output_rg',
                        default=r'..\example_input_files\output_rg_sim_matrix.txt',
                        type=str,
                        help='Name of the output of the smiles rg file')
    parser.add_argument('-p',
                        dest='representationtype',
                        default='rg',
                        type=str,
                        help='The type of representation used to cluster: both smiles rg')
    parser.add_argument('-t',
                        dest='separator',
                        default='\t',
                        help='File separator in the input file')
    parser.add_argument('-v',
                        dest='verbose',
                        action='store_true',
                        help='Switches on verbose mode')

    args = parser.parse_args()

    return args

def read_in_file():
    ''' Function: read_in_file
    This function reads in the file and generates list and dictionary to allow
    the rest of the document to work on the data
    Input: None
    Output: smiles_list - [] of str of SMILES
            reduced_graph_list - [] of str of RG
            smiles_dict - df of key: SMILES and value: RG
    '''
    input_dataframe = pd.read_csv(args.input, sep=args.separator, header=0)
    smiles_list = input_dataframe['SMILES'].tolist()
    reduced_graph_list = input_dataframe['RG'].tolist()

    return(smiles_list, reduced_graph_list)

def tanimoto_mcs(a,
                 b,
                 ab):
    ''' Function: tanimoto_mcs
    This function calculates the tanimoto similarity
    Input: a - int from atoms in molecule 1
           b - int from atoms in molecule 1
           ab - int from atoms in both molecule 1 and molecule 2
    Output: tcms - float of the tanimoto similarity
    '''
    tcmcs = (ab)/(a + b - ab)

    return tcmcs

def establishing_smiles_mcs_search(smiles,
                                   smiles2,
                                   completion_smiles2,
                                   file_object):
    ''' Function: establishing_smiles_mcs_search
    This function establishes whether a mcs search needs to be done and if
    so which one for the SMILES
    Input: smiles - SMILES of the initial molecule being searched
           smiles2 - SMILES of the molecule being compared too
           completion_smiles2 - [] of the completed SMILES
           file_object - the file which writing too
    Output: None
    '''
    if smiles == smiles2:
        completion_smiles2.append(smiles)
        file_object.write("\t1")
    else:
        mol1 = Chem.MolFromSmiles(smiles)
        mol2 = Chem.MolFromSmiles(smiles2)
        mcs = rdFMCS.FindMCS([mol1, mol2],
                             bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                             timeout=1).smartsString
        mcs = Chem.MolFromSmarts(mcs)
        tcmcs = tanimoto_mcs(mol1.GetNumAtoms(),
                             mol2.GetNumAtoms(),
                             mcs.GetNumAtoms())
        completion_smiles2.append(smiles)
        file_object.write("\t%f" % round(tcmcs, 3))

    return None

def smiles_mcs(smiles_list):
    ''' Function: smiles_mcs
    This function performs operations to find all the SMILES mcs'
    Input: smiles_list - [] of str of SMILES
    Output: None
    '''
    try:
        smiles_df = pd.read_csv(args.output_smiles, sep='\t', header=0)
        complete_smiles = smiles_df['Unnamed: 0'].tolist()
        file_object = open(args.output_smiles, 'a')
    except FileNotFoundError:
        complete_smiles = []
        file_object = open(args.output_smiles, 'a')
        for smile in smiles_list:
            file_object.write("\t%s" % smile)
    completion_smiles2 = []

    length_of_reduced_graph_list = len(smiles_list)
    begin_smiles = time.time()
    for index, smiles in enumerate(smiles_list):
        print(index+1, ' out of ', length_of_reduced_graph_list, ' molecules')
        if smiles not in complete_smiles:
            file_object.write("\n%s" % smiles)
            for smiles2 in smiles_list:
                if smiles2 not in completion_smiles2 and smiles2 not in complete_smiles:
                    establishing_smiles_mcs_search(smiles,
                                                   smiles2,
                                                   completion_smiles2,
                                                   file_object)
                else:
                    file_object.write("\t")

    print('Time taken for SMILES to find similarity: ', time.time() - begin_smiles)
    file_object.close()

    return None

def establishing_rg_mcs_search(smiles,
                               smiles2,
                               completion_smiles2,
                               file_object):
    ''' Function: establishing_rg_mcs_search
    This function establishes whether a mcs search needs to be done and if
    so which one for the RG
    Input: smiles - SMILES of the initial molecule being searched
           smiles2 - SMILES of the molecule being compared too
           completion_smiles2 - [] of the completed SMILES
           file_object - the file which writing too
    Output: None
    '''
    params = Chem.SmilesParserParams()
    params.removeHs = True
    params.sanitize = False
    if smiles == smiles2:
        completion_smiles2.append(smiles)
        file_object.write("\t1")
    else:
        mol1 = Chem.MolFromSmiles(smiles, params)
        mol2 = Chem.MolFromSmiles(smiles2, params)
        mcs_smarts = rdFMCS.FindMCS([mol1, mol2],
                                    bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                                    timeout=1).smartsString
        smiles = ''.join(char for char in smiles if not char.isdigit())
        smiles2 = ''.join(char for char in smiles2 if not char.isdigit())
        smiles_nodes = smiles.replace(']', '').replace('=', '').replace('(', '').replace(')', '').replace('.', '').split('[')
        smiles2_nodes = smiles2.replace(']', '').replace('=', '').replace('(', '').replace(')', '').replace('.', '').split('[')
        smiles_nodes.remove('')
        smiles2_nodes.remove('')
        mcs_num_atoms = Chem.MolFromSmarts(mcs_smarts).GetNumAtoms()
        if mcs_num_atoms == 0:
            if bool([node for node in smiles_nodes if node in smiles2_nodes]) is True:
                mcs_num_atoms = 1
        tcmcs = tanimoto_mcs(mol1.GetNumAtoms(),
                             mol2.GetNumAtoms(),
                             mcs_num_atoms)
        completion_smiles2.append(smiles)
        file_object.write("\t%f" % round(tcmcs, 3))
    return None

def filling_matrix(output_file):
    df1 = pd.read_csv(output_file, sep='\t', header=0)
    ## reindexes the similarity matrix
    df2 = df1.set_index('Unnamed: 0')
    reindex_list = list(df2.columns.values)
    df3 = df2.reindex(reindex_list)
    
    ## need to fill the half of the matrix that is empty
    for idx in reindex_list:
        for idx2 in reindex_list:
            if pd.isnull(df3[idx][idx2]) == True:
                if idx == idx2:
                    df3[idx][idx2] = 1
                else:
                    df3[idx][idx2] = df3[idx2][idx]

    df3.to_csv(output_file, sep=args.separator)
    
    return None

def rg_mcs(reduced_graph_list):
    ''' Function: rg_mcs
    This function performs operations to find all the RG mcs'
    Input: reduced_graph_list - [] of str of RG
    Output: None
    '''
    reduced_graph_list = list(set(reduced_graph_list))
    try:
        smiles_df = pd.read_csv(args.output_rg, sep='\t', header=0)
        complete_smiles = smiles_df['Unnamed: 0'].tolist()
        file_object = open(args.output_rg, 'a')
    except FileNotFoundError:
        complete_smiles = []
        file_object = open(args.output_rg, 'a')
        for smile in reduced_graph_list:
            file_object.write("\t%s" % smile)
    completion_smiles2 = []

    length_of_reduced_graph_list = len(reduced_graph_list)
    begin_mcs = time.time()
    for index, smiles in enumerate(reduced_graph_list):
        print(index+1, ' out of ', length_of_reduced_graph_list, ' molecules')
        if smiles not in complete_smiles:
            file_object.write("\n%s" % smiles)
            for smiles2 in reduced_graph_list:
                if smiles2 not in completion_smiles2 and smiles2 not in complete_smiles:
                    establishing_rg_mcs_search(smiles,
                                               smiles2,
                                               completion_smiles2,
                                               file_object)
                else:
                    completion_smiles2.append(smiles)
                    file_object.write("\t")

    print('Time taken for RGs to find similarity: ', time.time() - begin_mcs)
    file_object.close()

    return None


if __name__ == '__main__':
    args = input_args()
    smiles_list, reduced_graph_list = read_in_file()

    if args.representationtype == 'smiles' or args.representationtype == 'both':
        smiles_mcs(smiles_list)
        filling_matrix(args.output_smiles)
    if args.representationtype == 'rg' or args.representationtype == 'both':
        rg_mcs(reduced_graph_list)
        filling_matrix(args.output_rg)
