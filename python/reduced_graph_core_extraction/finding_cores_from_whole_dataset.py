##############################################################################
# finding_cores_from_whole_dataset_with_full_checks
#
#
# Jess Stacey
##############################################################################

### find all cores from the whole data file 
### However doesn't check the size of the initial mcs

import argparse
import pandas as pd
import copy
from rdkit import Chem
from rdkit.Chem import rdFMCS

def input_args():
    parser = argparse.ArgumentParser(description='Arguments for the Reduced Graph Code')
    parser.add_argument('-r',
                        dest='reduced_graphs',
                        default=r'..\example_input_files\output_rg_file.txt',
                        help='Name of the input rg file')
    parser.add_argument('-m',
                        dest='mcs',
                        default=r'..\example_input_files\output_rg_rg_sim_matrix.txt',
                        help='Name of the input mcs or similarity matrix file')
    parser.add_argument('-o',
                        dest='output',
                        default=r'..\example_input_files\output_rg_core_extraction.txt',
                        type=str,
                        help='Name of the file containing all the information and values')
    parser.add_argument('-i',
                        dest='minsim',
                        default=0.1,
                        type=float,
                        help='Minimum similarity')
    parser.add_argument('-d',
                        dest='minnodes',
                        default=2,
                        type=int,
                        help='Minimum number of nodes as the core')
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

def finding_centroid(unrepresented_rgs,
                     rg_set,
                     similarity_matrix):
    ''' Function: finding_centroid
    This function takes a list of unrepresented_rgs and finds the rg that
    has the most neighbours and therefore is identified as the centroid.
    Input: unrepresented_rgs - list of rgs
        rg_set - set of the rg within the dataset
        similarity_matrix - pandas df of the similarity matrix
    Output: centroid - string of the RG that is the centroid
        mcs_smarts - string of SMARTS that represents the mcs
    '''
    print(len(unrepresented_rgs), ' Unrepresented rgs')
    params = Chem.SmilesParserParams()
    params.removeHs = True
    params.sanitize = False

    neighbour_dictionary = {}
    for rg in unrepresented_rgs:
        rg_neighbours = similarity_matrix[rg].tolist()
        rg_nearest_neighbours = [rg_n for rg_n in rg_neighbours if rg_n > args.minsim]
        
        if len(rg_nearest_neighbours) in neighbour_dictionary.keys():
            neighbour_dictionary[len(rg_nearest_neighbours)].append(rg)
        else:
            neighbour_dictionary[len(rg_nearest_neighbours)] = [rg]
    ### make a list of these centroids in the order
    ordered_list_of_centroids = []
    for key in sorted(neighbour_dictionary.keys(), reverse=True):
        ordered_list_of_centroids.extend(neighbour_dictionary[key])
    
    centroid_iteration = 0
    ### select the iteration number
    mcs_smarts_found = False
    neighbour_iteration_check = False
    while mcs_smarts_found == False:
        mcs_smarts_dict = {}
        neighbour_iteration = 0
        ### select the centroid
        centroid = ordered_list_of_centroids[centroid_iteration]
        ### find the furthest neighbour
        rg_neighbours = similarity_matrix[similarity_matrix[centroid] > args.minsim][centroid]
        ordered_rg_neighbours = rg_neighbours.sort_values(ascending=False).drop(centroid)
        ordered_rg_neighbours = ordered_rg_neighbours.to_frame()
        
        while neighbour_iteration_check == False:
            if neighbour_iteration == len(ordered_rg_neighbours):
                if len(ordered_rg_neighbours) == 0:
                    mcs_smarts = centroid
                    neighbour_iteration_check = True
                    mcs_smarts_found = True
                else:
                    centroid_iteration += 1
                    neighbour_iteration_check = True
            else:
                new_nearest_neighbour = list(ordered_rg_neighbours.index)[neighbour_iteration]
                
                mol1 = Chem.MolFromSmiles(centroid, params)
                mol2 = Chem.MolFromSmiles(new_nearest_neighbour, params)
                mcs_smarts = rdFMCS.FindMCS([mol1, mol2],
                                            bondCompare=rdFMCS.BondCompare.CompareOrderExact).smartsString
                if mcs_smarts == '':
                    neighbour_iteration += 1 
                ### check the size of the mcs_smarts
                elif Chem.MolFromSmarts(mcs_smarts).GetNumAtoms() >= args.minnodes:
                    neighbour_iteration_check = True
                    mcs_smarts_found = True
                else:
                    neighbour_iteration += 1
                    if mcs_smarts in mcs_smarts_dict.keys():
                        mcs_smarts_dict[mcs_smarts] += 1
                    else:
                        mcs_smarts_dict[mcs_smarts] = 1
            
            ### if a mcs cant be found for a centroid with certain length use 
            ### the one that has the most hits if all are blank go to the next centroid
            if mcs_smarts_found == False and neighbour_iteration == len(ordered_rg_neighbours):
                    print('None can be found for this centroid')
                    print('Currently just use the one that has the most hits')

                    core_nodes = ''
                    number_of_nodes = 0
                    for k, v in mcs_smarts_dict.items():
                        if k.count('[') > core_nodes.count('['):
                            core_nodes = k
                            number_of_nodes = v
                        elif k.count('[') >= core_nodes.count('['):
                            if v > number_of_nodes:
                                core_nodes = k
                                number_of_nodes = v
                    
                    mcs_smarts = core_nodes
                    ### if no mcs is found at all use the centroid
                    if mcs_smarts == '':
                        mcs_smarts = centroid
                    mcs_smarts_found = True    

    return(centroid, mcs_smarts)

def finding_core_rg(centroid,
                    min_neighbour_rg,
                    params):
    '''Function: finding_core_rg
    This function looks to find the MCS between the centroid and the nearest
    neighbour
    Input: centroid - str of SMILES of the centroid (RG)
           min_neighbour_rg - str of SMILES of the nearest neighbour RG
           params - rdkit parameters object
    Output: mcs_smarts - str of MCS SMARTS 
    '''
    ### Find current RG core
    mol1 = Chem.MolFromSmiles(centroid, params)
    mol2 = Chem.MolFromSmiles(min_neighbour_rg, params)
    mcs_smarts = rdFMCS.FindMCS([mol1, mol2],
                                bondCompare=rdFMCS.BondCompare.CompareOrderExact).smartsString

    return mcs_smarts

def checking_core(rg_set,
                  mcs_smarts,
                  centroid,
                  params):
    '''Function: checking_core
    This function searchs the rest of the set to see whether the mcs_smarts is
    present if not makes a subsmart of the existing smart if it is greater than
    or equal to the minnodes
    Input: rg_set - set of the rg within the dataset
           mcs_smarts - str of MCS SMARTS 
           centroid - str of SMILES of the centroid (RG)
           params - rdkit parameters object
    Output: Chem.MolToSmarts(mcs_smarts_mol) - str of MCS SMARTS
    '''
    mcs_smarts_mol = Chem.MolFromSmarts(mcs_smarts)
    for rg in rg_set:
        rg_mol = Chem.MolFromSmiles(rg, params)
        if rg_mol.HasSubstructMatch(mcs_smarts_mol):
            continue
        else:
            potential_new_mcs_smart = finding_core_rg(centroid, rg, params)
            ### checks that the MCS is above the minnodes number
            if Chem.MolFromSmarts(potential_new_mcs_smart).GetNumAtoms() >= args.minnodes:
                ### checks to see if the potential new MCS is a substructure of the old MCS
                if mcs_smarts_mol.HasSubstructMatch(Chem.MolFromSmarts(potential_new_mcs_smart)):
                    new_mcs_smart = mcs_smarts_mol.GetSubstructMatch(Chem.MolFromSmarts(potential_new_mcs_smart))
                    new_mcs_smart_mol = Chem.RWMol(mcs_smarts_mol)
                    for atom in range(new_mcs_smart_mol.GetNumAtoms()-1, -1, -1):
                        if atom not in new_mcs_smart:
                            new_mcs_smart_mol.RemoveAtom(atom)
                    if new_mcs_smart_mol.GetNumAtoms() >= args.minnodes:
                        mcs_smarts_mol = new_mcs_smart_mol
                ### find the MCS between the potential new MCS and the old MCS
                elif rdFMCS.FindMCS([mcs_smarts_mol, Chem.MolFromSmarts(potential_new_mcs_smart)], bondCompare=rdFMCS.BondCompare.CompareOrderExact).smartsString:
                    new_mcs_smart = rdFMCS.FindMCS([mcs_smarts_mol, Chem.MolFromSmarts(potential_new_mcs_smart)], bondCompare=rdFMCS.BondCompare.CompareOrderExact).smartsString
                    new_mcs_smart_mol = Chem.MolFromSmarts(new_mcs_smart)
                    if new_mcs_smart_mol.GetNumAtoms() >= args.minnodes:
                        mcs_smarts_mol = new_mcs_smart_mol
                ### is neither a substructure or contains a MCS
                else:
                    continue
            else:
                ### this MCS isn't big enough
                continue

    return Chem.MolToSmarts(mcs_smarts_mol)

def all_possible_rg_reps(core_list,
                         rg_set,
                         params):
    '''Function: all_possible_rg_reps
    Identifies which rg core are not represented by the cores within core_list 
    Input: core_list - list of rg cores
        rg_set - set of the rg
        params - rdkit parameters object
    Output: unrepresented_rgs - a list of the remaining RGs that do not contain a core within core_list 
    '''
    unrepresented_rgs = []
    for rg in rg_set:
        rg_mol = Chem.MolFromSmiles(rg, params)
        all_possible_rgs = []
        for core in core_list:
            if rg_mol.HasSubstructMatch(Chem.MolFromSmarts(core)):
                all_possible_rgs.append('Yes')
            else:
                all_possible_rgs.append('No')
        if 'Yes' not in set(all_possible_rgs):
            unrepresented_rgs.append(rg)
    return unrepresented_rgs

def core_cluster_generator(core_df,
                           rg_set,
                           similarity_matrix):
    '''Function: core_cluster_generator
    This function a core list and core dictionary of the cores that are present
    within the dataset. It also calls the other functions in order to help
    generate these.
    Input: core_df - currently an empty pandas df with columns RG and core
           rg_set - set of the rg within the dataset
           similarity_matrix - pandas df of the similarity matrix
    Output: core_df - pandas df of the RG and its corresponding core
    '''
    core_list = []
    core_dictionary = {}

    params = Chem.SmilesParserParams()
    params.removeHs = True
    params.sanitize = False

    unrepresented_rgs = copy.deepcopy(rg_set)

    while len(unrepresented_rgs) > 0:
        centroid, mcs_smarts = finding_centroid(unrepresented_rgs,
                                                      rg_set,
                                                      similarity_matrix)
        
        mcs_smarts = checking_core(unrepresented_rgs,
                                   mcs_smarts,
                                   centroid,
                                   params)
        core_list.append(mcs_smarts)

        new_mcs_smart = Chem.MolFromSmarts(mcs_smarts)
        for atom_index in range(0, new_mcs_smart.GetNumAtoms()):
            atom = new_mcs_smart.GetAtomWithIdx(atom_index)
            atom.SetNoImplicit(True)
        if new_mcs_smart not in core_dictionary.keys():
            core_dictionary[new_mcs_smart] = [mol for mol in unrepresented_rgs if Chem.MolFromSmiles(mol, params).HasSubstructMatch(Chem.MolFromSmarts(mcs_smarts))]
        unrepresented_rgs = all_possible_rg_reps(core_list,
                                                 unrepresented_rgs,
                                                 params)
        ### Filter the represented_rgs out of the similarity matrix (columns and rows)
        similarity_matrix = similarity_matrix[list(unrepresented_rgs)]
        similarity_matrix = similarity_matrix.filter(items=list(unrepresented_rgs), axis=0)

    for core, value in core_dictionary.items():
        for rg in value:
            core_df.loc[len(core_df)] = [rg, Chem.MolToSmiles(core)]
    
    return core_df

def core_generator(chemical_map_data):
    ''' Function: core_generator
    This function provides the correct information into the core_clsuter_generation
    function and merges together two dataframes
    Input: chemical_map_data - pandas df of all the molecules and their RG
    Output: core_df - pandas df of all the molecules and their RG and core RG
    '''
    similarity_matrix = pd.read_csv(args.mcs, sep='\t', header=0, index_col=0)
    core_df = pd.DataFrame(columns=['RG', 'core'])
    
    rg_list = chemical_map_data['RG'].dropna().tolist()
    rg_set = set(rg_list)
    core_df = core_cluster_generator(core_df,
                                     rg_set,
                                     similarity_matrix)

    return core_df

if __name__ == '__main__':
    args = input_args()
    
    merge_df = pd.read_csv(args.reduced_graphs, sep=args.separator, header=0)
    final_df = core_generator(merge_df)
    final_df.to_csv(args.output, sep=args.separator, index=False)
    print('All Cores: ', len(set(final_df['core'].tolist())))
