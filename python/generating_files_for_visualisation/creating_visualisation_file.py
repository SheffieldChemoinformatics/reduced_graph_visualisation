##############################################################################
# creating_visualisation_file_NEW
#
#
# Jess Stacey
##############################################################################
#### This finds the matching of all the possible cores for the molecules
### creating a new algorithm for establishing which atoms are part of the core or additional

import re
import argparse
import zlib
import base64
import ast
import pandas as pd
import numpy as np
from rdkit import Chem
import itertools
from rdkit.Chem import rdFMCS
from copy import deepcopy
from rdkit.Chem import Descriptors

class molecules_with_multiple_possible_mappings:
    def __init__(self):
        self.list_of_mols_part1 = []
        self.list_of_mols_part2 = []
        self.list_of_mols_part3 = []
        self.list_of_mols_part4 = []
        
        self.molecules_with_multiple_mappings = []

def input_args():
    parser = argparse.ArgumentParser(description='Arguments for the Creating Visualisation File')
    parser.add_argument('-r',
                        dest='reduced_graphs',
                        default=r'..\example_input_files\output_reduced_graphs_file.txt',
                        help='Name of the input rg file')
    parser.add_argument('-s',
                        dest='sdf',
                        default=r'..\example_input_files\output_reduced_graphs_file.sdf',
                        help='Name of the input sdf file')
    parser.add_argument('-a',
                        dest='activities',
                        default=r'..\example_input_files\pIC50_data.txt',
                        help='Name of the input file that contains the activities of the molecules')
    parser.add_argument('-c',
                        dest='cores',
                        default=r'..\example_input_files\output_rg_core_extraction.txt',
                        help='Name of the input clusters file')
    parser.add_argument('-o',
                        dest='output',
                        default=r'..\example_input_files\output_node_information.txt',
                        type=str,
                        help='Name of the file containing all the information and values')
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


def generate_compressed_sdf():
    molecules_dict = {}
    with open(args.sdf) as file:
        molecule = []
        for line in file:
            if '$$$$' not in line:
                molecule.append(line)
            else:
                molecule.append(line)
                mol_block = ''.join(molecule)
                mol_block_compress = zlib.compress(mol_block.encode('utf-8'))
                molecules_dict[molecule[0].replace('\n','')] = str(base64.b64encode(mol_block_compress), 'utf-8')
                molecule = []

    return molecules_dict

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
            
def extracting_all_possible_cores(cores,
                                  new_df):
    ###extract out all the cores and then iterate through the cores to see whether any are present within each RG
    list_of_cores = set(cores['core'].tolist())
    params = Chem.SmilesParserParams()
    params.removeHs = True
    params.sanitize = False
    all_rg_core_matches = {}
    for i, row in cores.iterrows():
        rg_mol = Chem.MolFromSmiles(row['RG'], params)
        core_matches = []
        for core in list_of_cores:
            core_mol = Chem.MolFromSmiles(core, params)
            if rg_mol.HasSubstructMatch(core_mol):
                core_matches.append(core)
        all_rg_core_matches[row['RG']] = core_matches
    
    final_df = pd.DataFrame(columns=['SMILES', 'ID', 'RG', 'SDF', 'pIC50', 'core'])
    
    for i, row in new_df.iterrows():
        try:
            core = all_rg_core_matches[row['RG']]
            for c in core:
                new_row = row.to_dict()
                new_row['core'] = c
                ### row to append
                final_df = final_df.append(new_row, ignore_index=True)
        except:
            new_row = row.to_dict()
            new_row['core'] = row['RG']
            final_df = final_df.append(new_row, ignore_index=True)
    
    return final_df
            
def read_in_files():
    ### reduced graphs
    reduced_graphs = pd.read_csv(args.reduced_graphs, header=0, sep=args.separator)
    ### sdf
    sdf_molecules_dict = generate_compressed_sdf()
    
    for r, v in reduced_graphs.iterrows():
        if v['SMILES'] in sdf_molecules_dict.keys():
            reduced_graphs.at[r, 'SDF'] = sdf_molecules_dict[v['SMILES']]
        else:
            print(v['SMILES'], v['ID'])
            print('ERROR! SMILES not present in SDF')

    activities = pd.read_csv(args.activities, header=0, sep=args.separator)
    ###filter so just ID and activities as pIC50s
    activities = activities[['ID', 'pIC50']]
    new_df = pd.merge(reduced_graphs, activities, on=['ID'], how='left') 
    
    ### cores
    cores = pd.read_csv(args.cores, header=0, sep=args.separator)
    
    ###remove cluster and core_x if there and then rename the columns
    if 'cluster' in cores.columns:
        del cores['cluster']
    if 'core_x' in cores.columns:
        del cores['core_x']
        cores.rename(columns={'core_y':'core'}, inplace=True)
    if 'number_of_heavy_atoms' in cores.columns:
        del cores['number_of_heavy_atoms']
    
    final_df = extracting_all_possible_cores(cores,
                                             new_df)   

    return final_df

def core_edges_and_isotopes(core_list,
                            params):
    core_numbered = {}
    core_edges = {}
    for core in core_list:
        core_rg_mol = Chem.MolFromSmiles(core, params)
        core_rg_mol_edited = Chem.Mol(core_rg_mol)
        
        separated_rg = {item: '{0}{1}'.format(item, core_rg_mol_edited.GetAtomWithIdx(item).GetSmarts().replace('[', '').replace(']', '')) for item in range(0, core_rg_mol_edited.GetNumAtoms())}
        core_numbered[core] = separated_rg
        
        edges = []
        for atom in range(0, core_rg_mol_edited.GetNumAtoms()):
            for neigh in [bond.GetOtherAtomIdx(atom) for bond in core_rg_mol_edited.GetAtomWithIdx(atom).GetBonds()]:
                if neigh > atom:
                    bond_type = core_rg_mol_edited.GetBondBetweenAtoms(neigh, atom).GetBondType()
                    if bond_type == Chem.rdchem.BondType.SINGLE:
                        edges.append({"source": atom,
                                      "target": neigh,
                                      "weight": 1})
                    elif bond_type == Chem.rdchem.BondType.DOUBLE:
                        edges.append({"source": atom,
                                      "target": neigh,
                                      "weight": 2})
                    elif bond_type == Chem.rdchem.BondType.TRIPLE:
                        edges.append({"source": atom,
                                      "target": neigh,
                                      "weight": 3})
        core_edges[core] = edges
    return core_numbered, core_edges

def r_group_assignment(fragmented_removed_core):
    r_group_dictionary = {}
    if len(fragmented_removed_core) == 1:
        r_group_dictionary[Chem.MolToSmiles(fragmented_removed_core[0])] = 1
    else:
        r_group_neighbours = {}
        for frag in fragmented_removed_core:
            neighbour = [int(neigh.strip('[').strip('*]')) for neigh in re.findall(r'\d+\*', Chem.MolToSmiles(frag))]
            ### * isn't getting recognised in this so need to add a zero if it is present
            if '*[' in Chem.MolToSmiles(frag):
                neighbour.append(0)
            for neigh in neighbour:
                if neigh not in r_group_neighbours.keys():
                    r_group_neighbours[neigh] = [Chem.MolToSmiles(frag)]
                else:
                    r_group_neighbours[neigh].append(Chem.MolToSmiles(frag))
            if not neighbour:
                if 'none' not in r_group_neighbours:
                    r_group_neighbours['none'] = [Chem.MolToSmiles(frag)]
                else:
                    r_group_neighbours['none'].append(Chem.MolToSmiles(frag))

        ### inspect the values of the dictionary if any are unique assign 1 if not need to establish a hierarchy
        for key, value in r_group_neighbours.items():
            if len(value) == 1:
                r_group_dictionary[value[0]] = 1
            else:
                for group in value:
                    r_group_dictionary[group] = value.index(group)+1
        
    return r_group_dictionary
    
def extracting_r_groups(rg_core_indexes,
                        rg,
                        core,
                        smiles,
                        fragmented_removed_core):
    params = Chem.SmilesParserParams()
    params.removeHs = True
    params.sanitize = False
    r_groups = {}
    r_group_indexed = r_group_assignment(fragmented_removed_core)
    for frag in fragmented_removed_core:
        r_group_smarts = []
        r_group_atom_indexes = []
        r_group_dictionary = {}
        for f in range(0, frag.GetNumAtoms()):
            if frag.GetAtomWithIdx(f).HasProp('rg_idx'):
                rg_idx = int(frag.GetAtomWithIdx(f).GetProp('rg_idx'))
                atom = rg.GetAtomWithIdx(rg_idx)
                if atom.HasProp('molFileValue'):
                    node_smarts = atom.GetProp('molFileValue').split(' ')
                    indexes = node_smarts[1:-1]
                    smiles_atom_indexes = [int(s.strip('[').strip(']').strip("'").strip("',")) for s in indexes]
                    r_group_atom_indexes.extend(smiles_atom_indexes)
                    r_group_smarts.append(node_smarts[-1])
                r_group_dictionary[(f, rg.GetAtomWithIdx(rg_idx).GetSmarts())] = node_smarts[-1]
            else:
                r_group_dictionary[(f, '[*]')] = '[*]'
        
        r_groups[tuple(r_group_atom_indexes)] = [Chem.MolToSmiles(frag), r_group_indexed[Chem.MolToSmiles(frag)], r_group_smarts, r_group_dictionary]
    
    return r_groups

def fragmentation(smiles,
                  core_atom_indexes):
    ### fragment the molecule into core and R groups with the corresponding wild card atom indexes
    ### find the bond between the core and R groups and fragment the molecule
    fragment_on_bonds = []
    for atom in core_atom_indexes:
        neighbours = [neigh.GetIdx() for neigh in smiles.GetAtomWithIdx(atom).GetNeighbors()  if neigh.GetIdx() not in core_atom_indexes]
        ###Get Bond between the two atoms
        for atom2 in neighbours:
            bond_idx = smiles.GetBondBetweenAtoms(atom, atom2)
            fragment_on_bonds.append(bond_idx.GetIdx())  
    ### need to set the fragment_on_bonds too!
    ### if dotted smiles and the core is the whole of one dot then fragmented_smiles doesnt work
    if fragment_on_bonds:
        fragmented_smiles = Chem.FragmentOnBonds(smiles, set(fragment_on_bonds))
    else:
        fragmented_smiles = smiles
    
    ###remove the core from the fragmented_smiles
    remove_core_indexes = []
    for atom in core_atom_indexes:
        neighbours = [neigh.GetIdx() for neigh in fragmented_smiles.GetAtomWithIdx(atom).GetNeighbors()  if neigh.GetIdx() not in core_atom_indexes]
        remove_core_indexes.append(atom)
        remove_core_indexes.extend(neighbours)
    fragmented_smiles = Chem.RWMol(fragmented_smiles)
    
    ### need to make a set so that if something has multiple neighbours 
    remove_core_indexes = set(remove_core_indexes)
    for atom in sorted(remove_core_indexes, reverse=True):
        fragmented_smiles.RemoveAtom(atom)

    ### fragment the fragments
    try:
        fragments_mols = Chem.GetMolFrags(fragmented_smiles, asMols=True)
    except:
        fragments_mols = Chem.GetMolFrags(fragmented_smiles, asMols=True, sanitizeFrags=False)
        
    return fragments_mols

def processing_additional_information(frag,
                                      frag_index,
                                      rg_r_groups):
    group_dict = {}
    group_dict['RG'] = rg_r_groups[frag_index][0]
    group_dict['SMARTS'] = Chem.MolToSmiles(frag)
    ### do i need to add a 0 
    group_dict['neighbours'] = [int(neigh.strip('[').strip('*]')) for neigh in re.findall(r'\d+\*', rg_r_groups[frag_index][0])]
    group_dict['atom_indexes'] = list(frag_index)
    group_dict['r_group'] = rg_r_groups[frag_index][1]
    group_dict['name'] = ''.join([str(group_dict['neighbours']), str(group_dict['r_group'])])
    group_dict['individual_smarts'] = rg_r_groups[frag_index][2]
    group_dict['breakdown'] = rg_r_groups[frag_index][3]
    
    return group_dict

def canonicalise_smarts(smarts):
    mol = Chem.MolFromSmiles(smarts, sanitize=False)

    return Chem.MolToSmiles(mol)

def extracting_more_meta_data(smiles,
                              molecule,
                              core_smiles,
                              core,
                              row):
    ### Gives core smarts
    core_smarts = Chem.ReplaceSidechains(smiles, core_smiles)
    removalsmiles = Chem.RWMol(smiles)

    ### add meta data into the removalsmiles of the atom
    for atom_idx in range(0, removalsmiles.GetNumAtoms()):
        atom = removalsmiles.GetAtomWithIdx(atom_idx)
        atom.SetProp('atom_idx', str(atom_idx))

    additionals_list=[]    
    ###If the core does not equal the full RG find the additional groups
    if row['core'] != row['RG']:
        fragments_mols_r_groups = fragmentation(removalsmiles,
                                                row['core_atom_indexes'])#smiles_core_atom_indexes
        
        fragments_rg_r_groups = fragmentation(molecule,
                                              row['core_node_indexes'])#rg_core_atom_indexes
        ### Gives R groups      
        rg_r_groups = extracting_r_groups(row['core_node_indexes'],#rg_core_atom_indexes
                                          molecule,
                                          core,
                                          molecule,
                                          fragments_rg_r_groups)
        
        for frag in fragments_mols_r_groups:
            fragment_atom_indexes = []
            ### select out the fragment indexes in the original molecule
            for atom in range(0, frag.GetNumAtoms()):
                if frag.GetAtomWithIdx(atom).GetSymbol() != '*':
                    fragment_atom_indexes.append(int(frag.GetAtomWithIdx(atom).GetProp('atom_idx')))
            
                      
            ### search for key in the fragments list
            fragment_indexes = [idxs for idxs in rg_r_groups.keys() if set(fragment_atom_indexes).issuperset(set(idxs))]
            if len(fragment_indexes) == 1:
                frag_index = fragment_indexes[0]
                group_dict = processing_additional_information(frag,
                                                               frag_index,
                                                               rg_r_groups)
                additionals_list.append(group_dict)
            elif len(fragment_indexes) == 0:
                ### this happens when the node hasn't been defined in the parameters of the RG
                print('not adding as is not necessary for this parameterisation')
            else: 
                print('ERROR')
                print(fragment_indexes)
                print(fragment_atom_indexes, 'frag indexes')
                print(rg_r_groups.keys(), 'keys')
                
    return core_smarts, additionals_list

def extracting_core_meta_data_and_lengths(smiles,
                                          molecule,
                                          rg_core_atom_indexes):
    core_smiles_atom_indexes = []
    core_info = {}
    ### create a list of atom indexes and the node that they lie in
    molecule_node_atom_indexes = [-1]*smiles.GetNumAtoms()
    core_node_idx = {}
    
    for core_idx, atom_idx in enumerate(rg_core_atom_indexes):
        atom = molecule.GetAtomWithIdx(atom_idx)
        node_smarts = atom.GetProp('molFileValue').split(' ')
        
        indexes = node_smarts[1:-1]
        smiles_atom_indexes = [int(s.strip('[').strip(']').strip("'").strip("',")) for s in indexes]
        core_smiles_atom_indexes.extend(smiles_atom_indexes)
        core_info[core_idx] = node_smarts[-1]
        
        core_node_idx[core_idx] = smiles_atom_indexes
        ### understanding distances to other core nodes
        for val in smiles_atom_indexes:
            molecule_node_atom_indexes[val] = core_idx
    
    ##get subsitution nodes ones that aren't art of the core
    substitution_atom_indexes = []
    for val in core_smiles_atom_indexes:
        neighbours = [v.GetIdx() for v in smiles.GetAtomWithIdx(val).GetNeighbors() if v.GetIdx() not in core_smiles_atom_indexes]
        substitution_atom_indexes.extend(neighbours)
    
    ### understanding distances to other core nodes
    length_between_nodes = {}
    all_substitution_topological_distances = []
    for node_idx in range(0, len(rg_core_atom_indexes)):
        node_smiles_idx = core_node_idx[node_idx]
        for node_idx2 in range(node_idx+1, len(rg_core_atom_indexes)):
            shortest_path_length = None
            node_smiles_idx2 = core_node_idx[node_idx2]
            
            combo_list = itertools.product(node_smiles_idx, node_smiles_idx2)
            distance_matrix = Chem.GetDistanceMatrix(smiles)
            for combo in combo_list:
                sp = int(distance_matrix[combo[0]][combo[1]])
                if shortest_path_length == None:
                    shortest_path_length = sp
                elif sp < shortest_path_length:
                    shortest_path_length = sp
            if shortest_path_length == None:
                shortest_path_length = 0
            length_between_nodes['{0}-{1}'.format(node_idx, node_idx2)] = shortest_path_length
            
        subtitution_topological = []
        distance_matrix = Chem.GetDistanceMatrix(smiles)
        for sub_val in substitution_atom_indexes:
            shortest_substituiton_length = None
            for core_val in node_smiles_idx:
                sp = int(distance_matrix[core_val][sub_val])
                if shortest_substituiton_length == None:
                    shortest_substituiton_length = sp
                elif sp < shortest_substituiton_length:
                    shortest_substituiton_length = sp
            subtitution_topological.append(shortest_substituiton_length)
        all_substitution_topological_distances.append(subtitution_topological)
                        
    return length_between_nodes, core_smiles_atom_indexes, core_info, all_substitution_topological_distances

def process_single_mappings(one_mappings_df,
                            params):
    core_length_mappings = {}
    core_smarts_mappings = {}
    core_substitution_info = {}
    
    for row_idx, row in one_mappings_df.iterrows():
        smiles = Chem.MolFromSmiles(row['SMILES'])
        core = Chem.MolFromSmiles(row['core'], params)
    
        ### decompress the sdf
        mol_block = zlib.decompress(base64.b64decode(row['SDF']))
        molecule = Chem.MolFromMolBlock(mol_block, sanitize=False, strictParsing=False, removeHs=True)
    
        ### remove the hydrogens from each rg node and set meta data (indexes)
        for atomI in range(0, molecule.GetNumAtoms()):
            molecule.GetAtomWithIdx(atomI).SetNoImplicit(True)
            molecule.GetAtomWithIdx(atomI).SetProp('rg_idx', str(atomI))
    
        rg_core_atom_indexes = row['node_idx_of_core_potential_mappings'][0]
        length_between_nodes, core_smiles_atom_indexes, core_info, all_substitution_topological_distances = extracting_core_meta_data_and_lengths(smiles,
                                                                                                                                                  molecule,
                                                                                                                                                  rg_core_atom_indexes)
        
        one_mappings_df.at[row_idx, 'length_between_core_nodes'] = length_between_nodes
        if row['core'] not in core_length_mappings.keys():
            core_length_mappings[row['core']] = {str(length_between_nodes): 1}
            core_substitution_info[row['core']] = {str(length_between_nodes): [all_substitution_topological_distances]}
        else:
            if str(length_between_nodes) not in core_length_mappings[row['core']]:
                core_length_mappings[row['core']][str(length_between_nodes)] = 1
                core_substitution_info[row['core']][str(length_between_nodes)] = [all_substitution_topological_distances]
            else:
                core_length_mappings[row['core']][str(length_between_nodes)] += 1
                if all_substitution_topological_distances not in core_substitution_info[row['core']][str(length_between_nodes)]:
                    core_substitution_info[row['core']][str(length_between_nodes)].append(all_substitution_topological_distances)
                
            
        core_smiles = Chem.RWMol(smiles)
        not_core_atoms = [i for i in range(0, (core_smiles.GetNumAtoms())) if i not in core_smiles_atom_indexes]
        ### remove the atoms that aren't part of the core
        for atom_idx in sorted(not_core_atoms, reverse=True):
            core_smiles.RemoveAtom(atom_idx)
        
        one_mappings_df.at[row_idx, 'core_node_indexes'] = rg_core_atom_indexes
        one_mappings_df.at[row_idx, 'core_atom_indexes'] = core_smiles_atom_indexes
    
    
        core_smarts, additionals_list = extracting_more_meta_data(smiles,
                                                                  molecule,
                                                                  core_smiles,
                                                                  core,
                                                                  row)
        
        ### Need to canonicalise the core_smarts
        one_mappings_df.at[row_idx, 'core_smarts'] = canonicalise_smarts(Chem.MolToSmiles(core_smarts))
        one_mappings_df.at[row_idx, 'core_info'] = core_info
        one_mappings_df.at[row_idx, 'additional'] = additionals_list
        
        if row['core'] not in core_smarts_mappings.keys():
            core_smarts_mappings[row['core']] = {one_mappings_df.at[row_idx, 'core_smarts']: core_info}
        else:
            if one_mappings_df.at[row_idx, 'core_smarts'] not in core_smarts_mappings[row['core']].keys():
                core_smarts_mappings[row['core']][one_mappings_df.at[row_idx, 'core_smarts']] = core_info
    return core_length_mappings, core_smarts_mappings, core_substitution_info

class extracting_mapped_rg_core():
    def __init__(self):
        self.multi_mappings_df = None
        self.core_length_mappings = []
        self.core_substitution_info = []
        self.core_smarts_mapping = []
        
        self.row_idx = ''
        self.row = ''
        self.smiles = ''
        self.molecule = ''
        self.core = ''
        self.all_lengths_between_nodes = []
        self.all_core_smiles= []
        self.all_sub_topo_distances = []

        self.params = Chem.SmilesParserParams()
        self.params.removeHs = True
        self.params.sanitize = False
        
        self.to_further_examine = 0
        self.topological_map_resolved = 0
        self.topological_sub_resolved = 0
        self.check_chemical_graph_resolved = 0
        self.hac_resolved = 0
        self.mw_resolved = 0
        self.sub_top_resolved = 0
        self.sub_top_avg_resolved = 0
        self.arbitary_resolved = 0
        
        self.no_initial_edit_distance = 0
        self.no_initial_avg_edit_distance = 0
        self.no_initial_topological_sub_resolved = 0
        self.no_initial_check_chemical_graph_resolved = 0
        self.no_initial_hac_resolved = 0
        self.no_initial_mw_resolved = 0
        self.no_initial_sub_top_resolved = 0
        self.no_initial_sub_top_avg_resolved = 0
        self.no_initial_arbitary_resolved = 0
        
        self.check_chemical_graph_resolved2 = 0
        self.no_initial_check_chemical_graph_resolved2 = 0
        
        return None
    
    def set_value_of_single_mappings(self,
                                     multi_mappings_df,
                                     core_length_mappings,
                                     core_substitution_info,
                                     core_smarts_mappings):
        self.multi_mappings_df = multi_mappings_df
        self.core_length_mappings = core_length_mappings
        self.core_substitution_info = core_substitution_info
        self.core_smarts_mappings = core_smarts_mappings
        return None
        
    def set_values(self,
                   row,
                   smiles,
                   molecule,
                   core,
                   all_lengths_between_nodes,
                   all_core_smiles,
                   all_sub_topo_distances):
        self.row = row
        self.smiles = smiles
        self.molecule = molecule
        self.core = core
        self.all_lengths_between_nodes = all_lengths_between_nodes
        self.all_core_smiles = all_core_smiles
        self.all_sub_topo_distances = all_sub_topo_distances
        return None
    
    def print_stats(self):
        print('With An Initial Topological Node Map Match:')
        print('Topo Node Match Resolved: ', self.topological_map_resolved)
        print('Topo Substitution Resolved: ', self.topological_sub_resolved)
        print('CG Resolved: ', self.check_chemical_graph_resolved)
        print('HAC Resolved: ', self.hac_resolved)
        print('MW Resolved: ', self.mw_resolved)
        print('Sub Top Resolved: ', self.sub_top_resolved)
        print('Sub Top Avg Resolved: ', self.sub_top_avg_resolved)
        print('Second CG Resolved: ', self.check_chemical_graph_resolved2)
        print('Arbitary Choice: ', self.arbitary_resolved)
        
        print('\nWith No Initial Topological Node Map Matches:')
        print('Topo Node Min Edit Distance Resolved: ', self.no_initial_edit_distance)
        print('Topo Node MinAvg Edit Distance Resolved: ', self.no_initial_avg_edit_distance)
        print('Topo Substitution Resolved: ', self.no_initial_topological_sub_resolved)
        print('CG Resolved: ', self.no_initial_check_chemical_graph_resolved)
        print('HAC Resolved: ', self.no_initial_hac_resolved)
        print('MW Resolved: ', self.no_initial_mw_resolved)
        print('Sub Top Resolved: ', self.no_initial_sub_top_resolved)
        print('Sub Top Avg Resolved: ', self.no_initial_sub_top_avg_resolved)
        print('Second CG Resolved: ', self.no_initial_check_chemical_graph_resolved2)
        print('Arbitary Choice: ', self.no_initial_arbitary_resolved)
        
        print('\nProblem Molecules: ', self.to_further_examine)       
        return None
       
    def compare_structs(self,
                        substructures_to_compare,
                        substructures_to_examine,
                        match_coming_from):
        ### Examine if structures are the same 
        if len(substructures_to_compare) == 2:
            if (Chem.MolFromSmiles(substructures_to_compare[0], self.params).HasSubstructMatch(Chem.MolFromSmiles(substructures_to_compare[1], self.params)) and Chem.MolFromSmiles(substructures_to_compare[1], self.params).HasSubstructMatch(Chem.MolFromSmiles(substructures_to_compare[0], self.params))) == True:
                idx_to_use = substructures_to_examine[0]
                problem_label = 'LastValsSame'
            else:
                if match_coming_from == 'topo_node_match':
                    self.arbitary_resolved += 1
                    problem_label = 'SortedViaArbitaryChose'
                elif match_coming_from == 'no_topo_node_match':
                    self.no_initial_arbitary_resolved += 1
                    problem_label = 'SolvedNoThenSortedViaArbitaryChose'
                idx_to_use = substructures_to_examine[0]
        else:
            idx_to_use = 0
            problem_label = 'OneToCompare'
            
        return idx_to_use, problem_label
    
    def compare_chemical_graph2(self,
                                substructures_to_compare,
                                substructures_to_examine,
                                match_coming_from):    
        if len(set(substructures_to_compare)) == 1:
            if match_coming_from == 'topo_node_match':
                problem_label = 'SortedViaCG2'
                self.check_chemical_graph_resolved2 +=1 
            elif match_coming_from == 'no_topo_node_match':
                problem_label = 'SolvedNoThenSortedViaCG2'
                self.no_initial_check_chemical_graph_resolved2 +=1
            idx_to_use = substructures_to_examine[0]
        else:
            substructures_to_compare2 = []
            substructures_to_examine2 = []
            
            for substruc in set(substructures_to_compare):
                idx = substructures_to_compare.index(substruc)
                substructures_to_compare2.append(substruc)
                substructures_to_examine2.append(substructures_to_examine[idx])
            
            idx_to_use, problem_label = self.compare_structs(substructures_to_compare2,
                                                             substructures_to_examine2,
                                                             match_coming_from)
        
        return idx_to_use, problem_label

    def comparing_substructures(self,
                                substructures_to_compare,
                                substructures_to_examine,
                                match_coming_from):
        max_mcs = []
        average_mcs = []
        for sub_idx, substructure1 in enumerate(substructures_to_compare):                                        
            core_smiles_atom_indexes = self.all_core_smiles[substructures_to_examine[sub_idx]]
            core_smiles = Chem.RWMol(self.smiles)
            not_core_atoms = [i for i in range(0, (core_smiles.GetNumAtoms())) if i not in core_smiles_atom_indexes]
            ### remove the atoms that aren't part of the core
            for atom_idx in sorted(not_core_atoms, reverse=True):
                core_smiles.RemoveAtom(atom_idx)
            core_smarts, additionals_list = extracting_more_meta_data(self.smiles,
                                                                      self.molecule,
                                                                      core_smiles,
                                                                      self.core,
                                                                      self.row)
            
            mcs_list = []
            for substructure2 in self.core_smarts_mappings[self.row['core']].keys():
                substructure2_mol = Chem.MolFromSmiles(substructure2, self.params)
                mcs = rdFMCS.FindMCS([core_smarts, substructure2_mol],
                                     bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                                     timeout=1).smartsString
                mcs = Chem.MolFromSmarts(mcs)
                tcmcs = tanimoto_mcs(core_smarts.GetNumAtoms(),
                                     substructure2_mol.GetNumAtoms(),
                                     mcs.GetNumAtoms())
                mcs_list.append(tcmcs)
            max_mcs.append(max(mcs_list))
            average_mcs.append(np.mean(mcs_list))

        if max_mcs.count(max(max_mcs)) == 1:
            idx_to_use = substructures_to_examine[max_mcs.index(max(max_mcs))]
            if match_coming_from == 'topo_node_match':
                problem_label = 'SortedViaTopSub'
                self.sub_top_resolved +=1
            elif match_coming_from == 'no_topo_node_match':
                problem_label = 'SolvedNoThenSortedViaTopSub'
                self.no_initial_sub_top_resolved +=1
        elif average_mcs.count(max(average_mcs)) == 1:
            idx_to_use = substructures_to_examine[average_mcs.index(max(average_mcs))]
            if match_coming_from == 'topo_node_match':
                problem_label = 'SortedViaTopAvgSub'
                self.sub_top_avg_resolved +=1
            elif match_coming_from == 'no_topo_node_match':
                problem_label = 'SolvedNoThenSortedViaTopAvgSub'
                self.no_initial_sub_top_avg_resolved +=1
        else:
            substructures_to_compare = [substructures_to_compare[i] for i, v in enumerate(average_mcs) if v == max(average_mcs)]
            substructures_to_examine = [substructures_to_examine[i] for i, v in enumerate(average_mcs) if v == max(average_mcs)]
            idx_to_use, problem_label = self.compare_chemical_graph2(substructures_to_compare,
                                                                     substructures_to_examine,
                                                                     match_coming_from)
            
        return idx_to_use, problem_label

    def comparing_molecular_weight(self,
                                   substructures_to_compare,
                                   substructures_to_examine,
                                   match_coming_from):
        mw_values = [Descriptors.ExactMolWt(Chem.MolFromSmiles(substructure, self.params)) for substructure in substructures_to_compare]
        max_mw_list = [(t == max(mw_values)) for t in mw_values]
        if max_mw_list.count(True) == 1:
            if match_coming_from == 'topo_node_match':
                problem_label = 'SortedViaMW'
                self.mw_resolved +=1
            elif match_coming_from == 'no_topo_node_match':
                problem_label = 'SolvedNoThenSortedViaMW'
                self.no_initial_mw_resolved +=1
            idx_to_use = substructures_to_examine[mw_values.index(max(mw_values))]
        else:
            ### look at substructure
            substructures_to_compare = [substructures_to_compare[i] for i, v in enumerate(max_mw_list) if v == True]
            substructures_to_examine = [substructures_to_examine[i] for i, v in enumerate(max_mw_list) if v == True]
            idx_to_use, problem_label = self.comparing_substructures(substructures_to_compare,
                                                                     substructures_to_examine,
                                                                     match_coming_from)
            
        return idx_to_use, problem_label

    def comparing_heavy_atom_count(self,
                                   substructures_to_compare,
                                   substructures_to_examine,
                                   match_coming_from):
        hac_values = [Chem.MolFromSmiles(substructure, self.params).GetNumAtoms() for substructure in substructures_to_compare]
        max_hac_list = [(t == max(hac_values)) for t in hac_values]
        if max_hac_list.count(True) == 1:
            if match_coming_from == 'topo_node_match':
                problem_label = 'SortedViaHAC'
                self.hac_resolved +=1
            elif match_coming_from == 'no_topo_node_match':
                problem_label = 'SolvedNoThenSortedViaHAC'
                self.no_initial_hac_resolved +=1
            idx_to_use = substructures_to_examine[hac_values.index(max(hac_values))]
        else:
            ### Look at MW
            ### Filter ones so only examining ones with highest HAC
            substructures_to_compare = [substructures_to_compare[i] for i, v in enumerate(max_hac_list) if v == True]
            substructures_to_examine = [substructures_to_examine[i] for i, v in enumerate(max_hac_list) if v == True]
            idx_to_use, problem_label = self.comparing_molecular_weight(substructures_to_compare,
                                                                        substructures_to_examine,
                                                                        match_coming_from)
        
        return idx_to_use, problem_label

    def comparing_chemical_graph(self,
                                 substructures_to_examine,
                                 match_coming_from):
        substructures_to_compare = []
        for sub in substructures_to_examine:
            core_smiles_atom_indexes = self.all_core_smiles[sub]
            core_smiles = Chem.RWMol(self.smiles)
            not_core_atoms = [i for i in range(0, (core_smiles.GetNumAtoms())) if i not in core_smiles_atom_indexes]
            ### remove the atoms that aren't part of the core
            for atom_idx in sorted(not_core_atoms, reverse=True):
                core_smiles.RemoveAtom(atom_idx)
            substructures_to_compare.append(Chem.MolToSmiles(core_smiles, 1))
     
        if len(set(substructures_to_compare)) == 1:
            if match_coming_from == 'topo_node_match':
                problem_label = 'SortedViaCG'
                self.check_chemical_graph_resolved +=1 
            elif match_coming_from == 'no_topo_node_match':
                problem_label = 'SolvedNoThenSortedViaCG'
                self.no_initial_check_chemical_graph_resolved +=1 
            idx_to_use = substructures_to_examine[0]
        else:
            ### HAC
            idx_to_use, problem_label = self.comparing_heavy_atom_count(substructures_to_compare,
                                                                        substructures_to_examine,
                                                                        match_coming_from)
        return idx_to_use, problem_label

    def topological_substitutions(self,
                                  indexes_to_compare,
                                  match_coming_from):
        all_number_of_changes_to_make = []
        number_of_changes_to_make = []
        for compare_idx in indexes_to_compare:
            this_substitution = self.all_sub_topo_distances[compare_idx]
            if match_coming_from == 'topo_node_match':
                substitutions_to_compare_to = self.core_substitution_info[self.row['core']][str(self.all_lengths_between_nodes[compare_idx])]
            elif match_coming_from == 'no_topo_node_match':
                previous_seen = list(self.core_substitution_info[self.row['core']].values())
                substitutions_to_compare_to = [item for sub in previous_seen for item in sub]
            
            all_number_of_changes = []
            lowest_number_of_changes = None
            for sub in substitutions_to_compare_to:
                changes_to_make = 0
                for node_idx in range(0, len(sub)):
                    a=this_substitution[node_idx]
                    b=sub[node_idx]
                    nodes_changes_to_make = 0
                    for s in set(a):
                        count_a = a.count(s)
                        count_b = b.count(s)
                        if count_a != count_b:
                            nodes_changes_to_make += abs(count_a - count_b)
                    for s in set(b):
                        if s not in a:
                            nodes_changes_to_make += b.count(s)
                    changes_to_make += nodes_changes_to_make
                    all_number_of_changes.append(changes_to_make)
                if lowest_number_of_changes == None:
                    lowest_number_of_changes = changes_to_make
                elif lowest_number_of_changes > changes_to_make:
                    lowest_number_of_changes = changes_to_make
            number_of_changes_to_make.append(lowest_number_of_changes)
            all_number_of_changes_to_make.append(all_number_of_changes)
        min_value_topo = min(number_of_changes_to_make)
        min_topo_list = [(t == min_value_topo) for t in number_of_changes_to_make]
        if min_topo_list.count(True) == 1:
            if match_coming_from == 'topo_node_match':
                self.topological_sub_resolved += 1
                problem_label = 'SortedViaTopoSubstitution'
            elif match_coming_from == 'no_topo_node_match':
                self.no_initial_topological_sub_resolved += 1
                problem_label = 'SolvedNoThenTopoSubEdit'
            idx_to_use = indexes_to_compare[number_of_changes_to_make.index(min_value_topo)]
        else:
            ### Examine each of the substructures as sometimes there is no difference
            substructures_to_examine = [indexes_to_compare[i] for i, v in enumerate(min_topo_list) if v == True]
            idx_to_use, problem_label = self.comparing_chemical_graph(substructures_to_examine,
                                                                      match_coming_from)
        return idx_to_use, problem_label

    def topological_distance_nodes(self):
        ### find where True
        true_list = []
        for value in self.all_lengths_between_nodes:
            if str(value) in self.core_length_mappings[self.row['core']].keys():
                true_list.append(self.core_length_mappings[self.row['core']][str(value)])
            else:
                true_list.append(0)
        ### find the max length
        max_value = max(true_list)
        max_count_list = [(t == max_value) for t in true_list]
        
        if max_count_list.count(True) == 1:
            problem_label = 'SortedViaTopoMap'
            self.topological_map_resolved += 1
            idx_to_use = true_list.index(max_value)
        else:
            ###get indexes of the ones to compare - only wnt to compare those top examples
            indexes_to_compare = [i for i, v in enumerate(max_count_list) if v == True]
            idx_to_use, problem_label = self.topological_substitutions(indexes_to_compare,
                                                                       'topo_node_match')
            
        return idx_to_use, problem_label
    
    def topological_edit_distance(self):
        ### Perform edit distance to current mappings
        all_number_of_changes_to_make = []
        number_of_changes_to_make = []
        for mapping in self.all_lengths_between_nodes:
            vals = list(mapping.values())
            all_number_of_changes = []
            lowest_number_of_changes = None
            for sub in self.core_length_mappings[self.row['core']]:
                sub = list(ast.literal_eval(sub).values())
                changes_to_make = 0
                for node_idx in range(0, len(sub)):
                    a=vals[node_idx]
                    b=sub[node_idx]
                    changes_to_make += abs(a-b)
                all_number_of_changes.append(changes_to_make)
                if lowest_number_of_changes == None:
                    lowest_number_of_changes = changes_to_make
                elif lowest_number_of_changes > changes_to_make:
                    lowest_number_of_changes = changes_to_make
            number_of_changes_to_make.append(lowest_number_of_changes)
            all_number_of_changes_to_make.append(all_number_of_changes)
        
        if number_of_changes_to_make.count(min(number_of_changes_to_make)) == 1:
            index_to_use = number_of_changes_to_make.index(min(number_of_changes_to_make))
            problem_label = 'SolvedNoThenTopoEditDistance'
            self.no_initial_edit_distance += 1
        else:
            ### look into the average edit distance
            indexes_to_compare = [i for i, v in enumerate(number_of_changes_to_make) if v == min(number_of_changes_to_make)]
            all_averages = [np.mean(v) for v in all_number_of_changes_to_make]
            averages = [np.mean(v) for i, v in enumerate(all_number_of_changes_to_make) if i in indexes_to_compare]
            if averages.count(min(averages)) == 1:
                index_to_use = all_averages.index(min(averages))
                problem_label = 'SolvedNoThenAvgTopEditDistance'
                self.no_initial_avg_edit_distance += 1
            else:
                ### select the indexes that are the minimum changes
                index_to_use, problem_label = self.topological_substitutions(indexes_to_compare,
                                                                             'no_topo_node_match')
        return index_to_use, problem_label
    
def no_single_mappings_function(no_single_mappings_df,
                                no_single_mappings,
                                params,
                                multi_mappings_df,
                                row_indexes_in_main_df):
    core_not_present_in_overlap = len(no_single_mappings_df)
    number_of_cores = len(set(no_single_mappings_df['core']))
    fixed_by_max = 0
    fixed_by_cg = 0
    fixed_by_hac = 0
    fixed_by_mw = 0
    fixed_by_max_core_seen = 0
    fixed_by_max_core_seen_arbit = 0
    fixed_by_arbit = 0
    
    for core in set(no_single_mappings_df['core']):
        filtered_df = no_single_mappings_df[no_single_mappings_df['core'] == core]
        node_topo_map = []
        for row_idx, row in filtered_df.iterrows():
            node_topo_map.append(no_single_mappings[row_idx][0])
        
        node_topo_map = [str(y) for x in node_topo_map for y in x]
        topo_vals = {tm:node_topo_map.count(tm) for tm in set(node_topo_map)}
        
        for row_idx, row in filtered_df.iterrows():
            smiles = Chem.MolFromSmiles(row['SMILES'])
            core = Chem.MolFromSmiles(row['core'], params)
        
            ### decompress the sdf
            mol_block = zlib.decompress(base64.b64decode(row['SDF']))
            molecule = Chem.MolFromMolBlock(mol_block, sanitize=False, strictParsing=False, removeHs=True)
        
            ### remove the hydrogens from each rg node and set meta data (indexes)
            for atomI in range(0, molecule.GetNumAtoms()):
                molecule.GetAtomWithIdx(atomI).SetNoImplicit(True)
                molecule.GetAtomWithIdx(atomI).SetProp('rg_idx', str(atomI))
            
            all_core_smiles = no_single_mappings[row_idx][1]
            all_core_info = no_single_mappings[row_idx][3]
            
            topos = no_single_mappings[row_idx][0]
            indexes = [topo_vals[str(k)] for k in topos]
            max_indexes = [i for i, v in enumerate(indexes) if v == max(indexes)]
            if len(max_indexes) == 1:
                index_to_use = max_indexes[0]
                problem_label = 'CoreNotPresentInSingleMaxOverlap'
                fixed_by_max+=1
            elif len(max_indexes) == 0:
                print('cmpare edit distance')
                ### use the one that has the lowest edit distance to the main one
                sys.exit()
            else:
                substructures_to_compare = []
                for sub in max_indexes:
                    core_smiles_atom_indexes = all_core_smiles[sub]
                    core_smiles = Chem.RWMol(smiles)
                    not_core_atoms = [i for i in range(0, (core_smiles.GetNumAtoms())) if i not in core_smiles_atom_indexes]
                    ### remove the atoms that aren't part of the core
                    for atom_idx in sorted(not_core_atoms, reverse=True):
                        core_smiles.RemoveAtom(atom_idx)
                    substructures_to_compare.append(Chem.MolToSmiles(core_smiles, 1))
                if len(set(substructures_to_compare)) == 1:
                    index_to_use = max_indexes[0]
                    problem_label = 'CoreNotPresentInSingleCG'
                    fixed_by_cg+=1
                else:
                    ##HAC and MW
                    hac_values = [Chem.MolFromSmiles(substructure, params).GetNumAtoms() for substructure in substructures_to_compare]
                    max_hac_list = [(t == max(hac_values)) for t in hac_values]
                    if max_hac_list.count(max(max_hac_list)) == 1:
                        index_to_use = max_indexes[hac_values.index(max(hac_values))]
                        problem_label = 'CoreNotPresentInSingleHAC'
                        fixed_by_hac+=1
                    else:
                        substructures_to_compare = [substructures_to_compare[i] for i, v in enumerate(hac_values) if v == max(hac_values)]
                        substructures_to_examine = [max_indexes[i] for i, v in enumerate(hac_values) if v == max(hac_values)]
                        
                        mw_values = [Descriptors.ExactMolWt(Chem.MolFromSmiles(substructure, params)) for substructure in substructures_to_compare]
                        max_mw_list = [(t == max(mw_values)) for t in mw_values]
                        if max_mw_list.count(True) == 1:
                            fixed_by_mw+=1
                            problem_label = 'CoreNotPresentInSingleMW'
                            index_to_use = substructures_to_examine[mw_values.index(max(mw_values))]
                        else:
                            ### look at substructure
                            substructures_to_compare = [substructures_to_compare[i] for i, v in enumerate(max_mw_list) if v == True]
                            substructures_to_examine = [substructures_to_examine[i] for i, v in enumerate(max_mw_list) if v == True]
            
                            set_counts = {structure:substructures_to_compare.count(structure) for structure in set(substructures_to_compare)}
                            if list(set_counts.values()).count(max(set_counts.values())) == 1:
                                fixed_by_max_core_seen+=1
                                structure_to_use = max(set_counts, key= lambda x: set_counts[x])
                                index_to_use = substructures_to_examine[substructures_to_compare.index(structure_to_use)]
                                problem_label = 'CoreNotPresentInSingleOneCore/MostSeenCore'
                            elif list(set_counts.values()).count(max(set_counts.values())) > 1:
                                fixed_by_max_core_seen_arbit+=1
                                structure_to_use = max(set_counts, key= lambda x: set_counts[x])
                                index_to_use = substructures_to_examine[substructures_to_compare.index(structure_to_use)]
                                problem_label = 'CoreNotPresentInSingleOneCore/MostSeenCoreArbitrary'
                            else:
                                fixed_by_arbit+=1
                                index_to_use = substructures_to_compare[0]
                                problem_label = 'CoreNotPresentInSingleArbitrary'

            core_to_use = []
            core_to_use_smarts = []
            core_to_use_add_info = []
            
            core_smiles_atom_indexes = all_core_smiles[index_to_use]
            core_smiles = Chem.RWMol(smiles)
            not_core_atoms = [i for i in range(0, (core_smiles.GetNumAtoms())) if i not in core_smiles_atom_indexes]
            ### remove the atoms that aren't part of the core
            for atom_idx in sorted(not_core_atoms, reverse=True):
                core_smiles.RemoveAtom(atom_idx)
            core_smarts, additionals_list = extracting_more_meta_data(smiles,
                                                                      molecule,
                                                                      core_smiles,
                                                                      core,
                                                                      row)
            core_to_use = core_smiles_atom_indexes
            core_to_use_smarts = canonicalise_smarts(Chem.MolToSmiles(core_smarts))
            core_to_use_add_info = additionals_list
    
            exist_row_idx = row_indexes_in_main_df[row_idx]
            multi_mappings_df.at[exist_row_idx, 'ProblemMol'] = problem_label
            multi_mappings_df.at[exist_row_idx, 'core_node_indexes'] = row['node_idx_of_core_potential_mappings'][index_to_use]
            multi_mappings_df.at[exist_row_idx, 'core_atom_indexes'] = core_to_use
            multi_mappings_df.at[exist_row_idx, 'core_smarts'] = core_to_use_smarts
            multi_mappings_df.at[exist_row_idx, 'core_info'] = all_core_info[index_to_use]
            multi_mappings_df.at[exist_row_idx, 'additional'] = core_to_use_add_info
            
    print('Number of Molecules with no single mappings: ', core_not_present_in_overlap, ' across ', number_of_cores, ' cores')
    print('Fixed by selecting max overlap: ', fixed_by_max)
    print('Fixed by CG: ', fixed_by_cg)
    print('Fixed by HAC: ', fixed_by_hac)
    print('Fixed by MW: ', fixed_by_mw)
    print('Fixed by max core seen: ', fixed_by_max_core_seen)
    print('Fixed by max core seen arbitrary: ', fixed_by_max_core_seen_arbit)
    print('Fixed by Arbitrary: ', fixed_by_arbit)
    return multi_mappings_df

def processing_multiple_mappings(multi_mappings_df,
                                 params,
                                 core_length_mappings,
                                 core_substitution_info,
                                 core_smarts_mappings):
    zero_count = 0
    multi_count = 0
    core_not_present_in_1 = 0
    single_topo_match = 0
    
    ###set the class
    emRGc = extracting_mapped_rg_core()
    ###set the info for all
    emRGc.set_value_of_single_mappings(multi_mappings_df,
                                       core_length_mappings,
                                       core_substitution_info,
                                       core_smarts_mappings)
    
    no_single_mappings_df = pd.DataFrame()
    no_single_mappings = []
    no_single_mappings_row_idxs = []
    
    print('Number of Molecules That Have Multiple Mappings: ', len(multi_mappings_df))
    for row_idx, row in multi_mappings_df.iterrows():
        smiles = Chem.MolFromSmiles(row['SMILES'])
        core = Chem.MolFromSmiles(row['core'], params)
    
        ### decompress the sdf
        mol_block = zlib.decompress(base64.b64decode(row['SDF']))
        molecule = Chem.MolFromMolBlock(mol_block, sanitize=False, strictParsing=False, removeHs=True)
    
        ### remove the hydrogens from each rg node and set meta data (indexes)
        for atomI in range(0, molecule.GetNumAtoms()):
            molecule.GetAtomWithIdx(atomI).SetNoImplicit(True)
            molecule.GetAtomWithIdx(atomI).SetProp('rg_idx', str(atomI))
        
        all_lengths_between_nodes = []
        all_core_smiles = []
        all_core_info = []
        all_sub_topo_distances = []
        for mcs_smarts in row['node_idx_of_core_potential_mappings']:
            length_between_nodes, core_smiles_atom_indexes, core_info, all_substitution_topological_distances = extracting_core_meta_data_and_lengths(smiles,
                                                                                                                                                      molecule,
                                                                                                                                                      mcs_smarts)
            if length_between_nodes != None:
                all_lengths_between_nodes.append(length_between_nodes)
            elif length_between_nodes == None:
                all_lengths_between_nodes.append(0)
            all_core_smiles.append(core_smiles_atom_indexes)
            all_core_info.append(core_info)
            all_sub_topo_distances.append(all_substitution_topological_distances)
        emRGc.set_values(row,
                         smiles,
                         molecule,
                         core,
                         all_lengths_between_nodes,
                         all_core_smiles,
                         all_sub_topo_distances)    
            
        if row['core'] in core_length_mappings.keys():
            list_of_t_f = [(str(value) in core_length_mappings[row['core']].keys()) for value in all_lengths_between_nodes]
        else:
            core_not_present_in_1 += 1
            list_of_t_f = []
        if list_of_t_f.count(True) == 1:
            single_topo_match +=1
            index_to_use = list_of_t_f.index(True)
            problem_label = ''
            
        elif len(list_of_t_f) == 0:
            index_to_use = 0
            problem_label = 'CoreNotPresentInSingle'
            no_single_mappings_df = no_single_mappings_df.append(row, ignore_index=True)
            no_single_mappings.append([all_lengths_between_nodes,
                                       all_core_smiles,
                                       all_sub_topo_distances,
                                       all_core_info])
            no_single_mappings_row_idxs.append(row_idx)
            
        else:                        
            if list_of_t_f.count(True) == 0:  
                zero_count += 1
                index_to_use, problem_label = emRGc.topological_edit_distance() 
            else:
                multi_count += 1
                index_to_use, problem_label = emRGc.topological_distance_nodes()
                                        
        core_to_use = []
        core_to_use_smarts = []
        core_to_use_add_info = []
        core_smiles_atom_indexes = all_core_smiles[index_to_use]
        core_smiles = Chem.RWMol(smiles)
        not_core_atoms = [i for i in range(0, (core_smiles.GetNumAtoms())) if i not in core_smiles_atom_indexes]
        ### remove the atoms that aren't part of the core
        for atom_idx in sorted(not_core_atoms, reverse=True):
            core_smiles.RemoveAtom(atom_idx)
        core_smarts, additionals_list = extracting_more_meta_data(smiles,
                                                                  molecule,
                                                                  core_smiles,
                                                                  core,
                                                                  row)
        core_to_use = core_smiles_atom_indexes
        core_to_use_smarts = canonicalise_smarts(Chem.MolToSmiles(core_smarts))
        core_to_use_add_info = additionals_list

        multi_mappings_df.at[row_idx, 'ProblemMol'] = problem_label
        multi_mappings_df.at[row_idx, 'core_node_indexes'] = row['node_idx_of_core_potential_mappings'][index_to_use]
        multi_mappings_df.at[row_idx, 'core_atom_indexes'] = core_to_use
        multi_mappings_df.at[row_idx, 'core_smarts'] = core_to_use_smarts
        multi_mappings_df.at[row_idx, 'core_info'] = all_core_info[index_to_use]
        multi_mappings_df.at[row_idx, 'additional'] = core_to_use_add_info
        multi_mappings_df.at[row_idx, 'length_between_core_nodes'] = all_lengths_between_nodes[index_to_use]

    if len(no_single_mappings_df) > 0:
        multi_mappings_df = no_single_mappings_function(no_single_mappings_df,
                                                        no_single_mappings,
                                                        params,
                                                        multi_mappings_df,
                                                        no_single_mappings_row_idxs)
    print(zero_count, multi_count, core_not_present_in_1)
    print('Single Match to Topological Mapping: ', single_topo_match)
    emRGc.print_stats()
    
    return multi_mappings_df

def establishing_core(initial_df):
    params = Chem.SmilesParserParams()
    params.removeHs = True
    params.sanitize = False
    core_numbered, core_edges = core_edges_and_isotopes(set(initial_df['core'].tolist()),
                                                        params)

    ### add the new columns 
    initial_df['number_of_core_potential_mappings'] = ""
    initial_df['core_numbered'] = ""
    initial_df['core_node_indexes'] = ""
    initial_df['core_atom_indexes'] = ""
    initial_df['core_edges'] = ""
    initial_df['core_smarts'] = ""
    initial_df['core_info'] = ""
    initial_df['additional'] = ""
    initial_df['length_between_core_nodes'] = ""
    initial_df['node_idx_of_core_potential_mappings'] = ""
    
    for idx, row in initial_df.iterrows():
        initial_df.at[idx, 'core_numbered'] = core_numbered[row['core']]
        initial_df.at[idx, 'core_edges'] = core_edges[row['core']]
        core = Chem.MolFromSmiles(row['core'], params)
        mol_block = zlib.decompress(base64.b64decode(row['SDF']))
        molecule = Chem.MolFromMolBlock(mol_block, sanitize=False, strictParsing=False, removeHs=True)
        mcs_smarts = molecule.GetSubstructMatches(core)
        initial_df.at[idx, 'number_of_core_potential_mappings'] = len(mcs_smarts)
        initial_df.at[idx, 'node_idx_of_core_potential_mappings'] = list(mcs_smarts)
    
    one_mappings_df = initial_df[initial_df['number_of_core_potential_mappings'] == 1]
    core_length_mappings, core_smarts_mappings, core_substitution_info = process_single_mappings(one_mappings_df,
                                                                                                 params)

    multi_mappings_df = deepcopy(initial_df[initial_df['number_of_core_potential_mappings'] != 1])
    multi_mappings_df['ProblemMol'] = None
    
    multi_mappings_df = processing_multiple_mappings(multi_mappings_df,
                                                     params,
                                                     core_length_mappings,
                                                     core_substitution_info,
                                                     core_smarts_mappings)
    
    combined_df = one_mappings_df.append(multi_mappings_df, ignore_index=True, sort=False)

    return combined_df
        
if __name__ == '__main__':
    args = input_args()
    mmpm = molecules_with_multiple_possible_mappings()
    
    initial_df = read_in_files()

#    ### establish the cores 
    complete_df = establishing_core(initial_df)
    
    complete_df.to_csv(args.output, sep=args.separator, index=False)