##############################################################################
# Reduced graph program
#
#
#
# Jess Stacey - modified Eleanor Gardiner C++ code
##############################################################################

from itertools import chain
import copy
import os
import time
import rdkit
from rdkit import Chem

import canonicalisor
from setting_nodes import defining_neighbours, set_nodes, set_basic_node, set_unassigned_node, reassigning_node, unknown_set_node, set_combined_node, set_halogen_node
from setting_ring_nodes import establishing_ring_nodes
from reduced_graph_arguments import input_args
from connecting_reduced_graph import generating_rg_from_scratch 

max_nodes = 200

class ReadingFiles:
    ''' Class: ReadingFiles
    This class does error handling of the smarts by reading in either smarts
    file and canonicalises them before saving them in a dictionary with there
    respective node type
    '''
    def __init__(self):
        self.smarts_dictionary = {}
        self.smiles_dictionary = {}
        self.smiles_unable_to_read = []
        self.duplicate_smiles = []

    def smarts_file(self):
        ''' Function: smarts_file
        This function reads the smarts file and canonicalises them and stores
        them in the smarts a {}
        input: self
        output: smarts_dictionary - {} of all the different node types smarts
        '''
        try:
            if args.smarts == 'smarts.smt':
                smart_file = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), args.smarts), 'r')
            else:
                smart_file = open(args.smarts, 'r')
        except:
            print(args.smarts)
            raise ValueError('Cannot open SMARTS file')
        for line in smart_file:
            smarts_line = line.replace("\n", "").replace("\r", "").split('\t')
            can_smarts, mol = canonicalisor.canonicalise_smarts(smarts_line[0], 0)
            self.smarts_dictionary[can_smarts] = {}
            self.smarts_dictionary[can_smarts]['Type'] = smarts_line[1]
            self.smarts_dictionary[can_smarts]['mol'] = mol

        return self.smarts_dictionary

    def smarts_isolating_file(self):
        ''' Function: smarts_isolating_file
        This function reads the isolating smarts file and canonicalises them
        and stores them in the smarts a {}
        input: self
        output: smarts_dictionary - {} of all the different node types smarts
        '''
        try:
            if args.carbon == 'isolating_carbon.smt':
                isolating_carbon_file = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), args.carbon), 'r')
            else:
                isolating_carbon_file = open(args.carbon, 'r')
        except:
            raise ValueError('Cannot open isolating carbons file')
        for line in isolating_carbon_file:
            can_smarts, mol = canonicalisor.canonicalise_smarts(line.replace("\n", ""), 1)
            self.smarts_dictionary[can_smarts] = {}
            self.smarts_dictionary[can_smarts]['Type'] = 'isolating carbon'
            self.smarts_dictionary[can_smarts]['mol'] = mol

        return self.smarts_dictionary

class SubStructureSearching:
    ''' Class: SubStructureSearching
    This class does a substructure search on the respective smarts type to
    see whether there are any present within the smiles_mol object, all
    return a list of the indexes that match the substructure search
    '''
    def __init__(self, can_smiles, smiles_mol):
        number_atoms = smiles_mol.GetNumAtoms()
        self.assigned_atoms = [-1]*number_atoms
        self.metal_atoms = [-1]*number_atoms
        self.halogens_atoms = [-1]*number_atoms
        self.linker_atoms = [-1]*number_atoms
        self.dictionary_of_nodes = {}

        if args.verbose:
            print(can_smiles + ' contains ' + str(number_atoms) + ' (non-H) atoms')
            list_atomic_numbers(smiles_mol)

    def find_metal_atoms(self,
                         smiles_mol):
        ''' Function: find_metal_atoms
        This function substructure searches for the metal atoms that are within
        the metal list that gets imported
        input: smiles_mol - mol object of the smiles that is being searched
        ouput: matches - [] of tuples of the atom indexes of the atoms that match
                        the metals in the metal list
        '''
        from metals import metals
        metal_matches = []
        metal_smarts = ','.join(metals)
        metal_smarts = '[' + metal_smarts + ']'
        metals_mol = Chem.MolFromSmarts(metal_smarts)

        if smiles_mol.HasSubstructMatch(metals_mol, useChirality=True):
            metal_matches = list(smiles_mol.GetSubstructMatches(metals_mol,
                                                                useChirality=True))
            for metal in metal_matches:
                for atom in metal:
                    self.metal_atoms[atom] = 1
                    self.assigned_atoms[atom] = 1
                if args.verbose:
                    print('Number of Matches ' + str(len(metal_matches)))
                    if len(metal_matches) >= 1:
                        print('Indexes that match: ',
                              [index for tup in metal_matches for index in tup])
        else:
            if args.verbose:
                print('No metal matches')

        return metal_matches

    def find_atoms(self,
                   smarts_dicitionary,
                   smiles_mol,
                   node_type):
        ''' Function: find_atoms
        This function substructure searches for the atoms that are defined in
        the smarts dictionary that correspond to the node type
        input: smarts_dicitionary - {} of all the different node types smarts
               smiles_mol - mol object of the smiles that is being searched
               node_type - string of the type of node that is being searched for
        ouput: matches - [] of tuples of the atom indexes of the atoms that match the node type
        '''
        node_matches = []
        for value in smarts_dicitionary.values():
            if value['Type'] == node_type:
                smarts_mol = value['mol']
                if smiles_mol.HasSubstructMatch(smarts_mol, useChirality=True):
                    node_sub_matches = list(smiles_mol.GetSubstructMatches(smarts_mol,
                                                                           useChirality=True))
                    node_matches.extend(node_sub_matches)
                    for match in node_sub_matches:
                        for atom in match:
                            if node_type != 'isolating carbon':
                                self.assigned_atoms[atom] = 1     
        if args.verbose:
            print('Number of Matches ' + str(len(set(chain.from_iterable(node_matches)))))
            if len(node_matches) >= 1:
                print('Indexes that match: ',
                      [index for tup in node_matches for index in tup])
            else:
                print(' '.join(['No', node_type, 'matches']))

        return list(set(node_matches))

    def assign_ring_nodes(self,
                          list_of_rings):
        ''' Function: assign_ring_nodes
        This function makes the atom indexes of the atoms that are contain within a ring be
        classed as assigned (1) in the list assigned_atoms
        input: list_of_rings - [] of tuples of the atom indexes that are within rings
        output: None
        '''
        setting_ring_nodes = set(list(chain.from_iterable(list_of_rings)))
        for atom in setting_ring_nodes:
            self.assigned_atoms[atom] = 1

        return None

def find_halogens(smiles_mol):
    ''' Function: find_halogens
    This function substructure searches for the atoms that are halogens and
    returns a list of the atom indexes that match
    input: smiles_mol - mol object of the smiles that is being searched
    output: matches_list - [] of tuples of atom indexes that are halogens
    '''
    matches_list = []
    matches_list.extend(list(smiles_mol.GetSubstructMatches(Chem.MolFromSmarts('F'))))
    matches_list.extend(list(smiles_mol.GetSubstructMatches(Chem.MolFromSmarts('Cl'))))
    matches_list.extend(list(smiles_mol.GetSubstructMatches(Chem.MolFromSmarts('Br'))))
    matches_list.extend(list(smiles_mol.GetSubstructMatches(Chem.MolFromSmarts('I'))))

    return matches_list

def linker_function(atom_index,
                    list_of_neighbours,
                    neighbour_atom_index,
                    potential_iso_dict,
                    neighbours,
                    potential_iso_dict_keys_to_delete):
    ''' Function: linker_function
    This function looks for linker nodes
    input: atom_index - integer of the atom index
           list_of_neighbours - [] of the neighbours
           neighbour_atom_index - integer of the atom index of a neighbour
           potential_iso_dict - {} of the atom indexes (keys) and their neighbours (values)
           neighbours - {} of the atom indexes (keys) and their neighbours (values)
           potential_iso_dict_keys_to_delete - [] of the dictionary keys that will be deleted from potential_iso_dict
    output: None
    '''
    if [item for item in set(potential_iso_dict.keys()) if atom_index not in item] != []:
        if [item for item in list(potential_iso_dict.values()) if atom_index in item] != []:
            neighbours_index_combine = [item for item in list(potential_iso_dict.values()) if atom_index in item]
            key_list_index_combine = [atom_index]
            list_of_keys = [list(potential_iso_dict.keys())[list(potential_iso_dict.values()).index(neighbours_atom_index)] for neighbours_atom_index in neighbours_index_combine]
            key_list_index_combine.extend(list(chain.from_iterable(list_of_keys)))
            neighbours_index_combine = list(chain.from_iterable(neighbours_index_combine))
            if tuple(set(key_list_index_combine)) in list(potential_iso_dict.keys()):
                potential_iso_dict_keys_to_delete.append(tuple(set(key_list_index_combine)))

            neighbours_index_combine.extend(neighbours[atom_index])
            neighbours_index_combine = list(set(neighbours_index_combine))

            neighbours_index_combine = [key_combo for key_combo in neighbours_index_combine if key_combo not in key_list_index_combine]
            for key_combo in set(key_list_index_combine):
                if type(key_combo) == int:
                    key_combo = (key_combo,)
                elif type(key_combo) == list:
                    key_combo = tuple(key_combo)
                potential_iso_dict_keys_to_delete.append(key_combo)

            new_key = tuple(set(key_list_index_combine))
            new_key = sorted(new_key)
            if tuple(new_key) not in potential_iso_dict.keys():
                potential_iso_dict[tuple(new_key)] = neighbours_index_combine

        else:
            key_to_combine = list(neighbours.keys())[list(neighbours.values()).index(neighbour_atom_index)]
            list_index = [atom_index,
                          key_to_combine]
            neighbours_index = [list_of_neighbours,
                                neighbours[key_to_combine]]
            neighbours_index = list(chain.from_iterable(neighbours_index))
            neighbours_index = [key_combo for key_combo in neighbours_index if key_combo not in list_index]
            new_key = tuple(set(list_index))
            new_key = sorted(new_key)
            if tuple(new_key) not in potential_iso_dict.keys():
                potential_iso_dict[tuple(new_key)] = neighbours_index

    return None

def potential_hydrophobic_node(sss,
                               smiles_mol,
                               carbon):
    ''' Function: potential_hydrophobix_node
    This function searches to see whether the node that has been passed to this
    function should be classed as inert or hydrophobic based on the atoms,
    hybridization and valence of the atoms within this node
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           carbon - () of atom indexes
    output: True (bool)
    '''
    atoms_list = [smiles_mol.GetAtomWithIdx(atom_index).GetSmarts() for atom_index in carbon]
    hybridization_list = [smiles_mol.GetAtomWithIdx(atom_index).GetHybridization() for atom_index in carbon]
    valence_list = [smiles_mol.GetAtomWithIdx(atom_index).GetImplicitValence() for atom_index in carbon]

    if (set(atoms_list) == {'C'} and
            [set(hybridization_list), set(valence_list)] == [{Chem.rdchem.HybridizationType.SP3}, {2, 3}] or
            [set(hybridization_list), set(valence_list)] == [{Chem.rdchem.HybridizationType.SP3}, {3}]):
        set_basic_node(sss,
                       smiles_mol,
                       carbon,
                       'inert')
    else:
        set_basic_node(sss,
                       smiles_mol,
                       carbon,
                       'hydrophobic')

    return True

def differentiate_linker_isolating_carbon(sss,
                                          smiles_mol,
                                          potential_iso_dict):
    ''' Function: differentiate_linker_isolating_carbon
    This function differentiates between linker carbons and terminal carbons
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           potential_iso_dict - {} of the atom indexes (keys) and their neighbours (values)
    output: None
    '''
    for carbon, neighbours in potential_iso_dict.items():
        set_node = False
        linker_nodes = [key for key, value in sss.dictionary_of_nodes.items() if value['node_type'] == 'acyclic inert']
        if len(neighbours) == 1:
            if sss.halogens_atoms[neighbours[0]] == 1:
                node_tuple = [node for node, node_value in sss.dictionary_of_nodes.items() if neighbours[0] in node]
                for atom_indexes in node_tuple:
                    del sss.dictionary_of_nodes[atom_indexes]
                    carbon = carbon + atom_indexes

            if args.hydrophobic:
                set_node = potential_hydrophobic_node(sss,
                                                      smiles_mol,
                                                      carbon)
            elif args.terminallinker:
                set_node = True
                set_basic_node(sss,
                               smiles_mol,
                               carbon,
                               'inert')
            else:
                if sss.halogens_atoms[neighbours[0]] == 1:
                    set_node = True
                    set_basic_node(sss,
                                   smiles_mol,
                                   carbon,
                                   'inert')
        elif len(neighbours) > 1:
            for neighbour_index in neighbours:
                for index in [key for key in linker_nodes if neighbour_index in key]:
                    del sss.dictionary_of_nodes[index]
                    carbon = carbon + index
            set_node = True
            set_basic_node(sss,
                           smiles_mol,
                           carbon,
                           'inert')
        else:
            set_node = True
            set_basic_node(sss,
                           smiles_mol,
                           carbon,
                           'inert')
        if set_node is True:
            for carbon_index in carbon:
                sss.assigned_atoms[carbon_index] = 1
                sss.linker_atoms[carbon_index] = 1

    return None

def establishing_isolating_c_nodes(sss,
                                   smiles_mol,
                                   smarts_dictionary):
    ''' Function: establishing_isolating_c_nodes
    This function finds the isolating carbon nodes and then sees whether they should be a linker or not
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           smarts_dictionary - {} of all the different node types smarts
    output: None
    '''
    if args.verbose:
        print('Establishing Isolating Carbons')

    matches = sss.find_atoms(smarts_dictionary,
                             smiles_mol,
                             'isolating carbon')

    matches = list(chain.from_iterable(matches))
    potential_iso_c = [atom_number for atom_number in range(0, smiles_mol.GetNumAtoms()) if (sss.assigned_atoms[atom_number] == -1 and atom_number in matches)]
    neighbours_dictionary = {}
    potential_iso_dict = {}
    potential_iso_dict_keys_to_delete = []

    for index in potential_iso_c:
        neighbours_dictionary[index] = [bond.GetOtherAtomIdx(index) for bond in smiles_mol.GetAtomWithIdx(index).GetBonds()]

    for key, value in neighbours_dictionary.items():
        for atom_index in neighbours_dictionary.values():
            if key in atom_index:
                linker_function(key,
                                value,
                                atom_index,
                                potential_iso_dict,
                                neighbours_dictionary,
                                potential_iso_dict_keys_to_delete)
            else:
                potential_iso_dict[(key,)] = value

    pid_keys = list(set(potential_iso_dict.keys()))
    ### This adds the keys that are subsets of larger keys to the list
    potential_iso_dict_keys_to_delete.extend([node_key for index, node_key in enumerate(pid_keys) if any(set(node_key).issubset(set(smaller_node_key)) for smaller_node_key in (pid_keys[:index] + pid_keys[index+1:]))])

    for dict_key in set(potential_iso_dict_keys_to_delete):
        if dict_key in set(potential_iso_dict.keys()):
            del potential_iso_dict[dict_key]

    differentiate_linker_isolating_carbon(sss,
                                          smiles_mol,
                                          potential_iso_dict)

    if args.verbose:
        for atom_number in range(0, smiles_mol.GetNumAtoms()):
            if atom_number in matches:
                print('Atom ' + smiles_mol.GetAtomWithIdx(atom_number).GetSmarts() + ' is isolating Carbon')
            elif atom_number in matches:
                print('Atom ' + smiles_mol.GetAtomWithIdx(atom_number).GetSmarts() + ' is not isolating Carbon')

    return None

def setting_unassigned(sss,
                       smiles_mol,
                       unassigned_neighbours):
    ''' Function: setting_unassigned
    This function sets the unassigned atoms to existing nodes or creates new
    nodes depending on the arguments
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           unassigned_neighbours - {} of the atoms that are unassigned, key,
                                   and the value is their neighbours
    output: None
    '''
    rdkit_version = float('.'.join([rdkit.__version__.split('.')[0], rdkit.__version__.split('.')[1]]))
    for neighbour_index in unassigned_neighbours.values():
        acyclic_hba_hbd = [key for key, value in sss.dictionary_of_nodes.items() if 'acyclic' in value['node_type'] and ('HBA' in value['node_type'] or 'HBD' in value['node_type'])]
        linker_nodes = [key for key, value in sss.dictionary_of_nodes.items() if 'acyclic inert' in value['node_type']]
        set_node = False
        key = list(unassigned_neighbours.keys())[list(unassigned_neighbours.values()).index(neighbour_index)]
        
        if rdkit_version >= 2018.09:
            if type(key) == int:
                key_atoms = [smiles_mol.GetAtomWithIdx(key).GetSymbol()]
                key = (key,)
            elif type(key) == tuple:
                key_atoms = [smiles_mol.GetAtomWithIdx(atom_index).GetSymbol() for atom_index in key]
        else:
            if type(key) == int:
                key_atoms = [smiles_mol.GetAtomWithIdx(key).GetSmarts()]
                key = (key,)
            elif type(key) == tuple:
                key_atoms = [smiles_mol.GetAtomWithIdx(atom_index).GetSmarts() for atom_index in key]

        if set(key_atoms) != {'C'}:
            hba_hbd_all_matches = [item for neighbour in neighbour_index for item in acyclic_hba_hbd if neighbour in item]
            if hba_hbd_all_matches != []:
                new_key = [key]
                new_node_type = [sss.dictionary_of_nodes[match]['node_type'] for match in set(hba_hbd_all_matches)]
                new_key.extend([match for match in set(hba_hbd_all_matches)])

                for match in set(hba_hbd_all_matches):
                    del sss.dictionary_of_nodes[match]

                new_key = list(set(chain.from_iterable(new_key)))
                new_key.sort()
                new_key = tuple(new_key)

                if 'acyclic HBD HBA' in new_node_type:
                    node_type = 'HBD HBA'
                elif 'acyclic HBD' in new_node_type and 'acyclic HBA' in new_node_type:
                    node_type = 'HBD HBA'
                elif 'acyclic HBD' in new_node_type:
                    node_type = 'HBD'
                elif 'acyclic HBA' in new_node_type:
                    node_type = 'HBA'
                for atom_number in key:
                    sss.assigned_atoms[atom_number] = 1
                set_nodes(sss,
                          smiles_mol,
                          [new_key],
                          node_type)
                set_node = True

        if set_node is False:
            set_unassigned_node(args,
                                sss,
                                smiles_mol,
                                linker_nodes,
                                neighbour_index,
                                key)

    return None

def unassigned_charged(sss,
                       smiles_mol,
                       atom_index,
                       unassigned_neighbours,
                       items_to_delete):
    ''' Function: unassigned_charged
    This function assigns the unassigned_charged atoms in smiles_mol
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           atom_index - integer of the atom atom_index that is an unassigned
           unassigned_neighbours - {} of the unassigned atoms as key and their neighbours as values
           items_to_delete - [] list of the atoms index that need to be removed from unassigned_neighbours
    output: additional_other_atoms_ids - [] of the atom indexes that are neighbours
    '''
    initial_index = copy.deepcopy(atom_index)
    sss.assigned_atoms[initial_index] = 1
    atom = smiles_mol.GetAtomWithIdx(atom_index)
    other_atoms_ids = [bond.GetOtherAtomIdx(atom_index) for bond in atom.GetBonds()]
    unassigned_neighbours[atom_index] = other_atoms_ids
    items_to_delete.append(atom_index)
    if len(other_atoms_ids) == 0:
        set_basic_node(sss,
                       smiles_mol,
                       (atom_index,),
                       'inert')

        return other_atoms_ids

    else:
        additional_other_atoms_ids = copy.deepcopy(other_atoms_ids)
        delete_nodes = []
        linker_nodes = [key for other_atom in other_atoms_ids if smiles_mol.GetAtomWithIdx(other_atom).IsInRing() is False for key, value in sss.dictionary_of_nodes.items() if other_atom in key and (value['node_type'] == 'acyclic inert')]
        hba_hbd_nodes = [key for other_atom in other_atoms_ids if smiles_mol.GetAtomWithIdx(other_atom).IsInRing() is False for key, value in sss.dictionary_of_nodes.items() if other_atom in key and 'acyclic' in value['node_type'] and ('HBA' in value['node_type'] or 'HBD' in value['node_type'])]
        atom_index = (atom_index,)
        if hba_hbd_nodes != []:
            for node in hba_hbd_nodes:
                if tuple(set(node + atom_index)) != node:
                    delete_nodes.append(node)
                    atom_index = reassigning_node(args,
                                                  sss,
                                                  smiles_mol,
                                                  node,
                                                  atom_index)
        elif linker_nodes != []:
            for node in linker_nodes:
                atom_index += node
                atom_index = tuple(set(atom_index))
                if node != atom_index:
                    delete_nodes.append(node)
            set_basic_node(sss,
                           smiles_mol,
                           atom_index,
                           'inert')
        else:
            unknown_set_node(args,
                             sss,
                             smiles_mol,
                             atom_index,
                             other_atoms_ids)

        for delete in set(delete_nodes):
            if delete in set(sss.dictionary_of_nodes.keys()):
                del sss.dictionary_of_nodes[delete]

        additional_other_atoms_ids = [atom for atom in additional_other_atoms_ids if atom not in atom_index]

        return additional_other_atoms_ids

def assigning_unassigned(sss,
                         smiles_mol):
    ''' Function: assigning_unassigned
    This function assigns the current unassigned atoms within smiles_mol
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
    output: None
    '''
    unassigned_atoms = [atom_index for atom_index, assigned in enumerate(sss.assigned_atoms) if assigned == -1]
    unassigned_neighbours = {}
    items_to_delete = []
    linker_nodes = [key for key, value in sss.dictionary_of_nodes.items() if value['node_type'] == 'acyclic inert']
    for index in unassigned_atoms:
        if smiles_mol.GetAtomWithIdx(index).GetFormalCharge() != 0:
            additional_other_atoms_ids = unassigned_charged(sss,
                                                            smiles_mol,
                                                            index,
                                                            unassigned_neighbours,
                                                            items_to_delete)
        else:
            other_atoms_ids = [bond.GetOtherAtomIdx(index) for bond in smiles_mol.GetAtomWithIdx(index).GetBonds()]
            additional_other_atoms_ids = [bond.GetOtherAtomIdx(index) for bond in smiles_mol.GetAtomWithIdx(index).GetBonds()]
            unassigned_neighbours[index] = other_atoms_ids

            for other_atom in other_atoms_ids:
                if (smiles_mol.GetAtomWithIdx(other_atom).IsInRing() is False and
                        smiles_mol.GetBondBetweenAtoms(other_atom, index).GetBondType() == Chem.rdchem.BondType.DOUBLE and
                            set([smiles_mol.GetAtomWithIdx(other_atom).GetSmarts(), smiles_mol.GetAtomWithIdx(index).GetSmarts()]) != {'C'}):
                    if other_atom not in linker_nodes and other_atom in additional_other_atoms_ids:
                        additional_other_atoms_ids = unassigned_charged(sss,
                                                                        smiles_mol,
                                                                        index,
                                                                        unassigned_neighbours,
                                                                        items_to_delete)

    for item in items_to_delete:
        if item in set(unassigned_neighbours.keys()):
            del unassigned_neighbours[item]
    if len(unassigned_neighbours) == 1:
        setting_unassigned(sss,
                           smiles_mol,
                           unassigned_neighbours)
    elif len(unassigned_neighbours) > 1:
        unassigned_neighbours2 = {}
        for unassigned_index in set(unassigned_neighbours.keys()):
            for neighbour in unassigned_neighbours.values():
                if unassigned_index in neighbour:
                    other_key = list(unassigned_neighbours.keys())[list(unassigned_neighbours.values()).index(neighbour)]
                    if ((unassigned_index, other_key) in unassigned_neighbours2.keys()) is False and ((other_key, unassigned_index) in unassigned_neighbours2.keys()) is False:
                        if [key for key, value in unassigned_neighbours2.items() if unassigned_index in key] != [] or [key for key, value in unassigned_neighbours2.items() if other_key in key] != []:
                            
                            if unassigned_index in list(chain.from_iterable(unassigned_neighbours2.keys())) and other_key in list(chain.from_iterable(unassigned_neighbours2.keys())):
                                old_key1 = [key for key, value in unassigned_neighbours2.items() if unassigned_index in key][0]
                                old_key2 = [key for key, value in unassigned_neighbours2.items() if other_key in key][0]
                                new_neighbours = [unassigned_neighbours[unassigned_index],
                                              unassigned_neighbours[other_key],
                                              unassigned_neighbours2[old_key1],
                                              unassigned_neighbours2[old_key2]]
                                old_key = tuple(set(old_key1 + old_key2))
                                if old_key1 != old_key2:
                                    del unassigned_neighbours2[old_key1]
                                    del unassigned_neighbours2[old_key2]
                                else:
                                    del unassigned_neighbours2[old_key1]
                                
                            
                            if [key for key, value in unassigned_neighbours2.items() if unassigned_index in key] != []:
                                old_key = [key for key, value in unassigned_neighbours2.items() if unassigned_index in key][0]
                                new_neighbours = [unassigned_neighbours[unassigned_index],
                                              unassigned_neighbours[other_key],
                                              unassigned_neighbours2[old_key]]
                                del unassigned_neighbours2[old_key]
                            elif [key for key, value in unassigned_neighbours2.items() if other_key in key] != []:
                                old_key = [key for key, value in unassigned_neighbours2.items() if other_key in key][0]
                                new_neighbours = [unassigned_neighbours[unassigned_index],
                                              unassigned_neighbours[other_key],
                                              unassigned_neighbours2[old_key]]
                                del unassigned_neighbours2[old_key]
                            
                            new_neighbours = list(chain.from_iterable(new_neighbours))
                            new_neighbours = list(set(new_neighbours))
                            new_neighbours.remove(unassigned_index)
                            if other_key in new_neighbours:
                                new_neighbours.remove(other_key)
                            for ok in old_key:
                                if ok in new_neighbours:
                                    new_neighbours.remove(ok)
                            new_key = list(old_key)
                            new_key.append(unassigned_index)
                            new_key.append(other_key)
                            unassigned_neighbours2[tuple(set(new_key))] = new_neighbours
                        else:
                            new_neighbours = [unassigned_neighbours[unassigned_index],
                                              unassigned_neighbours[other_key]]
                            new_neighbours = list(chain.from_iterable(new_neighbours))
                            new_neighbours.remove(unassigned_index)
                            if other_key in new_neighbours:
                                new_neighbours.remove(other_key)
                            new_neighbours = list(set(new_neighbours))
                            unassigned_neighbours2[(unassigned_index, other_key)] = new_neighbours
                else:
                    if [key for key in unassigned_neighbours2.keys() if unassigned_index in key] == []:
                        unassigned_neighbours2[(unassigned_index,)] = unassigned_neighbours[unassigned_index]

        setting_unassigned(sss,
                           smiles_mol,
                           unassigned_neighbours2)

    return None

def establishing_smarts(ReadingFiles):
    ''' Function: establishing_smarts
    This function calls for the smarts dictionary to be made with it's classification
    if the program is in verbose mode the amount of each category is printed
    input: ReadingFiles - class
    output: None
    '''
    smarts_dictionary = ReadingFiles.smarts_file()
    donor_count = 0
    acceptor_count = 0
    acid_count = 0
    base_count = 0
    halo_count = 0
    charged_count = 0
    donor_acceptor_count = 0

    for value in smarts_dictionary.values():
        node_type = value['Type']
        if node_type == 'HBD':
            donor_count += 1
        elif node_type == 'HBA':
            acceptor_count += 1
        elif node_type == 'acid':
            acid_count += 1
        elif node_type == 'base':
            base_count += 1
        elif node_type == 'halo':
            halo_count += 1
        elif node_type == 'anion' or node_type == 'cation' or node_type == 'pos' or node_type == 'neg':
            charged_count += 1
        elif node_type == ('HBD HBA' or 'HBA HBD'):
            donor_acceptor_count += 1
        else:
            print(node_type)
            raise ValueError('Unknown smarts type')

    count_list = []
    count_list.append("Donors "+ str(donor_count))
    count_list.append("Acceptors "+ str(acceptor_count))
    count_list.append("Acid "+ str(acid_count))
    count_list.append("Bases "+ str(base_count))
    count_list.append("Halogens "+ str(halo_count))
    count_list.append("Anions or Cations "+ str(charged_count))
    count_list.append("Donor and Acceptors " + str(donor_acceptor_count))
    if args.verbose:
        print('\n'.join(count_list))

    return None

def establishing_isolating_c(ReadingFiles):
    ''' Function: establishing_isolating_c
    This function reads in the smarts_isolating_file and adds to the smarts dictionary
    input: ReadingFiles - class
    output: smarts_dictionary - {} dictionary of all the different node types smarts
    '''
    smarts_dictionary = ReadingFiles.smarts_isolating_file()

    return smarts_dictionary

def finding_nodes_to_combine(smiles_mol,
                             list_of_combination,
                             list_of_keys,
                             neighbour,
                             neighbour_index):
    '''Fucntion: finding_nodes_to_combine
    This function finds nodes that are neighbours to one another and puts them
    in the list_of_combination
    input: smiles_mol - mol object of the smiles that is being searched
           list_of_combination - [] of [] of tuples that need to be combined
           list_of_keys - [] of tuples of atom indexes of nodes
           neighbour - index of the neighbour atom
           neighbour_index - position of the neighbour within the
                             list_of_neighbour_atoms
    output: None
    '''
    matching_list = []
    if type(neighbour) == tuple:
        matching_list.append([item for atom in neighbour for item in list_of_keys if atom in item])
    elif type(neighbour) == int:
        matching_list.append([item for item in list_of_keys if neighbour in item])
    ### removing empty lists from matching list
    matching_list = [list1 for list1 in matching_list if list1]
    if len(matching_list) >= 1:
        if len(list_of_keys[neighbour_index]) == 1:
            if [smiles_mol.GetAtomWithIdx(atom_index).IsInRing() for
                    atom_index in list_of_keys[neighbour_index]] == [False]:
                matching_list = [item for sublist in matching_list for item in sublist]
                matching_list.append(list_of_keys[neighbour_index])
                if len(list_of_combination) == 0:
                    combo_tuples = matching_list
                elif len(list_of_combination) >= 1:
                    for combo in range(0, len(list_of_combination)):
                        for match in range(0, len(matching_list)):
                            if matching_list[match] in list_of_combination[combo]:
                                list_of_combination[combo] += matching_list
                                list_of_combination[combo] = list(set(list_of_combination[combo]))
                                list_of_combination[combo].sort()
                            else:
                                combo_tuples = matching_list
        else:
            matching_list.append(list_of_keys[neighbour_index])
            combo_tuples = []
            for combo in range(0, len(matching_list)):
                if type(matching_list[combo]) == list:
                    combo_tuples.append(matching_list[combo][0])
                elif type(matching_list[combo]) == tuple:
                    combo_tuples.append(matching_list[combo])
    combo_tuples = list(set(combo_tuples))
    combo_tuples.sort()
    if combo_tuples not in list_of_combination:
        further_combo = [tuples_to_combine for atom_indexes in combo_tuples for
                         tuples_to_combine in list_of_combination if atom_indexes in tuples_to_combine]
        if further_combo == []:
            list_of_combination.append(combo_tuples)
        else:
            new_combo_tuples = [combo_tuples]
            new_combo_tuples.extend([combination for combination in further_combo])
            new_combo_tuples = list(set(chain.from_iterable(new_combo_tuples)))
            list_of_combination = [tuple_ for tuple_ in list_of_combination if tuple_ not in further_combo]
            list_of_combination.append(new_combo_tuples)

    return list_of_combination

def checking_for_connected(sss,
                           smiles_mol):
    ''' Function: checking_for_connected
    This function check to see whether any HBA, HBD or HBD/HBA are connected
    to each other and then if so they are combined into the corresponding node
    type, i.e. carbonyl and N are combined to make an amide
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
    output: None
    '''
    ### Combines two keys together if they contain the same integer
    items_to_combine = []
    for node_key in sss.dictionary_of_nodes.keys():
        for atom_index in node_key:
            number_of_appearances_list = [key for key in sss.dictionary_of_nodes.keys()
                                          if atom_index in key]
            if len(number_of_appearances_list) > 1:
                if number_of_appearances_list not in items_to_combine:
                    items_to_combine.append(number_of_appearances_list)
    set_combined_node(sss,
              smiles_mol,
              items_to_combine)
    ### This function checks for things like amides
    list_of_neighbour_atoms = [value['neighbour_index'] for key, value in
                               sss.dictionary_of_nodes.items() if 'acyclic' in value['node_type']]
    list_of_keys = [key for key, value in sss.dictionary_of_nodes.items()
                    if 'acyclic' in value['node_type']]
    list_of_combination = []

    for neighbour_index in range(0, len(list_of_neighbour_atoms)):
        for neighbour in list_of_neighbour_atoms[neighbour_index]:
            try:
                list_of_combination = finding_nodes_to_combine(smiles_mol,
                                                               list_of_combination,
                                                               list_of_keys,
                                                               neighbour,
                                                               neighbour_index)
            except:
                pass

    list_of_combination = [m for i, m in enumerate(list_of_combination) if not
                           any(set(m).issubset(set(n)) for n in
                               (list_of_combination[:i] + list_of_combination[i+1:]))]
    ### this searches through the list_of_combination to see if any tuple appears
    ### multiple times if it does it then combines these instances
    list_of_combination2 = []
    for keys_to_combine in list_of_combination:
        new_key = []
        for intial_key in keys_to_combine:
            items_to_combine = [sublist for sublist in list_of_combination if intial_key in sublist]
            if len(items_to_combine) > 1:
                to_combine = list(chain.from_iterable(items_to_combine))
                to_combine = list(set(tuple(row) for row in to_combine))
                if to_combine not in list_of_combination2:
                    for tup in to_combine:
                        items_to_delete = [sublist for sublist in list_of_combination2 if tup in sublist]
                        if items_to_delete != []:
                            if len(items_to_delete) == 1:
                                list_of_combination2.remove(items_to_delete[0])
                            elif len(items_to_delete) > 1:
                                for index in range(0, len(items_to_delete)):
                                    list_of_combination2.remove(items_to_delete[index][0])
                            items_to_delete = list(set(tuple(row) for row in items_to_delete))
                            items_to_delete = list(chain.from_iterable(items_to_delete))
                            new_key += items_to_delete
                new_key += to_combine
            else:
                items_to_combine = list(chain.from_iterable(items_to_combine))
                new_key += items_to_combine
        new_key = list(set(tuple(row) for row in new_key))
        new_key.sort()
        new_key = [m for i, m in enumerate(new_key) if not any(set(m).issubset(set(n)) for n in (new_key[:i] + new_key[i+1:]))]
        list_of_combination2.append(new_key)

    set_combined_node(sss,
              smiles_mol,
              list_of_combination2)

    return None

def generate_reduced_graph(copy_mol):
    ''' Function: generate_reduced_graph
    This function makes all the atoms explicit and therefore gets no hydrogen
    in the smiles string and then generates the smiles string of the reduced graph
    input: copy_mol - mol object of the reduced graph
    output: reduced_graph - string of the smiles of the reduced graph
    '''
    ### Viewed structure in molblock and hydrogens are not present this means that
    ### the hydrogens only become present in the moltosmiles operation therefore
    ### need to sort this out via string
    for atom_index in range(0, copy_mol.GetNumAtoms()):
        atom = copy_mol.GetAtomWithIdx(atom_index)
        atom.SetNoImplicit(True)

    if args.verbose:
        print('Generating Reduced Graph')
        print(Chem.MolToSmiles(copy_mol, isomericSmiles=True))
    reduced_graph = Chem.MolToSmiles(copy_mol, isomericSmiles=True)

    return reduced_graph

def halogen_neighbours(sss,
                       smiles_mol,
                       list_of_halogen_nodes):
    ''' Function: halogen_neighbours
    This function searches to see if the halogen node is neighboured to a
    hbd or hbd or hbd/hba node, if halogen node is set(C,F) then add to that node
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           list_of_halogen_nodes - [] of tuples of halogen nodes
    output: list_of_combo - [] of tuples of atom indexes
            list_of_neighbour_atoms - [] of {} of the neighbours indexes and
                                    bond types
    '''
    list_of_neighbour_atoms = []
    list_to_delete = []
    for halogen_node in list_of_halogen_nodes:
        all_neighbours = {}
        all_neighbours_index = []
        for atom in halogen_node:
            neighbours, connecting_bonds_list, doubley_bonded_atoms_list = defining_neighbours(smiles_mol,
                                                                                               list(halogen_node),
                                                                                               atom)
            if neighbours != {}:
                all_neighbours.update(neighbours)
                if type(neighbours[atom]['neighbour']) == int:
                    all_neighbours_index.append(neighbours[atom]['neighbour'])
                else:
                    all_neighbours_index.extend(list(neighbours[atom]['neighbour']))
        hbahbd_keys = [key for key, value in sss.dictionary_of_nodes.items() if
                       ('HBA' in value['node_type'] or 'HBD' in value['node_type'])]
        if len(halogen_node) == 4 and len(all_neighbours) == 1 and [key for key in hbahbd_keys if key in all_neighbours.keys()] != []: #[0] in key
            additional_smarts = [smiles_mol.GetAtomWithIdx(atom).GetSmarts() for atom in halogen_node]
            if set(additional_smarts) == {'C', 'F'} or set(additional_smarts) == {'F', 'C'}:
                list_to_delete.append(list_of_halogen_nodes.index(halogen_node))
                key_to_combine = [key for key in hbahbd_keys if all_neighbours[0] in key]
                if type(key_to_combine[0]) == int:
                    new_key = (key_to_combine[0],) + halogen_node
                elif type(key_to_combine[0]) == tuple:
                    new_key = key_to_combine[0] + halogen_node
                sss.dictionary_of_nodes[new_key] = sss.dictionary_of_nodes[key_to_combine[0]]
                sss.dictionary_of_nodes[new_key]['atoms_index'] = new_key
                sss.dictionary_of_nodes[new_key]['atoms'].append(additional_smarts)
                sss.dictionary_of_nodes[new_key]['node_size'] = len(new_key)
                sss.dictionary_of_nodes[new_key]['nr_atoms'] = len(new_key)
                del sss.dictionary_of_nodes[key_to_combine[0]]
                all_neighbours = {}
                all_neighbours_index = []
                connecting_bonds = []
                double_bond_neighbours = []

                for atom_index in new_key:
                    neighbours, connecting_bonds_list, doubley_bonded_atoms_list = defining_neighbours(smiles_mol,
                                                                                                       new_key,
                                                                                                       atom_index)
                    if neighbours != {}:
                        all_neighbours.update(neighbours) #append
                        all_neighbours_index.append(neighbours[atom_index]['neighbour'])
                        connecting_bonds.extend(connecting_bonds_list)
                        double_bond_neighbours.extend(doubley_bonded_atoms_list)
                sss.dictionary_of_nodes[new_key]['connecting_bonds'] = connecting_bonds
                sss.dictionary_of_nodes[new_key]['double_bonds_connected'] = double_bond_neighbours
                sss.dictionary_of_nodes[new_key]['neighbours'] = all_neighbours
                sss.dictionary_of_nodes[new_key]['neighbour_index'] = all_neighbours_index
        list_of_neighbour_atoms.append(all_neighbours_index)
    list_to_delete.sort(reverse=True)

    for index in list_to_delete:
        del list_of_halogen_nodes[index]
        del list_of_neighbour_atoms[index]
    list_of_combo = []

    for item in range(0, len(list_of_neighbour_atoms)):
        for node in range(0, len(list_of_halogen_nodes)):
            for atom in list_of_neighbour_atoms[item]:
                if atom in list_of_halogen_nodes[node]:
                    if ([list_of_halogen_nodes[item],
                         list_of_halogen_nodes[node]] in list_of_combo) is False and ([list_of_halogen_nodes[node], list_of_halogen_nodes[item]] in list_of_combo) is False:
                        list_of_combo.append([list_of_halogen_nodes[node], list_of_halogen_nodes[item]])

    return(list_of_combo, list_of_neighbour_atoms)

def combination_of_tuples(list_of_combo):
    ''' Fucntion: combination_of_tuples
    This function searches through the list_of_combination to see if any tuple
    appears multiple times if it does then the instances get combined
    input: list_of_combo - [] of tuples of atom indexes
    output: list_of_combo2 - [] of [] of tuples of atom indexes
    '''
    list_of_combo2 = []
    for key in list_of_combo:
        new_key = []
        for intial_key in key:
            items_to_combine = [sublist for sublist in list_of_combo if intial_key in sublist]
            if len(items_to_combine) > 1:
                to_combine = list(chain.from_iterable(items_to_combine))
                to_combine = list(set(tuple(row) for row in to_combine))
                if to_combine not in list_of_combo2:
                    for tup in to_combine:
                        items_to_delete = [sublist for sublist in list_of_combo2 if tup in sublist]
                        if items_to_delete != []:
                            if len(items_to_delete) == 1:
                                list_of_combo2.remove(items_to_delete[0])
                            elif len(items_to_delete) > 1:
                                for index in range(0, len(items_to_delete)):
                                    list_of_combo2.remove(items_to_delete[index][0])
                            items_to_delete = list(set(tuple(row) for row in items_to_delete))
                            items_to_delete = list(chain.from_iterable(items_to_delete))
                            new_key += items_to_delete
                new_key += to_combine
            else:
                items_to_combine = list(chain.from_iterable(items_to_combine))
                new_key += items_to_combine
        new_key = list(set(tuple(row) for row in new_key))
        new_key = [node_key for index, node_key in enumerate(new_key) if not any(set(node_key).issubset(set(smaller_node_key)) for smaller_node_key in (new_key[:index] + new_key[index+1:]))]
        if new_key not in list_of_combo2:
            list_of_combo2.append(new_key)

    return list_of_combo2

def establishing_halogens(sss,
                          smiles_mol):
    ''' Function: establishing_halogens
    This function establishes if there are any halogens within the molecule and
    whether they are close to another halogens
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
    output: None
    '''
    matches_list = find_halogens(smiles_mol)
    current_keys = list(sss.dictionary_of_nodes.keys())
    current_keys = list(chain.from_iterable(current_keys))
    delete_list = [match_atom for match_tuple in matches_list
                   for match_atom in match_tuple if match_atom in current_keys or
                   smiles_mol.GetAtomWithIdx(match_atom).IsInRing() is True]
    matches_list = [match_tuple for match_tuple in matches_list if match_tuple[0] not in set(delete_list)]
    remove_list, list_of_halogen_nodes, halogens_neighbours = set_halogen_node(sss,
                                                                               smiles_mol,
                                                                               matches_list)

    matches_list = [halogen_index for halogen_index in matches_list if halogen_index not in remove_list]
    list_of_neighbour_atoms = []
    list_of_keys = [key for halogen_match in halogens_neighbours for key in halogen_match.keys()]

    for halogen in halogens_neighbours:
        for value in halogen.values():
            if type(value['neighbour']) == int:
                if smiles_mol.GetAtomWithIdx(value['neighbour']).IsInRing() is True:
                    list_of_neighbour_atoms.append([])
                else:
                    list_of_neighbour_atoms.append(value['neighbour'])
            elif type(value['neighbour']) == tuple:
                neighbour_list = [value['neighbour'][atom] for atom in range(0, len(value['neighbour'])) if
                                  smiles_mol.GetAtomWithIdx(value['neighbour'][atom]).IsInRing() is False]
                list_of_neighbour_atoms.append(list(set(neighbour_list)))

    for key in list_of_keys:
        list_of_nodes = [item for item in list_of_keys if item in list_of_neighbour_atoms]
        if list_of_nodes == []:
            list_of_halogen_nodes.append((key,))
        elif list_of_nodes != []:
            list_of_halogen_nodes.append(tuple(list_of_nodes))
    indexes_of_halogens_to_combine = []

    for neighbour_index in [x for x in list_of_neighbour_atoms if x]:
        indexes_of_halogens_to_combine.append([i for i, x in enumerate(list_of_neighbour_atoms) if x == neighbour_index])
    indexes_of_halogens_to_combine = [s for s in indexes_of_halogens_to_combine if len(s) > 1]
    indexes_of_halogens_to_combine = list(set(tuple(i_l) for i_l in indexes_of_halogens_to_combine))

    for list_index in indexes_of_halogens_to_combine:
        list_to_append = []
        for atom_index in list_index:
            list_to_append.append(list_of_keys[atom_index])
            list_to_append.append(list_of_neighbour_atoms[atom_index])
            sss.halogens_atoms[list_of_neighbour_atoms[atom_index]] = 1
            sss.assigned_atoms[list_of_neighbour_atoms[atom_index]] = 1
            key_tuple = (list_of_keys[atom_index],)
            list_of_halogen_nodes.remove(key_tuple)
        list_of_halogen_nodes.append(tuple(set(list_to_append)))
    list_of_combo, list_of_neighbour_atoms = halogen_neighbours(sss,
                                                                smiles_mol,
                                                                list_of_halogen_nodes)
    list_of_combo2 = combination_of_tuples(list_of_combo)
    for combo in list_of_combo2:
        new_key = tuple(chain.from_iterable(combo))
        list_of_halogen_nodes.append(new_key)
        list_of_halogen_nodes = [tuple_to_combine for tuple_to_combine in list_of_halogen_nodes if tuple_to_combine not in combo]

    for item in list_of_halogen_nodes:
        if args.hydrophobic:
            set_basic_node(sss,
                           smiles_mol,
                           item,
                           'hydrophobic')
        else:
            set_basic_node(sss,
                           smiles_mol,
                           item,
                           'inert')

    return None

def list_atomic_numbers(smiles_mol):
    ''' Function: list_atomic_numbers
    This function prints to screen each atom and it's corresponding atom index
    input: smiles_mol - mol object of the smiles that is being searched
    output: None
    '''
    atom_number = 0
    while atom_number < smiles_mol.GetNumAtoms():
        atom_mol = smiles_mol.GetAtomWithIdx(atom_number)
        atom = atom_mol.GetSmarts()
        print('Atom '+ str(atom_number) + ' is ' + atom)
        atom_number += 1

    return None

def finding_predefined_nodes(sss, smiles_mol):
    pre_defined_atoms = [-1]*smiles_mol.GetNumAtoms()
    if args.metal:
        matches = sss.find_metal_atoms(smiles_mol)
        for match in matches:
            set_basic_node(sss,
                           smiles_mol,
                           match,
                           'metal')
            
    matches = sss.find_atoms(smarts_dictionary,
                             smiles_mol,
                             'acid')
    if len(matches) >= 1:
        for match in matches:
            set_nodes(sss,
                      smiles_mol,
                      match,
                      'acid')
    matches = sss.find_atoms(smarts_dictionary,
                             smiles_mol,
                             'base')
    if len(matches) >= 1:
        for match in matches:
            set_nodes(sss,
                      smiles_mol,
                      match,
                      'base')
    matches = sss.find_atoms(smarts_dictionary,
                             smiles_mol,
                             'pos')
    matches_to_delete = [node_key for index, node_key in enumerate(matches) if any(set(node_key).issubset(set(smaller_node_key)) for smaller_node_key in (matches[:index] + matches[index+1:]))]
    matches = [match for match in matches if match not in matches_to_delete]
    if len(matches) >= 1:
        for match in matches:
            set_basic_node(sss,
                           smiles_mol,
                           match,
                           'positive')
    matches = sss.find_atoms(smarts_dictionary,
                             smiles_mol,
                             'neg')
    matches_to_delete = [node_key for index, node_key in enumerate(matches) if any(set(node_key).issubset(set(smaller_node_key)) for smaller_node_key in (matches[:index] + matches[index+1:]))]
    matches = [match for match in matches if match not in matches_to_delete]
    if len(matches) >= 1:
        for match in matches:
            set_basic_node(sss,
                           smiles_mol,
                           match,
                           'negative')
    for match in sss.dictionary_of_nodes.keys():
        for m in match:
            pre_defined_atoms[m] = 1
    
    matches = sss.find_atoms(smarts_dictionary,
                             smiles_mol,
                             'HBD')
    matches_to_delete = [node_key for index, node_key in enumerate(matches) if any(set(node_key).issubset(set(smaller_node_key)) for smaller_node_key in (matches[:index] + matches[index+1:]))]
    matches = [match for match in matches if match not in matches_to_delete]
    if len(matches) >= 1:
        for match in list(matches):
            if len(match) == 1:
                if pre_defined_atoms[match[0]] == 1:
                    matches.remove(match)
            elif len(match) > 1:
                for match_atom in match:
                    if pre_defined_atoms[match_atom] == 1:
                        node_keys = [key for key in sss.dictionary_of_nodes.keys() if match_atom in key]
                        if node_keys != []:
                            for item in node_keys:
                                if match == item:
                                    continue
                                elif set(item).issubset(match) is True:
                                    if item in sss.dictionary_of_nodes.keys():
                                        del sss.dictionary_of_nodes[item]
                                elif set(item).issuperset(match) is True:
                                    if match in matches:
                                        matches.remove(match)
                                else:
                                    if match in matches:
                                        matches.remove(match)
                                    match = reassigning_node(args,
                                                             sss,
                                                             smiles_mol,
                                                             item,
                                                             match)
        set_nodes(sss,
                  smiles_mol,
                  matches,
                  'HBD')

    matches = sss.find_atoms(smarts_dictionary,
                             smiles_mol,
                             'HBA')
    matches_to_delete = [node_key for index, node_key in enumerate(matches) if any(set(node_key).issubset(set(smaller_node_key)) for smaller_node_key in (matches[:index] + matches[index+1:]))]
    matches = [match for match in matches if match not in matches_to_delete]
    if len(matches) >= 1:
        new_matches = copy.deepcopy(matches)
        for match in new_matches:
            if len(match) == 1:
                if pre_defined_atoms[match[0]] == 1:
                    matches.remove(match)
            elif len(match) > 1:
                for match_atom in match:
                    if pre_defined_atoms[match_atom] == 1:
                        node_keys = [key for key in sss.dictionary_of_nodes.keys() if match_atom in key]
                        if node_keys != []:
                            for item in node_keys:
                                if match == item:
                                    continue
                                elif set(item).issubset(match) is True:
                                    if item in sss.dictionary_of_nodes.keys():
                                        del sss.dictionary_of_nodes[item]
                                elif set(item).issuperset(match) is True:
                                    if match in matches:
                                        matches.remove(match)
                                else:
                                    if match in matches:
                                        matches.remove(match)
                                    match = reassigning_node(args,
                                                             sss,
                                                             smiles_mol,
                                                             item,
                                                             match)
        set_nodes(sss,
                  smiles_mol,
                  matches,
                  'HBA')

    return None

def finding_feature_atoms(smarts_dictionary,
                          can_smiles,
                          rf):
    ''' Function: finding_feature_atoms
    This function is where most of the fucntions to find the nodes is called
    input: smarts_dictionary - {} dictionary of all the different node types smarts
           can_smiles - string of the canonicalised smiles
           rf - class
    output: None
    '''
    smiles_mol = Chem.MolFromSmiles(can_smiles)

    sss = SubStructureSearching(can_smiles,
                                smiles_mol)
    finding_predefined_nodes(sss, smiles_mol)
    checking_for_connected(sss,
                           smiles_mol)
    establishing_halogens(sss,
                          smiles_mol)
    establishing_ring_nodes(args,
                            sss,
                            smiles_mol)
    establishing_isolating_c_nodes(sss,
                                   smiles_mol,
                                   smarts_dictionary)
    assigning_unassigned(sss,
                         smiles_mol)
    copy_mol, mol_block = generating_rg_from_scratch(args,
                                                     sss.dictionary_of_nodes,
                                                     smiles_mol,
                                                     rf.smiles_dictionary[can_smiles]['Name'])

    reduced_graph = generate_reduced_graph(copy_mol)

    rf.smiles_dictionary[can_smiles]['reduced_graph'] = reduced_graph

    return mol_block

if __name__ == '__main__':
    start_time = time.time()
    args = input_args()
    ReadingFiles = ReadingFiles()
    establishing_smarts(ReadingFiles)
    smarts_dictionary = establishing_isolating_c(ReadingFiles)
    
    if args.additionalpIC50:
        try:
            pIC50_values = {}
            with open(args.pIC50inputfile, 'r') as pIC50line:
                for line in pIC50line:
                    pline = line.split(args.separator)
                    pIC50_values[pline[1]] = pline[-1]
        except:
            raise ValueError('Cannot open pIC50 file')
        
    try:
        smi = open(args.smiles, 'r')
    except:
        raise ValueError('Cannot open SMILES files')

    if args.verbose:
        num_lines = sum(1 for line in open(args.smiles))
        print(''.join(['Number of SMILES contained within the file: ',
                       str(num_lines)]))

    sdf_writer = Chem.SDWriter(args.sdf)
    output_file = open(args.output, 'w+')
    output_file.write("SMILES\tID\tRG\n")

    molecules_evaulated = 0
    for smile in smi:
        if smile != 'SMILES\tID\tpIC50' or smile != 'SMILES\tID\tpIC50\tRound':
            smiles_line = smile.replace("\n", "").split(args.separator)
            can_smiles = canonicalisor.canonicalise_smiles(smiles_line[0])
            if type(can_smiles) == tuple:
                ReadingFiles.smiles_unable_to_read.append(smiles_line[0])
            elif str(can_smiles) == 'nan' or str(can_smiles) == '':
                ### removes molecules that do not have any SMILES
                ReadingFiles.smiles_unable_to_read.append(smiles_line[1])
            elif can_smiles != False:
                if ReadingFiles.smiles_dictionary.get(can_smiles) is None:
                    ReadingFiles.smiles_dictionary[can_smiles] = {}
                    ReadingFiles.smiles_dictionary[can_smiles]['Name'] = smiles_line[1]
                    mol_block = finding_feature_atoms(smarts_dictionary,
                                                      can_smiles,
                                                      ReadingFiles)
    
                    if args.additionalpIC50:
                        mol_block.SetProp('pIC50', pIC50_values[smiles_line[1]])
                    sdf_writer.write(mol_block)
                    molecules_evaulated += 1
                    molecules_number = ''.join(['Molecules Evaulated ',
                                                str(molecules_evaulated)])
                    print(molecules_number)
                    output_file.write("%s\t%s\t%s\n" % (can_smiles,
                                                        ReadingFiles.smiles_dictionary[can_smiles]['Name'],
                                                        ReadingFiles.smiles_dictionary[can_smiles]['reduced_graph']))
                else:
                    ReadingFiles.duplicate_smiles.append(can_smiles)
                    print('ValueError '+ can_smiles +' appears twice will only use once')
            elif can_smiles is False:
                ReadingFiles.smiles_unable_to_read.append(smiles_line[0])#2

    sdf_writer.close()
    output_file.close()

    finish_time = time.time()
    print(finish_time - start_time,
          ' to complete ',
          len(ReadingFiles.smiles_dictionary),
          ' molecules')
    if len(ReadingFiles.duplicate_smiles) >= 1:
        print('Number of duplicate smiles within the dataset: ',
              len(ReadingFiles.duplicate_smiles))
    if len(ReadingFiles.smiles_unable_to_read) >= 1:
        print('Number of smiles within the dataset that couldn\'t be read: ',
              len(ReadingFiles.smiles_unable_to_read))

    print('Program Finished')
    