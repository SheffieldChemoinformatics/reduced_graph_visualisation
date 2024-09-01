##############################################################################
# Reduced graph program
# Setting Nodes
#
#
# Jess Stacey
##############################################################################

import copy
from itertools import chain
from rdkit import Chem

from node_codes import node_codes_dict as node_codes

def set_node_dictionary():
    ''' Function: set_node_dictionary
    This function creates the standard node dictionary with all values not set
    input: None
    output: node_dictionary - {} dictionary of node
    '''
    node_dictionary = {'aromatic' : False,
                       'atoms' : None,
                       'atoms_index' : None,
                       'connecting_bonds' : [],
                       'fused_neighbours' : None,
                       'neighbour_index': None,
                       'neighbours' : None,
                       'node_code' : 0,
                       'node_label' : None,
                       'node_size' : 0,
                       'node_type': -1,
                       'nr_fused_neighbours' : 0,
                       'nr_neighbours' : 0,
                       'nr_rings' : 0,
                       'representing_atom' : None,
                       'node_number' : None,
                       'double_bonds_connected' : []}

    return node_dictionary

def defining_neighbours(smiles_mol,
                        list_of_atoms,
                        atom_index):
    ''' Function: defining_neighbours
    This function finds the neighbours and bond types of an atom
    these neighbours are not part of the list of atoms within the node
    input: smiles_mol - mol object of the smiles that is being searched
           list_of_atoms - [] of atom indexes that are contained within the node
           atom_index - atom_index of the atom that is having it's neighbours
                        investigated
    output: neighbours_dictionary - {} of atom indexes that are neighbours
                                of the searched atom along with there bond type
            connecting_bonds_list - [] of bond indexes that are connected to the node
    '''
    neighbours_dictionary = {}
    other_atoms_ids = (bond.GetOtherAtomIdx(atom_index) for bond in \
                       smiles_mol.GetAtomWithIdx(atom_index).GetBonds())
    all_neighbours = []
    connecting_bonds_list = []
    doubley_bonded_atoms_list = []
    for other_id in other_atoms_ids:
        if other_id not in set(list_of_atoms):
            connecting_bonds_list.extend([item.GetIdx() for item in smiles_mol.GetAtomWithIdx(other_id).GetBonds()])
            if atom_index in set(neighbours_dictionary.keys()):
                all_neighbours.append(other_id)
                if type(neighbours_dictionary[atom_index]['neighbour']) == tuple:
                    new_neighbours = [neighbour for neighbour in neighbours_dictionary[atom_index]['neighbour']]
                    new_neighbours.append(other_id)
                    new_bond_type = [bond_type for bond_type in neighbours_dictionary[atom_index]['bond_type']]
                    new_bond_type.append(smiles_mol.GetBondBetweenAtoms(atom_index, other_id).GetBondType())
                    neighbours_dictionary[atom_index] = {'neighbour': tuple(new_neighbours),
                                                         'bond_type': tuple(new_bond_type)}
                elif type(neighbours_dictionary[atom_index]['neighbour']) == int:
                    new_neighbours = [neighbours_dictionary[atom_index]['neighbour'],
                                      other_id]
                    new_bond_type = [neighbours_dictionary[atom_index]['bond_type'],
                                     smiles_mol.GetBondBetweenAtoms(atom_index, other_id).GetBondType()]
                    neighbours_dictionary[atom_index] = {'neighbour': tuple(new_neighbours),
                                                         'bond_type': tuple(new_bond_type)}
            else:
                neighbours_dictionary[atom_index] = {'neighbour': other_id,
                                                     'bond_type': smiles_mol.GetBondBetweenAtoms(atom_index, other_id).GetBondType()}
            if smiles_mol.GetBondBetweenAtoms(atom_index, other_id).GetBondType() == Chem.BondType.DOUBLE:
                doubley_bonded_atoms_list.append(other_id)

    return(neighbours_dictionary, connecting_bonds_list, doubley_bonded_atoms_list)

def establishing_aromaticity(node_dictionary,
                             node_type,
                             aromaticity,
                             atom_ring_system):
    ''' Function: establishing_aromaticity
    This function establishes the aromaticity of the node
    input: node_dictionary - {} the node dictionary
           node_type - string of the type of node
           aromaticity - [] of each atom in the nodes aromaticity
           atom_ring_system - [] of whether each atom in the node is in a ring or not
    output: node_string - string of the node type which now includes the aromaticity
    '''
    if len(set(aromaticity)) == 1:
        node_dictionary['aromatic'] = aromaticity[0]
        if node_dictionary['aromatic'] is True:
            node_list = ['aromatic', node_type]
        elif node_dictionary['aromatic'] is False:
            if len(set(atom_ring_system)) == 1 and set(atom_ring_system) == {True}:
                node_list = ['aliphatic', node_type]
            else:
                node_list = ['acyclic', node_type]
    else:
        node_dictionary['aromatic'] = False
        if len(set(atom_ring_system)) == 1 and set(atom_ring_system) == {True}:
            node_list = ['aliphatic', node_type]
        else:
            node_list = ['acyclic', node_type]
    node_string = ' '.join(node_list)

    return node_string

def connection_of_interest(sss,
                           smiles_mol,
                           neighbours,
                           atom_index,
                           neighbour,
                           bond_type,
                           tuple_atoms_index,
                           neighbour_dict_list,
                           neighbour_list):
    '''Function: connection_of_interest
    This function looks to see if the atom is connected by doubley or tripley bonded to another atom
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           neighbours - {} of atom indexes that are neighbours of the searched atom
                   along with there bond type
           atom_index - integer of the atom index number
           neighbour - atom index of the neighbour
           bond_type - bond type of the bond that connects atom_index and neighbour
           tuple_atoms_index - () tuple of the atom indexes that are within that node
           neighbour_dict_list - [] of neighbour dictionaries
           neighbour_list -[] of atom indexes of the neighbours
    output: neighbours - {} of atom indexes that are neighbours of the searched atom
                   along with there bond type
            tuple_atoms_index - () tuple of the atom indexes that are within that node
            bool_set - true or false depending on whether the node has been set or not
    '''
    bool_set = False
    if (sss.metal_atoms[neighbour] == -1 and
            (bond_type == Chem.rdchem.BondType.DOUBLE or
             bond_type == Chem.rdchem.BondType.TRIPLE)):
        bool_set = True
        if smiles_mol.GetAtomWithIdx(neighbour).IsInRing() is False:
            sss.assigned_atoms[neighbour] = 1
            tuple_atoms_index = tuple_atoms_index + (neighbour,)
            additional_neighbours, connecting_bonds_list, doubley_bonded_atoms_list = defining_neighbours(smiles_mol,
                                                                                                          tuple_atoms_index,
                                                                                                          neighbour)
            if bool(additional_neighbours) is True:
                neighbour_dict_list.append(additional_neighbours)
                if type(additional_neighbours[neighbour]['neighbour']) == int:
                    neighbour_list.append(additional_neighbours[neighbour]['neighbour'])
                else:
                    neighbour_list.extend(list(additional_neighbours[neighbour]['neighbour']))
            if type(neighbour) == int:
                neighbour = {}
            elif len(neighbour) == 1:
                neighbours = {}
            else:
                neighbours[atom_index]['neighbour'] = list(neighbours[atom_index]['neighbour'])
                neighbours[atom_index]['bond_type'] = list(neighbours[atom_index]['bond_type'])
                neighbours[atom_index]['neighbour'].remove(neighbour)
                neighbours[atom_index]['bond_type'].remove(neighbour)
                if neighbours[atom_index]['neighbour'] == 1:
                    neighbours[atom_index]['neighbour'] = int(neighbours[atom_index]['neighbour'])
                    neighbours[atom_index]['bond_type'] = str(neighbours[atom_index]['bond_type'])
                else:
                    neighbours[atom_index]['neighbour'] = tuple(neighbours[atom_index]['neighbour'])
                    neighbours[atom_index]['bond_type'] = tuple(neighbours[atom_index]['bond_type'])

    return(neighbours, tuple_atoms_index, bool_set)

def potential_sulphoxide(sss,
                         smiles_mol,
                         neighbours,
                         atom_index,
                         tuple_atoms_index,
                         neighbour_dict_list,
                         neighbour_list):
    '''Function: potential_sulphoxide
    This function looks to see whether there are any sulphoroxides present in the molecule
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           neighbours - {} of atom indexes that are neighbours of the searched atom
                   along with there bond type
           atom_index - integer of the atom index number
           tuple_atoms_index - () tuple of the atom indexes that are within that node
           neighbour_dict_list - [] of neighbour dictionaries
           neighbour_list - [] of atom indexes of the neighbours
    output: tuple_atoms_index - () tuple of the atom indexes that are within that node
            neighbours - {} of atom indexes that are neighbours of the searched atom
                   along with there bond type
    '''
    sss.assigned_atoms[neighbours[atom_index]['neighbour']] = 1
    tuple_atoms_index = tuple_atoms_index + (neighbours[atom_index]['neighbour'],)
    more_neighbours, connecting_bonds_list, doubley_bonded_atoms_list = defining_neighbours(smiles_mol,
                                                                                            tuple_atoms_index,
                                                                                            neighbours[atom_index]['neighbour'])
    if type(more_neighbours[neighbours[atom_index]['neighbour']]['neighbour']) == int:
        if (more_neighbours[neighbours[atom_index]['neighbour']]['bond_type'] == Chem.rdchem.BondType.DOUBLE and
                smiles_mol.GetAtomWithIdx(more_neighbours[neighbours[atom_index]['neighbour']]['neighbour']).IsInRing() is False and
                smiles_mol.GetAtomWithIdx(more_neighbours[neighbours[atom_index]['neighbour']]['neighbour']).GetSmarts() == 'C'):
            tuple_atoms_index = tuple_atoms_index + (more_neighbours[neighbours[atom_index]['neighbour']]['neighbour'],)
            more_neighbours = {}
    elif type(more_neighbours[neighbours[atom_index]['neighbour']]['neighbour']) == tuple:
        for neighbour_index in range(0, len(more_neighbours[neighbours[atom_index]['neighbour']]['bond_type'])):
            if (more_neighbours[neighbours[atom_index]['neighbour']]['bond_type'][neighbour_index] == Chem.rdchem.BondType.DOUBLE and
                    smiles_mol.GetAtomWithIdx(more_neighbours[neighbours[atom_index]['neighbour']]['neighbour'][neighbour_index]).IsInRing() is False and
                    smiles_mol.GetAtomWithIdx(more_neighbours[neighbours[atom_index]['neighbour']]['neighbour'][neighbour_index]).GetSmarts() == 'C'):
                sss.assigned_atoms[more_neighbours[neighbours[atom_index]['neighbour']]['neighbour'][neighbour_index]] = 1
                tuple_atoms_index = tuple_atoms_index + (more_neighbours[neighbours[atom_index]['neighbour']]['neighbour'][neighbour_index],)
    if more_neighbours != {}:
        neighbour_dict_list.append(more_neighbours)
        if type(more_neighbours[neighbours[atom_index]['neighbour']]['neighbour']) == int:
            neighbour_list.append(more_neighbours[neighbours[atom_index]['neighbour']]['neighbour'])
        else:
            neighbour_list.extend(list(more_neighbours[neighbours[atom_index]['neighbour']]['neighbour']))
    neighbours = {}
    tuple_atoms_index = tuple(sorted(tuple_atoms_index))

    return(tuple_atoms_index, neighbours)

def neighbouring_double_bonds(sss,
                              smiles_mol,
                              neighbours,
                              atom_index,
                              tuple_atoms_index):
    ''' Function: neighbouring_double_bonds
    This function searches to see if there is a double or triple bond in the
    neighbours if so this atom now becomes part of the node and the bond gets
    appended to the bonds_to_reduce_list
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           neighbours - {} of atom indexes that are neighbours of the searched atom
                   along with there bond type
           atom_index - integer of the atom index number
           tuple_atoms_index - () tuple of the atom indexes that are within that node
    output: tuple_atoms_index - () tuple of the atom indexes that are within that node
            neighbour_dict_list - list of the dictionary of the neighbours
            neighbour_list - list of the neighbour atom indexes
    '''
    neighbour_dict_list = []
    neighbour_list = []

    if type(neighbours[atom_index]['neighbour']) == int:
        neighbours, tuple_atoms_index, bool_set = connection_of_interest(sss,
                                                                         smiles_mol,
                                                                         neighbours,
                                                                         atom_index,
                                                                         neighbours[atom_index]['neighbour'],
                                                                         neighbours[atom_index]['bond_type'],
                                                                         tuple_atoms_index,
                                                                         neighbour_dict_list,
                                                                         neighbour_list)
        if bool_set is False:
### new bit that have added to add S-OH together if S has four or six bonds
            if (smiles_mol.GetAtomWithIdx(atom_index).GetSmarts() == 'O' and
                smiles_mol.GetAtomWithIdx(neighbours[atom_index]['neighbour']).GetSmarts() == 'S' and
                smiles_mol.GetAtomWithIdx(neighbours[atom_index]['neighbour']).IsInRing() is False):
                tuple_atoms_index, neighbours = potential_sulphoxide(sss,
                                                                     smiles_mol,
                                                                     neighbours,
                                                                     atom_index,
                                                                     tuple_atoms_index,
                                                                     neighbour_dict_list,
                                                                     neighbour_list)

    elif type(neighbours[atom_index]['neighbour']) == tuple:
        for neighbour_index in range(len(neighbours[atom_index]['neighbour'])-1, -1, -1):
            neighbours, tuple_atoms_index, bool_set = connection_of_interest(sss,
                                                                             smiles_mol,
                                                                             neighbours,
                                                                             atom_index,
                                                                             neighbours[atom_index]['neighbour'][neighbour_index],
                                                                             neighbours[atom_index]['bond_type'][neighbour_index],
                                                                             tuple_atoms_index,
                                                                             neighbour_dict_list,
                                                                             neighbour_list)

    if neighbours != {}:
        neighbour_dict_list.append(neighbours)
        neighbour_list.extend((neighbours[atom_index]['neighbour'] if type(neighbours[atom_index]['neighbour']) is tuple else (neighbours[atom_index]['neighbour'],)))

    return(tuple_atoms_index,
           neighbour_dict_list,
           neighbour_list)

def set_basic_node(sss,
                   smiles_mol,
                   atoms_index,
                   node_type):
    ''' Function: set_basic_node
    This function create a node dictionary and assigns the values of the
    node dictionary
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           atoms_index - () tuple of the atom indexes that are contained within this node
           node_type - string of the type of node this wishes to be set as
    output: node_dictionary - {} of the node with values set
    '''
    node_dictionary = set_node_dictionary()
    all_neighbours = {}
    neighbours_index = []
    new_atoms_index = copy.deepcopy(atoms_index)
    connecting_bonds = []
    double_bond_neighbours = []
    ### creates a list of dictionaries of all the neighbours and bond types, filters out the blank
    for atom in atoms_index:
        sss.assigned_atoms[atom] = 1
        neighbours_dictionary, connecting_bonds_list, doubley_bonded_atoms_list = defining_neighbours(smiles_mol,
                                                                                                      atoms_index,
                                                                                                      atom)
        connecting_bonds.extend(connecting_bonds_list)
        double_bond_neighbours.extend(doubley_bonded_atoms_list)

        if bool(neighbours_dictionary) is True:
            if 'extra' in node_type:
                new_atoms_index, neighbour_dict_list, neighbour_list = neighbouring_double_bonds(sss,
                                                                                                 smiles_mol,
                                                                                                 neighbours_dictionary,
                                                                                                 atom,
                                                                                                 new_atoms_index)
                neighbours_index.extend(neighbour_list)
                for dictionary in neighbour_dict_list:
                    all_neighbours.update(dictionary)
                    if type(neighbours_dictionary[atom]['neighbour']) == int:
                        neighbours_index.append(neighbours_dictionary[atom]['neighbour'])
                    elif type(neighbours_dictionary[atom]['neighbour']) == tuple:
                        for neighbour in range(0, len(neighbours_dictionary[atom]['neighbour'])):
                            neighbours_index.append(neighbours_dictionary[atom]['neighbour'][neighbour])
            else:
                all_neighbours.update(neighbours_dictionary)
                if type(neighbours_dictionary[atom]['neighbour']) == int:
                    neighbours_index.append(neighbours_dictionary[atom]['neighbour'])
                elif type(neighbours_dictionary[atom]['neighbour']) == tuple:
                    for neighbour in range(0, len(neighbours_dictionary[atom]['neighbour'])):
                        neighbours_index.append(neighbours_dictionary[atom]['neighbour'][neighbour])

    atoms_index = new_atoms_index
    neighbours_index = list(set(neighbours_index))
    neighbours_index = [integer for integer in neighbours_index if integer not in atoms_index]
    node_dictionary['node_size'] = len(atoms_index)
    node_dictionary['atoms_index'] = atoms_index
    node_dictionary['neighbours'] = all_neighbours
    node_dictionary['neighbour_index'] = neighbours_index
    node_dictionary['nr_neighbours'] = len(neighbours_index)
    node_dictionary['atoms'] = [smiles_mol.GetAtomWithIdx(atom).GetSmarts() for atom in atoms_index]

    ### see whether any of the connecting bonds are actually contained within this node as then delete them from the list
    connecting_bonds_to_delete = []
    for atom1 in atoms_index:
        for atom2 in atoms_index:
            if atom1 != atom2 and smiles_mol.GetBondBetweenAtoms(atom1, atom2) != None:
                bond_idx = smiles_mol.GetBondBetweenAtoms(atom1, atom2)
                connecting_bonds_to_delete.append(bond_idx.GetIdx())  
    node_dictionary['connecting_bonds'] = [item for item in set(connecting_bonds) if item not in set(connecting_bonds_to_delete)]
    node_dictionary['double_bonds_connected'] = list(set(double_bond_neighbours))

    if 'extra' in node_type:
        node_type = node_type.replace('extra', '')
    if node_type == 'inert':
        node_string = 'acyclic inert'
    elif (node_type != 'acid' and
          node_type != 'base' and
          node_type != 'metal' and
          node_type != 'hydrophobic'):
        node_string = establishing_aromaticity(node_dictionary,
                                               node_type,
                                               [smiles_mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in atoms_index],
                                               [smiles_mol.GetAtomWithIdx(atom).IsInRing() for atom in atoms_index])
    else:
        node_string = node_type

    node_dictionary['node_type'] = node_string
    node_dictionary['node_code'] = node_codes[node_dictionary['node_type']]['code']
    node_dictionary['node_label'] = node_codes[node_dictionary['node_type']]['label']

    sss.dictionary_of_nodes[atoms_index] = node_dictionary
    return None

def addition_acyclic_inert(sss,
                           smiles_mol,
                           key,
                           neighbour_index):
    ''' Function: addition_acyclic_inert
    This function adds the key to any existing acyclic inert
    input: sss- class
           smiles_mol - mol object of the smiles that is being searched
           key - () tuple of the atom indexes
           neighbour_index - [] of the neighbours atom index
    output: None
    '''
    additions = []
    for unassigned in neighbour_index:
        for acyclic in [key for key, value in sss.dictionary_of_nodes.items() \
                        if unassigned in key and value['node_type'] == 'acyclic inert']:
            additions.append(acyclic)
    for old_key in set(additions):
        if type(old_key) == tuple:
            key = key + old_key
            del sss.dictionary_of_nodes[old_key]
        elif type(old_key) == list:
            for add_key in old_key:
                key = key + add_key
                del sss.dictionary_of_nodes[add_key]

    key = tuple(set(key))
    set_basic_node(sss,
                   smiles_mol,
                   key,
                   'inert')

    return None

def set_halogen_node(sss, smiles_mol, matches_list):
    ''' Function: set_halogen_node
    '''
    ### collect a list of the keys that are HBA as if next to a halogen then it is combined
    hba_list = [key for key, value in sss.dictionary_of_nodes.items() if
                'HBA' in value['node_type'] and 'acyclic' in value['node_type']]
    remove_list = []
    list_of_halogen_nodes = []
    halogens_neighbours = []

    for match in matches_list:
        for halogen in match:
            neighbours, connecting_bonds_list, doubley_bonded_atoms_list = defining_neighbours(smiles_mol,
                                                                                               [],
                                                                                               halogen)
            if neighbours != {}:
                hba_neighbours = [item for item in hba_list \
                                  if neighbours[halogen]['neighbour'] in item]
                if hba_neighbours != []:
                    new_key = hba_neighbours[0] + (halogen,)
                    if hba_neighbours[0] in sss.dictionary_of_nodes.keys():
                        sss.dictionary_of_nodes[new_key] = sss.dictionary_of_nodes[hba_neighbours[0]]
                        del sss.dictionary_of_nodes[hba_neighbours[0]]
                    else:
                        key_exists = [item for item in sss.dictionary_of_nodes.keys() \
                                      if set(item).issuperset(hba_neighbours[0]) is True]
                        if key_exists != []:
                            new_key = key_exists[0] + (halogen,)
                            sss.dictionary_of_nodes[new_key] = sss.dictionary_of_nodes[key_exists[0]]
                            del sss.dictionary_of_nodes[key_exists[0]]
                    sss.dictionary_of_nodes[new_key]['node_size'] = len(new_key)
                    sss.dictionary_of_nodes[new_key]['nr_atoms'] = len(new_key)
                    sss.dictionary_of_nodes[new_key]['atoms_index'] = new_key
                    sss.dictionary_of_nodes[new_key]['atoms'] = [smiles_mol.GetAtomWithIdx(item).GetSmarts() for item in new_key]
                    new_neighbours = {}
                    for atom in new_key:
                        atom_neighbour, connecting_bonds_list, doubley_bonded_atoms_list = defining_neighbours(smiles_mol,
                                                                                                               new_key,
                                                                                                               atom)
                        new_neighbours.update(atom_neighbour)
                    sss.dictionary_of_nodes[new_key]['neighbours'] = new_neighbours
                    remove_list.append(match)
                elif hba_neighbours == []:
                    halogens_neighbours.append(neighbours)
            elif neighbours == {}:
                list_of_halogen_nodes.append((halogen,))
            sss.halogens_atoms[halogen] = 1
            sss.assigned_atoms[halogen] = 1

    return(remove_list, list_of_halogen_nodes, halogens_neighbours)

def set_unassigned_node(args,
                        sss,
                        smiles_mol,
                        linker_nodes,
                        neighbour_index,
                        key):
    ''' Function: set_unassigned_node
    This function finds the current linker nodes within the mol and if these terminal linkers
    are next to then they combine or don't depending on the the argument of terminallinker
    input: args - arguments
           sss - class
           smiles_mol - mol object of the smiles that is being searched
           linker_nodes - [] of tuples of atom_indexes of the linker nodes
           neighbour_index - [] of atom indexes of neighbours to the node
           key - tuple of the atom_indexes of the node
    output: None
    '''
    all_matches = [item for neighbour in neighbour_index for item in \
                   linker_nodes if neighbour in item]

    for match in all_matches:
        if match in set(sss.dictionary_of_nodes.keys()):
            del sss.dictionary_of_nodes[match]

    all_matches = list(chain.from_iterable(all_matches))
    all_matches.extend(key)
    node_key = tuple(set(all_matches))

    if not args.terminallinker:
        atoms_list = [smiles_mol.GetAtomWithIdx(atom_index).GetSmarts() for atom_index in node_key]
        hybridization_list = [smiles_mol.GetAtomWithIdx(atom_index).GetHybridization() for atom_index in node_key]
        valence_list = [smiles_mol.GetAtomWithIdx(atom_index).GetImplicitValence() for atom_index in node_key]
        if (set(atoms_list) == {'C'} and
                [set(hybridization_list), set(valence_list)] == [{Chem.rdchem.HybridizationType.SP3}, {2, 3}] or
                [set(hybridization_list), set(valence_list)] == [{Chem.rdchem.HybridizationType.SP3}, {3}]):
            if len(neighbour_index) <= 1:
                for atom in node_key:
                    sss.assigned_atoms[atom] = -1
            else:
                addition_acyclic_inert(sss,
                                       smiles_mol,
                                       node_key,
                                       neighbour_index)
        else:
            addition_acyclic_inert(sss,
                                   smiles_mol,
                                   node_key,
                                   neighbour_index)
    elif args.terminallinker:
        addition_acyclic_inert(sss,
                               smiles_mol,
                               node_key,
                               neighbour_index)
    return None

def set_nodes(sss,
              smiles_mol,
              atoms_index,
              node_type):
    ''' Function: set_nodes
    This function looks to see if an node has already been defined then adds
    to if so if not sets the node
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           atoms_index - [] of tuples of atom indexes
           node_type - string of the type of node
    output: None
    '''
    for items in atoms_index:
        if items in sss.dictionary_of_nodes.keys():
            node_list = [sss.dictionary_of_nodes[items]['node_type'], node_type]
            new_node_type = ' '.join(node_list)
            if 'inert' in new_node_type and ('HBA' in new_node_type or 'HBD' in new_node_type or 'HBD HBA' in new_node_type) is True:
                new_node_type = new_node_type.replace('inert ', '')
            seen = set()
            new_node_type = [string for string in new_node_type.split(' ') if not (string in seen or seen.add(string))]
            new_node_type = ' '.join(new_node_type)
            sss.dictionary_of_nodes[items]['node_type'] = new_node_type
            sss.dictionary_of_nodes[items]['node_code'] = node_codes[new_node_type]['code']
            sss.dictionary_of_nodes[items]['node_label'] = node_codes[new_node_type]['label']

        elif [i for i in sss.dictionary_of_nodes.keys() if items[0] in i] != []:
            list_matches = [i for i in sss.dictionary_of_nodes.keys() if items[0] in i]
            if node_type not in sss.dictionary_of_nodes[list_matches[0]]['node_type']:
                node_list = [sss.dictionary_of_nodes[list_matches[0]]['node_type'], node_type]
                new_node_type = ' '.join(node_list)
                if 'inert' in new_node_type and ('HBA' in new_node_type or 'HBD' in new_node_type or 'HBD HBA' in new_node_type) is True:
                    new_node_type = new_node_type.replace('inert ', '')
                elif 'NHB' in new_node_type and ('HBA' in new_node_type or 'HBD' in new_node_type or 'HBD HBA' in new_node_type) is True:
                    new_node_type = new_node_type.replace('NHB ', '')
                seen = set()
                new_node_type = [string for string in new_node_type.split(' ') if not (string in seen or seen.add(string))]
                new_node_type = ' '.join(new_node_type)
            else:
                new_node_type = sss.dictionary_of_nodes[list_matches[0]]['node_type']
            sss.dictionary_of_nodes[list_matches[0]]['node_type'] = new_node_type
            sss.dictionary_of_nodes[list_matches[0]]['node_code'] = node_codes[new_node_type]['code']
            sss.dictionary_of_nodes[list_matches[0]]['node_label'] = node_codes[new_node_type]['label']

        else:
            set_basic_node(sss,
                           smiles_mol,
                           items,
                           ''.join(['extra', node_type]))

    return None

def reassigning_node(args,
                     sss,
                     smiles_mol,
                     old_key,
                     addition):
    ''' Function: reassigning_node
    This function adds atoms to an existing node
    input: args - arguments
           sss - class
           smiles_mol - mol object of the smiles that is being searched
           old_key - () tuple of the existing node
           addition - () tuple of the node that needs to add to it
    output: new_key - () tuple of the new node
    '''
    potential_key = [number for number in addition if number not in old_key]
    potential_key_list = []
    if tuple(potential_key) in sss.dictionary_of_nodes.keys():
        potential_key_list.append(sss.dictionary_of_nodes[tuple(potential_key)]['node_type'])
        del sss.dictionary_of_nodes[tuple(potential_key)]
    if type(addition) == int:
        sss.assigned_atoms[addition] = 1
        new_key = old_key + (addition,)
    elif type(addition) == tuple:
        new_key = old_key + addition
    new_key = tuple(set(new_key))
    if [k for k in sss.dictionary_of_nodes.keys() if set(k).issubset(set(addition)) is True] != []:
        for item in [k for k in sss.dictionary_of_nodes.keys() if set(k).issubset(set(addition)) is True]:
            if item != old_key:
                del sss.dictionary_of_nodes[item]

    sss.dictionary_of_nodes[new_key] = sss.dictionary_of_nodes[old_key]
    if addition not in sss.dictionary_of_nodes.keys():
        if not args.hydrophobic:
            set_basic_node(sss,
                           smiles_mol,
                           addition,
                           'inert')
        elif args.hydrophobic:
            atoms_list = [smiles_mol.GetAtomWithIdx(atom_idx).GetSmarts() for atom_idx in addition]
            if set(atoms_list) == {'C'}:
                set_basic_node(sss,
                               smiles_mol,
                               addition,
                               'inert')
            else:
                set_basic_node(sss,
                               smiles_mol,
                               addition,
                               'hydrophobic')

    node_type_list = [sss.dictionary_of_nodes[addition]['node_type'],
                      sss.dictionary_of_nodes[old_key]['node_type']]
    if potential_key_list != []:
        node_type_list.append(potential_key_list[0])
    if 'acyclic HBD HBA' in node_type_list:
        node_type = 'acyclic HBD HBA'
    elif 'acyclic HBD' in node_type_list and 'acyclic HBA' in node_type_list:
        node_type = 'acyclic HBD HBA'
    elif 'acyclic HBD' in node_type_list:
        node_type = 'acyclic HBD'
    elif 'acyclic HBA' in node_type_list:
        node_type = 'acyclic HBA'
    elif 'hydrophobic' in node_type_list:
        node_type = 'hydrophobic'
    elif 'acyclic positive' in node_type_list:
        node_type = 'acyclic positive'
    elif 'acyclic negative' in node_type_list:
        node_type = 'acyclic negative'
    elif 'acyclic inert' in node_type_list:
        node_type = 'acyclic inert'
    sss.dictionary_of_nodes[new_key]['node_type'] = node_type
    sss.dictionary_of_nodes[new_key]['node_label'] = node_codes[sss.dictionary_of_nodes[new_key]['node_type']]['label']
    sss.dictionary_of_nodes[new_key]['node_code'] = node_codes[sss.dictionary_of_nodes[new_key]['node_type']]['code']
    sss.dictionary_of_nodes[new_key]['atoms_index'] = new_key
    sss.dictionary_of_nodes[new_key]['atoms'] = [smiles_mol.GetAtomWithIdx(atom).GetSmarts() for atom in new_key]

    if type(addition) == tuple:
        if addition != new_key:
            del sss.dictionary_of_nodes[addition]

    sss.dictionary_of_nodes[new_key]['node_size'] = len(new_key)
    sss.dictionary_of_nodes[new_key]['nr_atoms'] = len(new_key)
    if old_key != new_key:
        del sss.dictionary_of_nodes[old_key]
    all_neighbours = {}
    all_neighbours_index = []
    connecting_bonds = []
    double_bond_neighbours = []

    for atom_index in new_key:
        neighbours, connecting_bonds_list, doubley_bonded_atoms_list = defining_neighbours(smiles_mol,
                                                                                           new_key,
                                                                                           atom_index)
        if neighbours != {}:
            all_neighbours.update(neighbours)
            all_neighbours_index.append(neighbours[atom_index]['neighbour'])
            connecting_bonds.extend(connecting_bonds_list)
            double_bond_neighbours.extend(doubley_bonded_atoms_list)
    sss.dictionary_of_nodes[new_key]['connecting_bonds'] = list(set(connecting_bonds))
    sss.dictionary_of_nodes[new_key]['double_bonds_connected'] = list(set(double_bond_neighbours))
    sss.dictionary_of_nodes[new_key]['neighbours'] = all_neighbours
    sss.dictionary_of_nodes[new_key]['neighbour_index'] = list(chain(*(i if isinstance(i, tuple) else (i,) for i in all_neighbours_index)))

    return new_key

def unknown_set_node(args, sss, smiles_mol, atom_index, other_atoms_ids):
    '''Function: unknown_set_node
    This function
    input: args - arguments
           sss - class
           smiles_mol - mol object of the smiles that is being searched
           atom_index - () tuple of the atom indexes
           other_atoms_ids - [] of the neighbour indexes
    output: None
    '''
    if not args.terminallinker:
        if atom_index not in sss.dictionary_of_nodes.keys():
            atoms_list = [smiles_mol.GetAtomWithIdx(index).GetSmarts() for index in atom_index]
            hybridization_list = [smiles_mol.GetAtomWithIdx(index).GetHybridization() \
                                  for index in atom_index]
            if (set(atoms_list) != {'C'} or
                    set(hybridization_list) != {Chem.rdchem.HybridizationType.SP3}):
                if args.hydrophobic:
                    hydrophobic_key = [key for other_atom in other_atoms_ids for \
                                       key, value in sss.dictionary_of_nodes.items() \
                                       if other_atom in key and (value['node_type'] == 'hydrophobic')]
                    if hydrophobic_key != []:
                        added_key = reassigning_node(args,
                                                     sss,
                                                     smiles_mol,
                                                     hydrophobic_key[0],
                                                     atom_index)
                        for item in range(1, len(hydrophobic_key)):
                            added_key = reassigning_node(args,
                                                         sss,
                                                         smiles_mol,
                                                         hydrophobic_key[item],
                                                         added_key)
                    else:
                        set_basic_node(sss,
                                       smiles_mol,
                                       atom_index,
                                       'hydrophobic')
                else:
                    set_basic_node(sss,
                                   smiles_mol,
                                   atom_index,
                                   'inert')

    if args.terminallinker:
        if atom_index not in sss.dictionary_of_nodes.keys():
            set_basic_node(sss,
                           smiles_mol,
                           atom_index,
                           'inert')

    return None

def set_combined_node(sss,
              smiles_mol,
              items_to_combine):
    ''' Function: set_combined_node
    This function combines the items into one node that are contained within
    the list items_to_combine
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           items_to_combine - [] of list of tuples of the node keys that need
                       to be combined each individual list needs to be combined
    output: None
    '''
    ### need to delete this key if it a subset within the list already
    new_items_to_combine = [m for i, m in enumerate(items_to_combine) if not \
                            any(set(m).issubset(set(n)) for n in \
                                (items_to_combine[:i] + items_to_combine[i+1:]))]
    items_to_delete = []
    for item in new_items_to_combine:
        list_of_atom_indexes = [key for key in item]
        items_to_delete.extend([key for key in item if tuple(set(key)) != tuple(set(chain.from_iterable(list_of_atom_indexes)))])
        node_type_list = []
        for key in item:
            node_type_string = sss.dictionary_of_nodes[key]['node_type'].split()
            node_type_string.pop(0)
            node_type_list.append(node_type_string)
        ### Try and except as if the list has a tuple that just contain one value
        ### then this fails i.e just an O or N
        try:
            list_of_atom_indexes = list(set([atom_index for list_ai in \
                                             list_of_atom_indexes for atom_index in list_ai]))
        except:
            list_of_atom_indexes_temp = []
            for list_ai in list_of_atom_indexes:
                if type(list_ai) == int:
                    list_of_atom_indexes_temp.append(list_ai)
                else:
                    list_of_atom_indexes_temp.extend([atom_index for atom_index in list_ai])
            list_of_atom_indexes = list(set(list_of_atom_indexes_temp))
        node_type = set(chain.from_iterable(node_type_list))
        if 'HBD' in node_type and 'HBA' in node_type:
            node_type = 'HBD HBA'
        elif 'HBD' in node_type and 'HBA' not in node_type:
            node_type = 'HBD'
        elif 'HBD' not in node_type and 'HBA' in node_type:
            node_type = 'HBA'
        else:
            node_type = list(node_type)[0]

        atom_in_ring_system = [smiles_mol.GetAtomWithIdx(match).IsInRing() for match in list_of_atom_indexes]

        if set(atom_in_ring_system) == {False}:
            set_basic_node(sss,
                           smiles_mol,
                           tuple(list_of_atom_indexes),
                           node_type)
        elif set(atom_in_ring_system) == {True, False} or set(atom_in_ring_system) == {False, True}:
            new_list_of_atom_indexes = [atom_item for atom_item in list_of_atom_indexes if smiles_mol.GetAtomWithIdx(atom_item).IsInRing() is False]
            if new_list_of_atom_indexes:
                set_basic_node(sss,
                               smiles_mol,
                               tuple(new_list_of_atom_indexes),
                               node_type)

    items_to_delete = list(set(items_to_delete))
    for delete in items_to_delete:
        del sss.dictionary_of_nodes[delete]

    return None
