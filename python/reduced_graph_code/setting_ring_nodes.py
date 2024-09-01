##############################################################################
# Reduced graph program
# Setting Ring Nodes
#
#
# Jess Stacey
##############################################################################

import copy
from itertools import chain, permutations

from setting_nodes import set_node_dictionary, defining_neighbours
from node_codes import node_codes_dict as node_codes


def establishing_number_of_rings_within_molecule(args,
                                                 smiles_mol,
                                                 list_of_rings):
    ''' Function: establishing_number_of_rings_within_molecule
    This function establishes the number of rings within the molecule
    by finding the number of ring atoms and ring bonds
    input: args - arguments
           smiles_mol - mol object of the smiles that is being searched
           list_of_rings - [] of tuples of the atom indexes within the rings
    output: number_of_rings - integer of the number of rings
    '''
    list_of_ring_bonds = [tup for tup in smiles_mol.GetRingInfo().BondRings() \
                          if len(tup) <= args.ring + 1]

    number_of_ring_atoms = len(set(chain.from_iterable(list_of_rings)))
    number_of_ring_bonds = len(set(chain.from_iterable(list_of_ring_bonds)))
    ### Calculating the number of rings that are present within the molecule,
    ### this equation is from Anantha Reddy paper c=b-a+1 b number of ring
    ### bonds a number of ring atoms
    number_of_rings = number_of_ring_bonds - number_of_ring_atoms + 1

    return number_of_rings

def function_fused_ring(atom_index,
                        neighbours_atoms_index,
                        list_of_rings):
    ''' Function: function_fused_rings
    This function find whether ring has any fused rings
    input: atom_index - integer of the atom index
           neighbours_atoms_index - integer or tuple of the atom indexes that
                                    are neighbours to i
           list_of_rings - [] of tuples of the atom indexes within the rings
    output: potential_fused_rings - [] of strings of the rings that are fused
    '''
    ### Finds whether the neighbours_atoms_index is in the list_of_rings and
    ### whether the atom_index is in the same ring if so it's a fused ring
    if type(neighbours_atoms_index) == int:
        tuple_neighbour = [element for element in list_of_rings \
                           if neighbours_atoms_index in element]
        tuple_ring = [item for item in tuple_neighbour if atom_index in item]

    elif type(neighbours_atoms_index) == tuple:
        tuple_ring = []
        for neighbour_index in neighbours_atoms_index:
            tuple_neighbour = [element for element in list_of_rings \
                               if neighbour_index in element]
            tuple_neighbour = [item for item in tuple_neighbour if atom_index in item]
            if tuple_neighbour != []:
                tuple_ring.append(tuple_neighbour)
        if len(tuple_ring) >= 1:
            tuple_ring = list(set(chain.from_iterable(tuple_ring)))

    potential_fused_rings = []
    ### creates a list of ring and index number based on positioning in the tuple_ring list
    if len(tuple_ring) >= 1:
        potential_fused_rings = [' '.join(['ring',
                                           str(list_of_rings.index(tuple_ring[tup]) + 1)]) \
            for tup in range(0, len(tuple_ring))]

    return potential_fused_rings

def set_ring_nodes(sss,
                   smiles_mol,
                   tuple_atoms_index,
                   aromatic_bool,
                   node_type,
                   list_of_rings):
    ''' Function: set_ring_nodes
    This function creates the ring node as a dictionary
    input: sss - class
           smiles_mol - mol object of the smiles that is being searched
           tuple_atoms_index - () tuple of atom indexes of the node
           aromatic_bool - string of where aromatic or aliphatic or acyclic
           node_type - string of the node type
           list_of_rings - [] of tuples of the atom indexes within the rings
    output: node_dictionary - {} dictionary of the node with values set
    '''
    node_dictionary = set_node_dictionary()
    node_type_string = []
    node_type_string.append(str(aromatic_bool))

    if node_type is None:
        node_type_string.append('NHB')
    elif node_type != None:
        node_type_string.append(node_type)

    node_type = ' '.join(node_type_string)
    ### Sets all the basics of the ring node
    node_dictionary['node_type'] = node_type
    node_dictionary['node_code'] = node_codes[node_type]['code']
    node_dictionary['node_label'] = node_codes[node_type]['label']
    node_dictionary['nr_rings'] = ' '.join(['ring',
                                            str(list_of_rings.index(tuple_atoms_index) + 1)])
    node_dictionary['node_size'] = len(tuple_atoms_index)
    node_dictionary['nr_atoms'] = len(tuple_atoms_index)
    node_dictionary['aromatic'] = aromatic_bool
    node_dictionary['atoms_index'] = tuple_atoms_index
    all_neighbours = {}
    neighbour_index = []
    potential_fused_neighbours = []
    connecting_bonds_list = []
    double_bond_neighbours = []

    for atom_index in tuple_atoms_index:
        neighbours, connecting_bonds, doubley_bonded_atoms_list = defining_neighbours(smiles_mol,
                                                                                      tuple_atoms_index,
                                                                                      atom_index)
        connecting_bonds_list.extend(connecting_bonds)
        double_bond_neighbours.extend(doubley_bonded_atoms_list)
        if bool(neighbours) is True:
            all_neighbours.update(neighbours)
            ###search to see if the neighbour is in ring if it is then change
            ### bond_type to double bond finding fused rings
            for value in neighbours.values():
                if type(value['neighbour']) == int:
                    neighbour_index.append(value['neighbour'])
                    value_atom = smiles_mol.GetAtomWithIdx(value['neighbour'])
                    if value_atom.IsInRing() is True:
                        potential_fused_neighbours.extend(function_fused_ring(atom_index,
                                                                              value['neighbour'],
                                                                              list_of_rings))
                elif type(value['neighbour']) == tuple:
                    for neighbour in value['neighbour']:
                        neighbour_index.append(neighbour)
                    neighbour_rings = [neighbour_id for neighbour_id in value['neighbour'] \
                                       if smiles_mol.GetAtomWithIdx(neighbour_id).IsInRing() \
                                       is True]
                    if len(neighbour_rings) >= 1:
                        potential_fused_neighbours.extend(function_fused_ring(atom_index,
                                                                              tuple(neighbour_rings),
                                                                              list_of_rings))
    potential_fused_neighbours = list(set(potential_fused_neighbours))

    node_dictionary['connecting_bonds'] = list(set(connecting_bonds_list))
    node_dictionary['double_bonds_connected'] = list(set(double_bond_neighbours))

    if potential_fused_neighbours != []:
        node_dictionary['fused_neighbours'] = potential_fused_neighbours
        node_dictionary['nr_fused_neighbours'] = len(potential_fused_neighbours)
    fused_ring_atoms = []

    for item in potential_fused_neighbours:
        string_ring = item.split(' ')
        ring_number = int(string_ring[1])
        fused_ring_atoms.append(list_of_rings[ring_number-1])
    fused_ring_atoms = list(chain.from_iterable(fused_ring_atoms))

    list_of_bonds_to_reduce = []
    for atom_index in set(all_neighbours.keys()):
        if (type(all_neighbours[atom_index]['neighbour']) == int and
                all_neighbours[atom_index]['neighbour'] in set(fused_ring_atoms)):
            neighbour_index.remove(all_neighbours[atom_index]['neighbour'])
            list_of_bonds_to_reduce.append([atom_index,
                                            all_neighbours[atom_index]['neighbour'],
                                            'remove'])
            del all_neighbours[atom_index]

        elif type(all_neighbours[atom_index]['neighbour']) == tuple:
            all_neighbours_copy = copy.deepcopy(all_neighbours[atom_index])
            for atom in all_neighbours_copy['neighbour']:
                if atom in set(fused_ring_atoms):
                    neighbour_index.remove(atom)
                    index_to_delete = all_neighbours[atom_index]['neighbour'].index(atom)
                    new_neighbour_list = list(all_neighbours[atom_index]['neighbour'])
                    new_bond_list = list(all_neighbours[atom_index]['bond_type'])
                    list_of_bonds_to_reduce.append([atom_index,
                                                    atom,
                                                    'remove'])
                    del new_neighbour_list[index_to_delete]
                    del new_bond_list[index_to_delete]
                    if new_neighbour_list == []:
                        del all_neighbours[atom_index]
                    else:
                        all_neighbours[atom_index]['neighbour'] = tuple(new_neighbour_list)
                        all_neighbours[atom_index]['bond_type'] = tuple(new_bond_list)

    node_dictionary['atoms'] = [smiles_mol.GetAtomWithIdx(atom_index).GetSmarts() \
                   for atom_index in tuple_atoms_index]
    node_dictionary['neighbours'] = all_neighbours
    node_dictionary['neighbour_index'] = neighbour_index
    node_dictionary['nr_neighbours'] = len(neighbour_index)
    sss.dictionary_of_nodes[tuple_atoms_index] = node_dictionary

    return None

def smallest_set_rings(args,
                       sss,
                       smiles_mol,
                       list_of_rings):
    ''' Function: smallest_set_rings
    This function establishes the smallest set of rings. This is done by
    ordering the list into ascending size order from here the rings that
    contain something interesting, i.e. heterocycles, in this case this
    means something such as a HBA, HBD or HBD/HBA then these are prioritised
    over just rings that contain carbon. From here the tuple of atom indexes
    is searched to see if all the atoms within this tuple are already defined
    as being in a ring if it is then this ring is removed
    input: args - arguments of the program
           sss - class
           smiles_mol - mol object of the smiles that is being searched
           list_of_rings - [] of tuples of the atom indexes of the rings
    output: list_of_rings - [] of tuples of the atom indexes of the rings
    '''
    original_list_of_rings = copy.deepcopy(list_of_rings)
    list_of_rings = [tuple(sorted(item)) for item in list_of_rings]
    list_of_rings.sort()
    list_of_rings = sorted(list_of_rings, key=lambda ring: len(ring))
    list_of_rings = [tup for tup in list_of_rings if len(tup) <= args.ring]
    number_of_rings = establishing_number_of_rings_within_molecule(args,
                                                                   smiles_mol,
                                                                   list_of_rings)

    if number_of_rings != len(list_of_rings):
        aliphatic_hba_hbd = [key for key, value in sss.dictionary_of_nodes.items() \
                             if 'aliphatic' in value['node_type']]
        aliphatic_hba_hbd = list(chain.from_iterable(aliphatic_hba_hbd))
        dictionary_ring_length = {}

        for item in list_of_rings:
            length_tup = len(item)
            if length_tup in set(dictionary_ring_length.keys()):
                dictionary_ring_length[length_tup].append(item)
            else:
                dictionary_ring_length[length_tup] = [item]
        list_of_rings_new = []

        for ring_atom_indexes in dictionary_ring_length.values():
            ### puts the rings with HBA and HBD groups at the beginning of the list of tuples
            if aliphatic_hba_hbd != []:
                for ring_atoms in ring_atom_indexes:
                    check_list = []
                    if len(aliphatic_hba_hbd) == 1:
                        check_list.append(aliphatic_hba_hbd[0] in ring_atoms)
                    elif len(aliphatic_hba_hbd) > 1:
                        check_list.extend([atom in ring_atoms for atom in aliphatic_hba_hbd])
                    if set(check_list) == {False}:
                        ring_atom_indexes.append(ring_atom_indexes.pop(ring_atom_indexes.index(ring_atoms)))
            new_value = copy.deepcopy(ring_atom_indexes)
            atoms_to_remain = []
            for item in ring_atom_indexes:
                if item[0] in set(atoms_to_remain) and aliphatic_hba_hbd == []:
                    new_value.append(new_value.pop(new_value.index(item)))
                elif item[0] not in set(atoms_to_remain):
                    atoms_to_remain.append(item[0])
            list_of_rings_new.append(new_value)

        list_of_rings_new = list(chain.from_iterable(list_of_rings_new))
        list_of_rings = []
        list_of_bonds = []

        for ring_atoms in list_of_rings_new:
            enumeration_number = [index for index, ring_atoms_indexes in \
                                  enumerate(original_list_of_rings, 0) if \
                                  tuple(sorted(ring_atoms_indexes)) == ring_atoms]
            bonds = smiles_mol.GetRingInfo().BondRings()[enumeration_number[0]]
            bond_ring_appearance = []
            for bond in bonds:
                if bond in set(list_of_bonds):
                    bond_ring_appearance.append(True)
                elif bond not in set(list_of_bonds):
                    bond_ring_appearance.append(False)
                    list_of_bonds.append(bond)
                if set(bond_ring_appearance) == {False} or \
                set(bond_ring_appearance) == {False, True} \
                and ring_atoms not in list_of_rings:
                    list_of_rings.append(ring_atoms)

    return list(set(list_of_rings))

def establishing_ring_nodes(args,
                            sss,
                            smiles_mol):
    ''' Function: establishing_ring_nodes
    This function establishes the ring nodes within the smiles_mol
    input: args - arguments of the program
           sss - class
           smiles_mol - mol object of the smiles that is being searched
    output: None
    '''
    if args.verbose:
        print('Establishing Ring Nodes')
    # may need to remove the number from the atom when a ring is formed
    ring_info = smiles_mol.GetRingInfo()
    ### list of tuples containing the index of the atoms that are contained within each ring
    list_of_rings = ring_info.AtomRings()
    
    if len(list_of_rings) > 0:
        aliphatic_aromatic_atoms = [key for key, value in sss.dictionary_of_nodes.items() \
                                    if ('aliphatic' in value['node_type'] or \
                                        'aromatic' in value['node_type'])]
        pos_neg_atom = [key for key, value in sss.dictionary_of_nodes.items() \
                        if ('positive' in value['node_type'] or 'negative' in value['node_type'])] 
        acyclic_atoms = [key for key, value in sss.dictionary_of_nodes.items() \
                                    if 'acyclic' in value['node_type']]
    
        if len(aliphatic_aromatic_atoms) == 0:
            max_tuple_length_atoms = 0
        else:
            max_tuple_length_atoms = max(map(len, aliphatic_aromatic_atoms))

        list_of_rings = smallest_set_rings(args,
                                           sss,
                                           smiles_mol,
                                           list_of_rings)
        dict_items_to_del = [(atom_index,) for ring in list_of_rings for atom_index
                             in ring if (atom_index,) in set(aliphatic_aromatic_atoms)]
        
        already_defined_rings = {}
        for ring in list_of_rings:
            for item in pos_neg_atom:
                if tuple(i for i in ring if i in item) == ring:
                    already_defined_rings[ring] = item
            for item in acyclic_atoms:
                if tuple(i for i in ring if i in item) == ring:
                    already_defined_rings[ring] = item
        
        for item in already_defined_rings.keys():
            list_of_rings.remove(item)
            list_of_rings.append(already_defined_rings[item])
        list_of_rings = list(set(list_of_rings))

        for ring in list_of_rings:
            aromatic_list = [smiles_mol.GetAtomWithIdx(atom_index).GetIsAromatic() \
                             for atom_index in ring]
            node_type = []
            ring_atom_index_that_already_defined = [atom for atom in ring if atom in set(chain.from_iterable(aliphatic_aromatic_atoms))]
            list_of_tuples_of_aliphatic_aromatic_atoms = [tuple(permutations(ring_atom_index_that_already_defined, 1))]
            for i in range(1, max_tuple_length_atoms+1):
                list_of_tuples_of_aliphatic_aromatic_atoms.append(tuple(permutations(ring_atom_index_that_already_defined, i)))
            list_of_tuples_of_aliphatic_aromatic_atoms = [item for sublist in list_of_tuples_of_aliphatic_aromatic_atoms for item in sublist]
            for atom_index in list_of_tuples_of_aliphatic_aromatic_atoms:
                if atom_index in set(aliphatic_aromatic_atoms):
                    dict_items_to_del.append(atom_index)
                    potential_node_type = sss.dictionary_of_nodes[atom_index]['node_type'].split()
                    potential_node_type.pop(0)
                    potential_node_type = ' '.join(potential_node_type)
                    node_type.append(potential_node_type)
                    
            aromatic_set = set(aromatic_list)
            ### see whether the aromatic set contains just true or false or
            ### both to set the aromatic_boolean
            if len(aromatic_set) == 1 and list(aromatic_set)[0] is True:
                aromatic_bool = 'aromatic'
            else:
                aromatic_bool = 'aliphatic'
            
            if ring in sss.dictionary_of_nodes.keys():
                potential_node_type = sss.dictionary_of_nodes[ring]['node_type'].split()
                potential_node_type.pop(0)
                potential_node_type = ' '.join(potential_node_type)
                node_type.append(potential_node_type)
                
                separated_rings_in_def = [subring for subring, node in already_defined_rings.items() if node == ring]
                aromaticity_list = []
                for itemy in separated_rings_in_def:
                    new_aromatic_list = [smiles_mol.GetAtomWithIdx(atom_index).GetIsAromatic() \
                                         for atom_index in itemy]
                    aromaticity_list.extend(set(new_aromatic_list))
                if True in aromaticity_list:
                    aromatic_bool = 'aromatic'
                else:
                    aromatic_bool = 'aliphatic'

            if len(set(node_type)) == 1:
                node_type = node_type[0]
            elif len(set(node_type)) > 1:
                node_type = 'HBD HBA'
            else:
                node_type = None
            ### creates node_dictionary for the ring, takes in the number
            ### of atoms within rings and whether the ring is aromatic or not
            set_ring_nodes(sss,
                           smiles_mol,
                           ring,
                           aromatic_bool,
                           node_type,
                           list_of_rings)
        ### this deletes the value from the dictionary_of_nodes if i.e an n is
        ### in ring delete separate n entry
        for delete_key in dict_items_to_del:
            if delete_key in list(sss.dictionary_of_nodes.keys()):
                del sss.dictionary_of_nodes[delete_key]
    sss.assign_ring_nodes(list_of_rings)

    return None
