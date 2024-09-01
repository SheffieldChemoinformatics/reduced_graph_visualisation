##############################################################################
# Reduced graph program
# Connecting Reduced Graph
#
#
# Jess Stacey
##############################################################################

from rdkit import Chem

def MolToMolBlock_WithAtomValues(mol, atom_values, mol_id):
    ''' Function: MolToMolBlock_WithAtomValues
    This function generates a mol block with the atom properties that are in
    a list in atom_values that are the number of atoms in that node and the
    smarts of that node
    input: mol - object mol
    output: mol_block - the mol block as a string
    '''
    try:
        mol_block = Chem.MolToMolBlock(mol).split("\n")
    except:
        params = Chem.SmilesParserParams()
        params.removeHs = True
        params.sanitize = False
        mol_block = Chem.MolToMolBlock(Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True), params)).split('\n')
    # Delete the "M  END" line.
    mol_block = mol_block[:-2]
    # Add appropriate "V" lines.
    for atom_idx in range(0, mol.GetNumAtoms()):
        prop_dict = mol.GetAtomWithIdx(atom_idx).GetPropsAsDict()
        mol_block.append(prop_dict['atomvalues'])
    mol_block.append("M  END")
    mol_block_sdf = Chem.MolFromMolBlock("\n".join(mol_block), sanitize=False, strictParsing=False)
    mol_block_sdf.SetProp('ID', str(mol_id))
    return mol_block_sdf

def generating_smarts(value, smiles_mol):
    ''' Function generating_smarts
    This function generates the smarts for that node by fragmenting the mol
    object and then deleting the atoms that are not part of that node
    input: value - {} of the nodes properties
           smiles_mol - mol object of the smiles that is being searched
    output: smarts - a string of the smarts
    '''
    if value['connecting_bonds'] != []:
        fragmented_mol = Chem.FragmentOnBonds(smiles_mol, value['connecting_bonds'])
        fragmented_mol = Chem.RWMol(fragmented_mol)
    else:
        fragmented_mol = Chem.RWMol(smiles_mol)
    atoms_to_delete = list(range(0, smiles_mol.GetNumAtoms()))
    atoms_to_delete = [atom for atom in atoms_to_delete if atom not in value['atoms_index']]
    atoms_to_delete.sort(reverse=True)
    for atom in atoms_to_delete:
        fragmented_mol.RemoveAtom(atom)
    smarts = Chem.MolToSmiles(fragmented_mol).split('.')
    smarts = [sm for sm in smarts if any(c.isalpha() for c in sm)]
    if len(smarts) == 1:
        smarts = smarts[0]
    else:
        print(smarts)
        print('Error Check The Subsmarts')
        smarts=smarts[0]
    
    return smarts

def generating_rg_from_scratch(args,
                               dictionary,
                               smiles_mol,
                               mol_id):
    '''Function: generating_rg_from_scratch
    This function does generates the RG from the node dictionary
    input: args - arguments
           dictionary - {} of all the nodes and properties to generate the RG
           smiles_mol - mol object of the smiles that is being searched
    output: new_mol - mol object of the reduced graph
            mol_block - string of the mol block
    '''
    new_mol = Chem.RWMol()
    atom_values = []
    for key, value in dictionary.items():
        new_atom = Chem.MolFromSmarts(value['node_label'])
        ### adding information into the atom so that know how many atoms
        ### and the other important information is retained so that it is possible
        ### to keep the relationship between the atoms and the node
        value['node_number'] = new_mol.AddAtom(new_atom.GetAtomWithIdx(0))
        
        smarts = generating_smarts(value, smiles_mol)
        atom_values.append("V  %3d %d %s %s" % (value['node_number'] + 1, value['node_size'], [str(v) for v in value['atoms_index']], smarts))
        new_mol.GetAtomWithIdx(value['node_number']).SetProp('atomvalues', "V  %3d %d %s %s" % (value['node_number'] + 1, value['node_size'], [str(v) for v in value['atoms_index']], smarts))
        for neighbour in value['neighbour_index']:
            tuples = [k for k, v in dictionary.items() if neighbour in k and v['node_number'] != None]
            for item in tuples:
                if new_mol.GetBondBetweenAtoms(dictionary[key]['node_number'], dictionary[item]['node_number']) == None and dictionary[key]['node_number'] != dictionary[item]['node_number']:
                    if args.doublebond and neighbour in value['double_bonds_connected']:
                        bond_number = new_mol.AddBond(dictionary[key]['node_number'], dictionary[item]['node_number'], Chem.rdchem.BondType.DOUBLE)
                        new_mol.GetBondWithIdx(bond_number - 1).SetProp('bond_type', 'nonring')
                    else:
                        new_mol.AddBond(dictionary[key]['node_number'], dictionary[item]['node_number'], Chem.rdchem.BondType.SINGLE)
        if value['fused_neighbours'] != None:
            for item in value['fused_neighbours']:
                rings = [k for k, v in dictionary.items() if v['nr_rings'] == item and v['node_number'] != None]
                for r in rings:
                    if new_mol.GetBondBetweenAtoms(dictionary[key]['node_number'], dictionary[r]['node_number']) == None and dictionary[key]['node_number'] != dictionary[r]['node_number']:
                        bond_number = new_mol.AddBond(dictionary[key]['node_number'], dictionary[r]['node_number'], Chem.rdchem.BondType.DOUBLE)
                        new_mol.GetBondWithIdx(bond_number - 1).SetProp('bond_type', 'ring')
    
    new_mol.SetProp('_Name', Chem.MolToSmiles(smiles_mol))
    ### Need someway of keeping all the order the same of the reduced graph 
    mol_block = MolToMolBlock_WithAtomValues(new_mol, atom_values, mol_id)

    return(new_mol, mol_block)
