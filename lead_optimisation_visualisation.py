#!/usr/bin/env python
# -*- coding: utf-8 -*-
##########################################################################################
# lead_optimisation_visulisation.py
#
# 
# Jess Stacey
#       August 2018
#
##########################################################################################

from __future__ import print_function
try:
    from urllib.parse import unquote
except ImportError:
    from urllib import unquote

from flask import Flask, render_template, request, send_file

from werkzeug.utils import secure_filename

from flask_script import Manager
app = Flask(__name__)
app.debug = False
manager=Manager(app)

import pandas as pd
import json
import os
import re
import ast
import time
import zlib
from rdkit import Chem
import base64
from rdkit.six import BytesIO
from rdkit.Chem import Draw


class static_variables:
    ''' Class: static_variables
    This class contains all the static variables
    '''
    def __init__(self):
        self.program_name = 'Lead Optimisation Tool'
        self.sheffield_logo = 'static/sheffield.png'
        self.gsk_logo = 'static/gsk.png'

class core_information:
    ''' Class: core_information
    This class contains all the dictionary variables of the core
    '''
    def __init__(self):
        self.core_dictionary = {}
        self.additional_dictionary = {}
        self.table_dictionary = {}

class running_new_dataset_class:
    ''' Class: running_new_dataset_class
    This class contains all the parameters for running the new dataset
    '''
    def __init__(self):
        self.initial_file = None
        self.minsim = 0.5
        self.minnodes = 4
        self.fingerprintbitlength = 2048
        self.rounddata = False
        self.roundnames = []
        self.coordinates_file = None
        self.node_information_file = None
        self.core_analysis_file = None
    
def generate_smarts_image(SMARTS):
    ''' Function: generate_smarts_image
    This function takes the inputted SMARTS and saves it into an image which
    is then transformed to a base64 string
    Input: SMARTS - str of SMARTS
    Output: data_uri - image of the SMARTS as a base64 string
    '''
    params = Chem.SmilesParserParams()
    params.removeHs = True
    params.sanitize = False
    rg_mol = Chem.MolFromSmiles(SMARTS, params)
    image_data = BytesIO()
    img = Draw.MolToImage(rg_mol, size=(150, 150), kekulize=False, wedgebonds=True, fitImage=False)
    img.save(image_data, format='PNG')
    image_data.seek(0)
    data_uri = base64.b64encode(image_data.read()).decode('ascii')

    return data_uri

def generate_image(SMARTS):
    ''' Function: generate_image
    This function take the inputted SMART and saves it into a image which is
    then transformed to a base64 string
    Input: SMARTS - str of SMARTS
    Output: data_uri - image of the SMARTS as a base64 string
    '''
    if SMARTS == 'nan' or pd.isnull(SMARTS):
        return ""
    else:
        rg_mol = Chem.MolFromSmiles(SMARTS)
        if rg_mol is None:
            rg_mol = Chem.MolFromSmarts(SMARTS)
        if rg_mol is None:
            params = Chem.SmilesParserParams()
            params.removeHs = True
            params.sanitize = False
            rg_mol = Chem.MolFromSmiles(SMARTS, params)
        image_data = BytesIO()
        img = Draw.MolToImage(rg_mol, size=(150, 150), kekulize=False, wedgebonds=True, fitImage=False)
        img.save(image_data, format='PNG')
        image_data.seek(0)
        data_uri = base64.b64encode(image_data.read()).decode('ascii')

        return data_uri
    
def generate_table_image(SMARTS):
    ''' Function: generate_image
    This function take the inputted SMART and saves it into a image which is
    then transformed to a base64 string
    Input: SMARTS - str of SMARTS
    Output: data_uri - image of the SMARTS as a base64 string
    '''
    if SMARTS == 'nan' or pd.isnull(SMARTS):
        return ""
    else:
        rg_mol = Chem.MolFromSmiles(SMARTS)
        if rg_mol is None:
            rg_mol = Chem.MolFromSmarts(SMARTS)
        if rg_mol is None:
            params = Chem.SmilesParserParams()
            params.removeHs = True
            params.sanitize = False
            rg_mol = Chem.MolFromSmiles(SMARTS, params)
        image_data = BytesIO()
        img = Draw.MolToImage(rg_mol, size=(500, 500), fitImage=True)
        img.save(image_data, format='PNG')
        image_data.seek(0)
        data_uri = base64.b64encode(image_data.read()).decode('ascii')

        return data_uri
    
def generate_highlighted_table_image(SMILES, core_atom_indexes, rg=False):
    ''' Function: generate_highlighted_table_image
    This function take the inputted SMILES and highlight the core atom indexes
    and saves it into a image which is then transformed to a base64 string
    Input: SMILES - str of SMILES of molecule or mol block of RG
           core_atom_indexes - list of atom indexes
           rg - bool as to whether the SMILES is a 
    Output: data_uri - image of the SMARTS as a base64 string
    '''
    if rg:
        mol_block = zlib.decompress(base64.b64decode(SMILES))
        smiles_mol = Chem.MolFromMolBlock(mol_block, sanitize=False, strictParsing=False, removeHs=True)
        for atom in smiles_mol.GetAtoms():
            atom.SetNumExplicitHs(0)
    else:
        smiles_mol = Chem.MolFromSmiles(SMILES)
        if smiles_mol == None:
            params = Chem.SmilesParserParams()
            params.removeHs = True
            params.sanitize = False
            smiles_mol = Chem.MolFromSmiles(SMILES, params)
    
    image_data = BytesIO()
    img = Draw.MolToImage(smiles_mol, size=(500, 500), fitImage=True, highlightAtoms=core_atom_indexes)
    img.save(image_data, format='PNG')
    image_data.seek(0)
    data_uri = base64.b64encode(image_data.read()).decode('ascii')

    return data_uri

@app.route('/')
def set_up():
    ''' Function: set_up
    Return a simple HTML file
    Input: None
    Output: render template function
    '''
    menu_options = dataset_json["datasets"]
    
    return render_template('layout.html', datasets=menu_options)
    
def cleaning_function(smarts):
    ''' Function: cleaning_function
    Cleans the SMARTS by removing any numbers on wild atoms and canonicalising them
    Input: smarts - str of SMARTS
    Output: SMARTS string
    '''
    ### remove numbered wild atoms to just wild atoms
    smarts_cleaned = re.sub(r'\[\d+\*\]', '[*]', smarts)
    ### canonicalise the smarts
    mol = Chem.MolFromSmarts(smarts_cleaned)
    canonical_smarts = Chem.MolToSmiles(mol, isomericSmiles=True)
    
    return canonical_smarts

@app.route('/<dataset>/<round_investigating>/multi_core_rg_generator/', methods=['GET'])
def multi_core_rg_generator(dataset, round_investigating):
    ''' Function: multi_core_rg_generator
    This function generators all the node and edges for all the cores in this
    dataset
    Input: dataset - str of the name of the dataset
           round_investigating - round that is being investigated
    Output: json of the cores data
    '''
    if dataset == 'new_dataset':
        key_information_data = pd.read_csv(rnd.coordinates_file, sep='\t', header=0)
        core_information_data = pd.read_csv(rnd.node_information_file, sep='\t', header=0)
        core_breakdown_data = pd.read_csv(rnd.core_analysis_file, sep='\t', header=0)
    else:
        key_information_data = pd.read_csv('datasets/{0}_coordinates.txt'.format(dataset), sep='\t', header=0)
        core_information_data = pd.read_csv('datasets/{0}_node_information.txt'.format(dataset), sep='\t', header=0)
        core_breakdown_data = pd.read_csv('datasets/{0}_core_analysis.txt'.format(dataset), sep='\t', header=0)
        
    previous_round_data = {}
    list_of_new_cores = []
    if round_investigating != 'norounds':
        ### find the index of the round 
        dataset_rounds = sorted(set(core_breakdown_data['Round'].tolist()), key=core_breakdown_data['Round'].tolist().index)
        index_of_round = dataset_rounds.index(round_investigating)
        previous_cores = []
        if index_of_round != 0:
            core_data_before = core_breakdown_data[core_breakdown_data['Round'] == dataset_rounds[index_of_round-1]]
            core_data_before['functional_groups_breakdown'] = core_data_before['functional_groups_breakdown'].apply(ast.literal_eval)
            for idx, row in core_data_before.iterrows():
                node_dict = {}
                for fg in row['functional_groups_breakdown']:
                    node_dict[fg['SMILES']] = fg['Occurrences']
                previous_round_data[(row['name'], row['core'])] = node_dict
            ### note the cores that are present in this round
            previous_cores = set(core_data_before['core'].tolist())
        
        ##filter
        key_information_data = key_information_data[key_information_data['Round'] == round_investigating]
        core_information_data = core_information_data[core_information_data['Round'] == round_investigating]
        core_breakdown_data = core_breakdown_data[core_breakdown_data['Round'] == round_investigating]
        ### note the cores that are new in this round
        current_cores = set(core_breakdown_data['core'].tolist())
        list_of_new_cores = [val for val in current_cores if val not in previous_cores]
            
    core_information_data['core_numbered'] = core_information_data['core_numbered'].apply(ast.literal_eval)
    core_information_data['core_info'] = core_information_data['core_info'].apply(ast.literal_eval)
    core_information_data['core_edges'] = core_information_data['core_edges'].apply(ast.literal_eval)
    core_information_data['core_node_indexes'] = core_information_data['core_node_indexes'].apply(ast.literal_eval)
    core_information_data['core_atom_indexes'] = core_information_data['core_atom_indexes'].apply(ast.literal_eval)
    core_breakdown_data['functional_groups_list'] = core_breakdown_data['functional_groups_list'].apply(ast.literal_eval)
    core_breakdown_data['functional_groups_breakdown'] = core_breakdown_data['functional_groups_breakdown'].apply(ast.literal_eval)
    
    edges = []
    core_info = {}
    table_info = {}
    # number_of_nodes = 0
    number_of_molecules = {}

    ### filter for each core
    node_list = []
    core_set_up_list= []
    
    ###Sets up the node dictionary so that the nodes can be generated
    ### Also establishes the node_dict of the nodes that are in each core
    image_created = {}
    for index, core_row in core_breakdown_data.iterrows():
        for item in core_row['functional_groups_breakdown']:
            ### clean smarts as well!
            if item['Image'] in image_created.keys():
                item['Image'] = image_created[item['Image']]
            else:
                image_created[item['Image']] = generate_smarts_image(cleaning_function(item['Image']))
                item['Image'] = image_created[item['Image']]
            ### add entry of how many extra from last round if round data
            if round_investigating != 'norounds':
                if previous_round_data:
                    if (core_row['name'], core_row['core']) in previous_round_data.keys():
                        list_of_previous_groups = previous_round_data[(core_row['name'], core_row['core'])]
                        if item['SMILES'] in list_of_previous_groups:
                            item['AddedExamples'] = item['Occurrences'] - list_of_previous_groups[item['SMILES']]
                            item['PreviousExamples'] = 'Old'
                        else:
                            item['AddedExamples'] = item['Occurrences'] - 0
                            item['PreviousExamples'] = 'NEWFRAG'
                    else:
                        item['PreviousExamples'] = 'NEWCORE'
        
        if core_row['core'] in core_info.keys():
            core_info[core_row['core']].append(core_row['name'])
        else:
            core_info[core_row['core']] = [core_row['name']]

        node_list.append(core_row.to_dict())
        if core_row['core'] not in core_set_up_list:
            core_set_up_list.append(core_row['core'])
        number_of_molecules[core_row['core']] = core_row['number_of_molecules']

    ### For each core has to find the corresponding molecules
    for core in core_set_up_list:
        print('Processing core: ', core)
        molecular_information = []
        table_info[core] = []
        core_nodes_table_info = []
        core_info_data_filtered = core_information_data[core_information_data.core == core]
        ### This generates the dictionary that allows the edges to be in the correct places
        for dictionary in core_info_data_filtered.loc[core_info_data_filtered.index[0], 'core_edges']:
            edges.append({'source': dictionary['source'],# + number_of_nodes,
                          'target': dictionary['target'],# + number_of_nodes,
                          'weight': dictionary['weight'],
                          'core': core})
        
        ### need a way of appenidng the number of nodes that are present!!! 
        # number_of_nodes = number_of_nodes + len(core_info[core])
        
        core_smarts_examined = []
        for index, (idx, row) in enumerate(core_info_data_filtered.iterrows()):
            mol_core_table_info = {}
            for key, value in row['core_info'].items():
                ###clean smarts!!!
                mol_core_table_info[row['core_numbered'][key]] = generate_smarts_image(cleaning_function(value))
            mol_core_table_info['combined'] = generate_smarts_image(cleaning_function(row['core_smarts']))
            mol_core_table_info['number of examples'] = core_info_data_filtered['core_smarts'].tolist().count(row['core_smarts'])
            core_smarts_examined.append(row['core_smarts'])
            core_nodes_table_info.append(mol_core_table_info)
            
            ### Sets up the molecular information when just exploring the one core
            ### generate molecular_information images in core_rg_generator
            molecular_information.append({'SMILES': row['SMILES'],
                                          'Image': generate_table_image(row['SMILES']),
                                          'ID': row['ID'],
                                          'rg_SDF': key_information_data[key_information_data['SMILES'] == row['SMILES']]['SDF'].tolist()[0],
                                          'pIC50': row['pIC50'],
                                          'Reduced Graph': row['RG'],
                                          'Core Node Indexes' : row['core_node_indexes'],
                                          'Core': row['core_smarts'],
                                          'Core Atom Indexes' : row['core_atom_indexes'],
                                          'Core Smarts': row['core_smarts']})
            
        ### Set the list of dictionaries
        core_nodes_table_info_set = []
        for n in core_nodes_table_info:
            if n not in core_nodes_table_info_set:
                core_nodes_table_info_set.append(n)
        for core_node_info_row in core_nodes_table_info_set:
            core_node_info_row['number of examples'] = core_nodes_table_info.count(core_node_info_row)
        ### Now order based upon the number of examples
        core_nodes_table_info = sorted(core_nodes_table_info_set, key=lambda d: d['number of examples'], reverse=True)            
 
        table_info[core].append(core_nodes_table_info)

        ci.core_dictionary[core] = {'nodes': [node for node in node_list if node['core'] == core],
                                    'edges': core_info_data_filtered.loc[core_info_data_filtered.index[0], 'core_edges'],
                                    'molecular_information': molecular_information}#

    ### iterate through to generate smarts_image
    ci.table_dictionary = table_info
    core_dictionary = {}
    core_dictionary['nodes'] = node_list
    core_dictionary['edges'] = edges
    core_dictionary['cores'] = core_info
    core_dictionary['table'] = table_info
    core_dictionary['number_of_molecules'] = number_of_molecules
    core_dictionary['text'] = [len(set(core_information_data['SMILES'].tolist())),
                               len(set(key_information_data['RG'].tolist())),
                               len(set(core_information_data['core'].tolist()))]
    core_dictionary['new_cores'] = list_of_new_cores

    return json.dumps(core_dictionary)

@app.route('/<dataset>/core_rg_generator/<core>/', methods=['GET'])
def core_rg_generator(dataset, core):
    ''' Function: core_rg_generator
    This function retrieves all the node and edges data for <core>
    Input: dataset - str of the name of the dataset
           core - str of core that has been clicked for further investigation
    Output: json of this cores data
    '''
    core_data = ci.core_dictionary[core]
    core_data_new = {}
    for key, value in core_data.items():
        if key != 'molecular_information':
            core_data_new[key] = value
        else:
            fg_breakdown = []
            for cell_data in value:
                new_dictionary = {'SMILES': cell_data['SMILES'],
                                  'Image': cell_data['Image'],
                                  'ID': cell_data['ID'],
                                  'pIC50': cell_data['pIC50'],
                                  'Reduced Graph': cell_data['Reduced Graph'],
                                  'Reduced Graph Image': generate_highlighted_table_image(cell_data['rg_SDF'], cell_data['Core Node Indexes'], rg=True),
                                  'Core': generate_highlighted_table_image(cell_data['SMILES'], cell_data['Core Atom Indexes'], rg=False),
                                  'Core Smarts': cell_data['Core Smarts']}
                fg_breakdown.append(new_dictionary)
            core_data_new['molecular_information'] = fg_breakdown  

    return json.dumps(core_data_new)

@app.route('/core_comparison/<core>/', methods=['GET'])
def core_comparison(core):
    ''' Function: core_comparison
    This function
    Input: core - str of the core
    Output: json of the core data
    '''
    core_data = ci.core_dictionary[core]
    core_data['table'] = {core: ci.table_dictionary[core]}

    return json.dumps(core_data)

@app.route('/core_rg_generator/<int:clusterID>/additional_groups/', methods=['GET'])
def core_rg_generator_additional_groups(clusterID):
    ''' Function: cluster_rg_generator_additional_groups
    This function return the additional node functional groups in a dictionary
    Input: clusterID - int of the cluster that is being looked at
    Output: json of the additional node data
    '''
    return json.dumps(ci.additional_dictionary)

@app.route('/<dataset>/rounds/', methods=['GET'])
def recognise_rounds(dataset):
    ''' Function recognise_rounds
    Deciphers the round data for a chosen dataset
    Input: dataset - str of name of dataset
    Output: json of all the round data for the chosen dataset
    '''
    if dataset == 'new_dataset':
        core_breakdown_data = pd.read_csv(rnd.core_analysis_file, sep='\t', header=0)
    else:
        core_breakdown_data = pd.read_csv('datasets/{0}_core_analysis.txt'.format(dataset), sep='\t', header=0)
        
        
    if 'Round' in core_breakdown_data.columns:
        ## would ideally like it ordered
        rounds = list(sorted(set(core_breakdown_data['Round'].tolist()),
                        key=core_breakdown_data['Round'].tolist().index))
        
        all_rounds_df = pd.DataFrame()
        
        for _round in rounds:
            ##filter
            core_round_df = core_breakdown_data[core_breakdown_data['Round'] == _round]
            core_round_df2 = core_round_df[['core','number_of_molecules']]
            core_round_df2 = core_round_df2.drop_duplicates()
            
            ## sort
            core_round_df2.sort_values(by=['number_of_molecules'], inplace=True, ascending=False)
            core_round_df2['row_index'] = core_round_df2['core']
            core_round_df2.set_index('row_index', inplace=True)

            ## Make multiindex - rename columns
            core_round_df2.rename(columns = {'core':'{0},core'.format(_round),
                                           'number_of_molecules':'{0},number_of_molecules'.format(_round)}, inplace = True)
            all_rounds_df = pd.concat([all_rounds_df, core_round_df2], axis=1, sort=False)
            if len(all_rounds_df.columns) > 2:
                ## make column that is how many introduced
                all_rounds_df['{0},difference'.format(_round)] = all_rounds_df.iloc[:, -1].fillna(0) - all_rounds_df.iloc[:, -3].fillna(0)
                all_rounds_df = pd.concat([all_rounds_df, core_round_df2], axis=1, sort=False)

        ## Make multiindex
        multicolums = all_rounds_df.columns.str.split(',', expand=True).values
        all_rounds_df.columns = pd.MultiIndex.from_tuples([x for x in multicolums])
        
        all_rounds_df = all_rounds_df.fillna(0)
        
        core_progression_dict = {}
        for idx, row in all_rounds_df.iterrows():
            row_data = {}
            for r in row.iteritems():
                if r[0][0] in row_data.keys():
                    row_data[r[0][0]][r[0][1]] = r[1]
                else:
                    row_data[r[0][0]] = {r[0][1] : r[1]}
            
            for _round, values in row_data.items():
                if _round in core_progression_dict.keys():
                    if len(values) == 2:
                        core_progression_dict[_round].append(values)
                    else:
                        core_progression_dict[_round].append({'Differences': values['difference'],
                                                              'Core': values['core'],
                                                              'Occurences': values['number_of_molecules']})
                else:
                    if len(values) == 2:
                        core_progression_dict[_round] = [values]
                    else:
                        core_progression_dict[_round] = [{'Differences': values['difference'],
                                                          'Core': values['core'],
                                                          'Occurences': values['number_of_molecules']}]
            
        return json.dumps({'Rounds': rounds,
                           'core_progression': core_progression_dict})
    else:
        return json.dumps({'Rounds': 'no rounds'})
    
@app.route('/<dataset>/<round_investigating>/chemical_map/')
def chemical_map(dataset, round_investigating):
    ''' Function: chemical_map
    This function return the chemical map data
    Input: dataset - name of dataset
           round_investigating - name of the rounds
    Output: json of the chemical map data
    '''
    if dataset == 'new_dataset':
        chemical_map_data = pd.read_csv(rnd.coordinates_file, sep='\t', header=0)
        core_data = pd.read_csv(rnd.node_information_file, sep='\t', header=0)
    else:
        chemical_map_data = pd.read_csv('datasets/{0}_coordinates.txt'.format(dataset), sep='\t', header=0)
        core_data = pd.read_csv('datasets/{0}_node_information.txt'.format(dataset), sep='\t', header=0)
    
    if round_investigating != 'norounds':
        core_data = core_data[core_data['Round'] == round_investigating]
        chemical_map_data = chemical_map_data[chemical_map_data['Round'] == round_investigating]
    
    core_data = core_data[['ID', 'core']]
    results = pd.merge(chemical_map_data, core_data, how='inner', on='ID')
    
    return results.to_json(orient='records')

@app.route('/molecule_core_finder/<dataset>/<molid>/', methods=['GET'])
def molecule_core_finder(dataset, molid):
    ''' Function: molecule_core_finder
    This function finds the cores that identify with a molecule (ID)
    Input: molidstr of the molecule ID
    Output: json of the core data
    '''
    if dataset == 'new_dataset':
        core_data = pd.read_csv(rnd.node_information_file, sep='\t', header=0)
    else:
        core_data = pd.read_csv('datasets/{0}_node_information.txt'.format(dataset), sep='\t', header=0)
    
    core_data_filtered = core_data[core_data['ID'] == molid]
    cores = core_data_filtered['core'].tolist()
    
    return json.dumps({molid: cores})

def generate_excel(file_name,
                   investigated_cores):
    ''' Function: generate_excel
    This function generates an excel file of the cores that are of interest
    Input: filename - str of excel file name to generate
           investigated_cores - list of string of the cores that are of interest
    Output: None
    '''
    overall_df = pd.DataFrame()
    
    core_table = ci.table_dictionary
    
    for core, table_info in core_table.items():
        if core in investigated_cores:
            core_df = pd.DataFrame.from_dict(table_info[0])
            core_df.columns = [[core]*len(core_df.columns), core_df.columns]
            overall_df = pd.concat([overall_df, core_df], axis=1)
    
    writer = pd.ExcelWriter(file_name, engine='xlsxwriter')
    overall_df.to_excel(writer, sheet_name='RG Core Visualisation')
    ### get an extra empty header row because of the multiindex
    writer.sheets['RG Core Visualisation'].set_row(2, None, None, {'hidden': True})
    workbook = writer.book
    worksheet = writer.sheets['RG Core Visualisation']
    
    tmpfnames = []
    for row_idx, row in enumerate(overall_df.iterrows()):
        for node_idx, node in enumerate(row[1]):
            if type(node) == str:
                ### image creation
                filename = 'img{0}{1}.jpg'.format(row_idx, node_idx)
                tmpfnames.append(filename)
                with open(filename, 'wb') as f:
                    ### remove text from cell
                    worksheet.write(row_idx+3, node_idx+1, ' ')
                    ### change bit64 to image 
                    f.write(base64.b64decode(node))
                    worksheet.insert_image(row_idx+3, node_idx+1, filename, {'positioning':1, 'x_offset':15, 'y_offset':5, 'x_scale':0.75, 'y_scale':0.75})

    worksheet.set_default_row(90)
    worksheet.set_column(0, len(overall_df.columns)+1, 20)
    workbook.close()
    writer.save()
    
    for fname in tmpfnames:
        if os.path.isfile(fname):
            os.remove(fname)
    
    return None

def generate_mol_excel(file_name,
                       mol_table):
    ''' Function: generate_excel
    This function generates an excel file of all the molecules within a RG core
    Input: filename - str of the file name
           mol_table - pandas df of the df to save
    Output: None
    '''
    writer = pd.ExcelWriter(file_name, engine='xlsxwriter')
    mol_table.to_excel(writer, sheet_name='RG Core Visualisation', index=False)
    workbook = writer.book
    worksheet = writer.sheets['RG Core Visualisation']
    
    tmpfnames = []
    for row_idx, row in mol_table.iterrows():
        node_idx= 1
        filename = 'img{0}{1}.jpg'.format(row_idx, node_idx)
        tmpfnames.append(filename)
        
        with open(filename, 'wb') as f:
            ### remove text from cell
            worksheet.write(row_idx+1, node_idx, ' ')
            ### change bit64 to image 
            f.write(base64.b64decode(row['Image']))
            worksheet.insert_image(row_idx+1, node_idx, filename, {'positioning':1, 'x_offset':15, 'y_offset':5, 'x_scale':0.75, 'y_scale':0.75})
    
    worksheet.set_default_row(300)
    worksheet.set_column("B:B", 60)
    workbook.close()
    writer.save()
    
    for fname in tmpfnames:
        os.remove(fname)
    
    return None

@app.route('/<dataset>/core/<core>/save_table/', methods=['GET','POST'])
def molecule_save_table(dataset, core):
    ''' Function: molecule_save_table
    This function downloads the table of the molecular information for core <core>
    Input: dataset - str of the name of the dataset of interest
           core - str of the core of interest
    Output: xlsx file
    '''
    file_name = ''.join([dataset, '_', core, '_table.xlsx'])
    
    ##need to make the table excel
    mol_table_info = ci.core_dictionary[core]['molecular_information']
    mol_table = pd.DataFrame.from_dict(mol_table_info)
    
    generate_mol_excel(file_name,
                       mol_table)
    
    return send_file(file_name,
                     mimetype='text/xlsx',
                     attachment_filename=file_name,
                     as_attachment=True)
    
@app.route('/<dataset>/comparator/save_table/', methods=['GET','POST'])
def comparator_save_table(dataset):
    ''' Function comparator_save_table
    This function downloads the table from the core comparator page
    Input: dataset - str of name of dataset
    Output: xlsx file
    '''    
    ### need to make the table excel
    ### need to filter the cores comparator
    cores = request.get_json().split('\t')
    
    cores_file_names = '_'.join(cores)
    file_name = ''.join([dataset, '_', cores_file_names, '_node_comparator_table.xlsx'])
     
    generate_excel(file_name,
                   cores)
    
    return send_file(file_name,
                     mimetype='text/xlsx',
                     attachment_filename=file_name,
                     as_attachment=True)

@app.route('/<dataset>/all_cores/save_table/', methods=['GET','POST'])
def all_save_table(dataset):
    ''' Function all_save_table
    This function downloads the table of all the RG core information
    Input: dataset - str of the name of the dataset of interest
    Output: xlsx file
    '''
    file_name = ''.join([dataset, '_node_table.xlsx'])
    cores = ci.table_dictionary.keys()
    
    generate_excel(file_name,
                   cores)
    
    return send_file(file_name,
                     mimetype='application/vnd.openxmlformatsofficedocument.spreadsheetml.sheet',
                     attachment_filename=file_name,
                     as_attachment=True)
    
def allowed_file(filename):
    '''
    This function establises whether a upload file is in the correct format
    '''
    allowed_extensions = set(['txt', 'smi'])
    return '.' in filename and filename.rsplit('.', 1)[1] in allowed_extensions
    
@app.route('/running_new_dataset/', methods=['POST'])
def running_new_dataset():
    ''' Function running_new_dataset
    This function uploads a file that the user would like to use in the RG Core visualisation
    '''
    if request.method == 'POST':
        f=request.files['file']
        if f and allowed_file(f.filename):
            f.save(secure_filename(f.filename))
            rnd.initial_file = f.filename
    return 'None'

@app.route('/running_new_dataset/create_rgs/', methods=['GET', 'POST'])
def create_rgs():
    ''' Function create_rgs
    This function sets the parameters and creates RGs for the uploaded file
    '''
    ### need to set parameters
    rnd.minsim = request.get_json()['minsim']
    rnd.minnodes = request.get_json()['minnodes']
    rnd.fingerprintbitlength = request.get_json()['fp_len']
    rnd.rounddata = request.get_json()['round_data']
    rnd.roundnames = [name for name in request.get_json()['round_data_names'].split(' ')]
    
    if rnd.rounddata == False:
        os.system('python python/reduced_graph_code/reduced_graph.py -i {0} -o datasets/output_reduced_graphs.txt -b datasets/output_reduced_graphs.sdf'.format(rnd.initial_file, ))
    else:
        whole_file = pd.read_csv(rnd.initial_file, sep='\t', header=0)
        filtered_df_to_append = pd.DataFrame()
        for roundn in rnd.roundnames:
            filtered_df = whole_file[whole_file['Round'] == roundn]
            filtered_df_to_append = filtered_df_to_append.append(filtered_df)
            filtered_df_to_append.to_csv(rnd.initial_file.replace('.txt', '_{0}.txt'.format(roundn)), sep='\t', index=False)
            
            os.system('python python/reduced_graph_code/reduced_graph.py -i {0} -o {1} -b {2}'.format(rnd.initial_file.replace('.txt', '_{0}.txt'.format(roundn)), 'datasets/output_reduced_graphs_{0}.txt'.format(roundn), 'datasets/output_reduced_graphs_{0}.sdf'.format(roundn)))

    return 'None'

@app.route('/running_new_dataset/create_mcs/', methods=['GET'])
def create_mcs():
    ''' Function create_mcs
    This function creates RG MCS for the RGs within uploaded file
    '''
    if rnd.rounddata == False:
        os.system('python python/MCS/mcs_similarity_matrix.py -i datasets/output_reduced_graphs.txt -r datasets/output_reduced_graphs_mcs.txt')
    else:
        for roundn in rnd.roundnames:
            os.system('python python/MCS/mcs_similarity_matrix.py -i {0} -r {1}'.format('datasets/output_reduced_graphs_{0}.txt'.format(roundn), 'datasets/output_reduced_graphs_mcs_{0}.txt'.format(roundn)))
            
    return 'None'

@app.route('/running_new_dataset/create_rg_cores/', methods=['GET'])
def create_rg_cores():
    ''' Function create_rg_cores
    This function creates RG Cores for the RGs within uploaded file
    '''
    if rnd.rounddata == False:
        os.system('python python/reduced_graph_core_extraction/finding_cores_from_whole_dataset.py -r datasets/output_reduced_graphs.txt -m datasets/output_reduced_graphs_mcs.txt -o datasets/output_rg_core_extraction.txt -i {0} -d {1}'.format(rnd.minsim, rnd.minnodes))
    else:
        for roundn in rnd.roundnames:
            os.system('python python/reduced_graph_core_extraction/finding_cores_from_whole_dataset.py -r {0} -m {1} -o {2} -i {3} -d {4}'.format('datasets/output_reduced_graphs_{0}.txt'.format(roundn), 'datasets/output_reduced_graphs_mcs_{0}.txt'.format(roundn), 'datasets/output_rg_core_extraction_{0}.txt'.format(roundn), rnd.minsim, rnd.minnodes))

    return 'None'

@app.route('/running_new_dataset/create_vis_files/', methods=['GET'])
def create_vis_files():
    ''' Function create_vis_files
    This function creates visualisation files for the RGs within uploaded file
    '''
    if rnd.rounddata == False:
        os.system('python python/generating_files_for_visualisation/generating_file_for_visualisation_coordinates.py -r datasets/output_reduced_graphs.txt -s datasets/output_reduced_graphs.sdf -o datasets/output_coordinates.txt -n {0}'.format(rnd.fingerprintbitlength))
        os.system('python python/generating_files_for_visualisation/creating_visualisation_file.py -r datasets/output_reduced_graphs.txt -s datasets/output_reduced_graphs.sdf -a {0} -c datasets/output_rg_core_extraction.txt -o datasets/output_node_information.txt'.format(rnd.initial_file))
        os.system('python python/generating_files_for_visualisation/creating_core_breakdown_analysis_file.py -i datasets/output_node_information.txt -o datasets/output_core_analysis.txt')
        
        rnd.coordinates_file = 'datasets/output_coordinates.txt'
        rnd.node_information_file = 'datasets/output_node_information.txt'
        rnd.core_analysis_file = 'datasets/output_core_analysis_file.txt' 
    else:
        coords_file = pd.DataFrame()
        node_info_file = pd.DataFrame()
        core_analysis_file = pd.DataFrame()
        
        for roundn in rnd.roundnames:
            os.system('python python/generating_files_for_visualisation/generating_file_for_visualisation_coordinates.py -r {0} -s {1} -o {2} -n {3}'.format('datasets/output_reduced_graphs_{0}.txt'.format(roundn), 'datasets/output_reduced_graphs_{0}.sdf'.format(roundn), 'datasets/output_coordinates_{0}.txt'.format(roundn), rnd.fingerprintbitlength))
            os.system('python python/generating_files_for_visualisation/creating_visualisation_file.py -r {0} -s {1} -a {2} -c {3} -o {4}'.format('datasets/output_reduced_graphs_{0}.txt'.format(roundn), 'datasets/output_reduced_graphs_{0}.sdf'.format(roundn), rnd.initial_file, 'datasets/output_rg_core_extraction_{0}.txt'.format(roundn), 'datasets/output_node_information_{0}.txt'.format(roundn)))
            os.system('python python/generating_files_for_visualisation/creating_core_breakdown_analysis_file.py -i {0} -o {1}'.format('datasets/output_node_information_{0}.txt'.format(roundn), 'datasets/output_core_analysis_{0}.txt'.format(roundn)))
        
            coord_df = pd.read_csv('output_coordinates_{0}.txt'.format(roundn), sep='\t', header=0)
            coord_df['Round'] = roundn
            node_df = pd.read_csv('output_node_information_{0}.txt'.format(roundn), sep='\t', header=0)
            node_df['Round'] = roundn
            core_analysis_df = pd.read_csv('output_core_analysis_file_{0}.txt'.format(roundn), sep='\t', header=0)
            core_analysis_df['Round'] = roundn
            
            coords_file = coords_file.append(coord_df)
            node_info_file = node_info_file.append(node_df)
            core_analysis_file = core_analysis_file.append(core_analysis_df)
        
        coords_file.to_csv('datasets/output_coordinates.txt', sep='\t', index=False)
        node_info_file.to_csv('datasets/output_node_information.txt', sep='\t', index=False)
        core_analysis_file.to_csv('datasets/output_core_analysis_file.txt', sep='\t', index=False)
        
        rnd.coordinates_file = 'datasets/output_coordinates.txt'
        rnd.node_information_file = 'datasets/output_node_information.txt'
        rnd.core_analysis_file = 'datasets/output_core_analysis_file.txt'        

    return 'None'

if __name__ == "__main__":
    sv = static_variables()
    ci = core_information()
    rnd = running_new_dataset_class()
    dataset_json = json.load(open('datasets/datasets.json'))
    
    manager.run()
