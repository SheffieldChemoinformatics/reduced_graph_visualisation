##############################################################################
# Reduced graph program to generate file for visualisation plots
#
#
# Jess Stacey 
##############################################################################

import argparse
import pandas as pd
import numpy as np
import zlib
import base64
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def input_args():
    parser = argparse.ArgumentParser(description='Arguments for the Reduced Graph Code')
    parser.add_argument('-r',
                        dest='reduced_graphs',
                        default=r'..\example_input_files\output_reduced_graphs_file.txt',
                        help='Name of the input rg file')
    parser.add_argument('-s',
                        dest='sdf',
                        default=r'..\example_input_files\output_reduced_graphs_file.sdf',
                        help='Name of the input sdf file')
    parser.add_argument('-o',
                        dest='output',
                        default=r'..\example_input_files\output_coordinates.txt',
                        type=str,
                        help='Name of the file containing all the information and values')
    parser.add_argument('-n',
                        dest='number_of_bits',
                        default=2048,
                        type=int,
                        help='Number of bit of the fingerprint')
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

def find_ecfp(smiles_mol):
    ecfp = AllChem.GetMorganFingerprintAsBitVect(smiles_mol, 2, nBits=args.number_of_bits)
    
    return ecfp

def generate_rg():
    reduced_graphs = pd.read_csv(args.reduced_graphs, header=0, sep=args.separator)
    
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
    fingerprint_list = []
    for r, v in reduced_graphs.iterrows():
        ecfp = find_ecfp(Chem.MolFromSmiles(v['SMILES']))
        fingerprint_list.append(ecfp)
        if v['SMILES'] in molecules_dict.keys():
            reduced_graphs.at[r, 'SDF'] = molecules_dict[v['SMILES']]
      
    return reduced_graphs, molecules_dict, fingerprint_list

def generate_coords(merge_data,
                    fingerprint_list):    
    fingerprint_list = [list(fingerprint_key.ToBitString()) for fingerprint_key in fingerprint_list]
    fingerprint_df = pd.DataFrame(fingerprint_list)
    
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(fingerprint_df)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['pc1', 'pc2'])
    coords = pd.concat([merge_data, principalDf], axis = 1)
    
    X = np.array(fingerprint_list)
    X_embedded = TSNE(n_components=2).fit_transform(X)
    principalDf = pd.DataFrame(data = X_embedded, columns = ['tsne1', 'tsne2'])
    
    coords = pd.concat([coords, principalDf], axis = 1)
    
    return coords
    
if __name__ == '__main__':
    args = input_args()
    
    reduced_graphs, sdf_file, fingerprint_list = generate_rg()
    rg_plus_coords_df = generate_coords(reduced_graphs, fingerprint_list)
    if 'number_of_heavy_atoms' in rg_plus_coords_df.columns:
        del rg_plus_coords_df['number_of_heavy_atoms']
    rg_plus_coords_df.columns = ['SMILES', 'ID', 'RG', 'SDF', 'x_pca', 'y_pca', 'x_tsne', 'y_tsne']

    rg_plus_coords_df.to_csv(args.output, sep='\t', index=False)