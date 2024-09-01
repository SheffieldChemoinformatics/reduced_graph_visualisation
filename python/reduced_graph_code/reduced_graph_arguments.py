##############################################################################
# Reduced graph program
# Arguments
#
#
# Jess Stacey
##############################################################################

import argparse

def input_args():
    parser = argparse.ArgumentParser(description='Arguments for the Reduced Graph Code')
    parser.add_argument('-i',
                        dest='smiles',
                        default=r'..\example_input_files\p2x7_gsk_smiles.txt',
                        type=str,
                        help='Name of the smiles input file')
    parser.add_argument('-s',
                        dest='smarts',
                        default='smarts.smt',
                        type=str,
                        help='Name of the smarts file')
    parser.add_argument('-c',
                        dest='carbon',
                        default='isolating_carbon.smt',
                        type=str,
                        help='Name of the isolating smarts file')
    parser.add_argument('-o',
                        dest='output',
                        default=r'..\example_input_files\output_reduced_graphs_p2x7_gsk.txt',
                        type=str, help='Name of the output file')
    parser.add_argument('-b',
                        dest='sdf',
                        default=r'..\example_input_files\outpt_reduced_graphs_p2x7_gsk.sdf',
                        type=str, help='Name of the sdf output file')
    parser.add_argument('-r',
                        dest='ring',
                        default=60,
                        type=int,
                        help='Maximum number of atoms in a ring node')
    parser.add_argument('-m',
                        dest='metal',
                        action='store_true',
                        help='Defines metal atoms as singular metal nodes')
    parser.add_argument('-l',
                        dest='terminallinker',
                        action='store_true',
                        help='Defines terminal carbon nodes as acyclic inert nodes')
    parser.add_argument('-p',
                        dest='hydrophobic',
                        action='store_true',
                        help='Defines anything that has not been previously defined \
                        that could be of interest, i.e. anything that isn\'t a straight chain of carbons')
    parser.add_argument('-d',
                        dest='doublebond',
                        action='store_true',
                        help='Defines whether want double bonds between nodes \
                        other than just fused rings to be defined')
    parser.add_argument('-t',
                        dest='separator',
                        default='\t',
                        help='File separator in the input smiles file')
    parser.add_argument('-a',
                        dest='additionalpIC50',
                        action='store_true',
                        help='Whether you want the pIC50s in the sdf file')
    parser.add_argument('-f',
                        dest='pIC50inputfile',
                        default=r'..\example_input_files\p2x7_gsk_pIC50_data.txt',
                        type=str, help='Name of the pIC50 input file')
    parser.add_argument('-v',
                        dest='verbose',
                        action='store_true',
                        help='Switches on verbose mode')

    args = parser.parse_args()
    return args
