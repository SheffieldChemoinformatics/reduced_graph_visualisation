# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 10:30:09 2017

@author: Jessica Stacey
"""

from rdkit import Chem

def canonicalise_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        canonical_smiles = Chem.MolToSmiles(mol, True)
    except:
        print('Bad SMILES unable to read ', smiles)
        return(False, False)
    
    return(canonical_smiles)
    
def canonicalise_smarts(smarts, n):
    try:
        mol = Chem.MolFromSmarts(smarts)
        canonical_smarts = Chem.MolToSmarts(mol)
    except:
        if n ==0:
            raise ValueError(smarts+' SMARTS cannot be read')
        if n==1:
            raise ValueError(smarts+' Isolating SMARTS cannot be read')
    return (canonical_smarts, mol)

def canonicalise_reducedgraph(reduced_graph):
    try:
        params = Chem.SmilesParserParams()
        params.removeHs=True
        params.sanitize=False
        mol = Chem.MolFromSmiles(reduced_graph, params)
        canonical_reduced_graph = Chem.MolToSmiles(mol, isomericSmiles=True)
    except:
        print('Bad SMILES unable to read reduced graph', reduced_graph)
        return(False, False)
    
    return(canonical_reduced_graph)