import random
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import glob


def interprete(seq):
    new_seq = ''
    seq = seq.split('-')
    interprete_dict = {'Arg': 'R', 'His': 'H', 'Lys': 'K', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Thr': 'T', 'Asn': 'N',
                       'Gln': 'Q', 'Cys': 'C', 'Sec': 'U', 'Gly': 'G', 'Pro': 'P', 'Ala': 'A', 'Ile': 'I', 'Leu': 'L',
                       'Met': 'M', 'Phe': 'F', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Dap': '1', 'Dab': '2',
                       'BOrn': '3', 'BLys': '4', 'Hyp': 'Z', 'Orn': 'O', 'bAla': '!', 'Gaba': '?', 'dDap': '5',
                       'dDab': '6',
                       'dBOrn': '7', 'dBLys': '8', 'dArg': 'r', 'dHis': 'h', 'dLys': 'k', 'dAsp': 'd', 'dGlu': 'e',
                       'dSer': 's',
                       'dThr': 't', 'dAsn': 'n', 'dGln': 'q', 'dCys': 'c', 'dSec': 'u', 'dGly': 'g', 'dPro': 'p',
                       'dAla': 'a',
                       'dIle': 'i', 'dLeu': 'l', 'dMet': 'm', 'dPhe': 'f', 'dTrp': 'w', 'dTyr': 'y', 'dVal': 'v',
                       'dHyp': 'z', 'dOrn': 'o', 'a5a': '=', 'a6a': '%', 'a7a': '$', 'a8a': '@', 'a9a': '#',
                       'Cys1': 'Ä', 'Cys2': 'Ö', 'Cys3': 'Ü', 'dCys1': 'ä', 'dCys2': 'ö', 'dCys3': 'ü',
                       'Ac': '&', 'NH2': '+', 'met': '-', 'cy': 'X'}
    for bb in seq:
        new_seq += interprete_dict[bb]
    seq = new_seq
    return seq




def clean_pdb(file, path):
    pdb = PDBFile(file)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, model='tip3p', padding=1*nanometer)
    modeller.deleteWater()
    with open(path, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    return path
        
    


        
        
def get_dihedrals(mol):
    dihedral = []
    for i, bond in enumerate(mol.GetBonds()):
        if bond.GetBondTypeAsDouble()==1 and bond.IsInRing()==False:
            atomB = bond.GetBeginAtom()
            atomE = bond.GetEndAtom()
        
            neighB = []
            neighE = []


            for atom in atomB.GetNeighbors():
                if atom.GetIdx() != atomE.GetIdx():
                    neighB.append(atom.GetIdx())
                
            for atom in atomE.GetNeighbors():
                if atom.GetIdx() != atomB.GetIdx():
                    neighE.append(atom.GetIdx())

            if neighB==[] or neighE==[]:
                continue
        
            if len(neighB)>len(neighE):
            
                for i in range(len(neighB)):
                    for j in range(len(neighE)):
                        dihedral.append([neighB[i], atomB.GetIdx(), atomE.GetIdx(), neighE[j]])
                    
            else:
                for i in range(len(neighE)):
                    for j in range(len(neighB)):
                        dihedral.append([neighB[j], atomB.GetIdx(), atomE.GetIdx(), neighE[i]])
    return dihedral        

        

    
    
    
def ffminimization(file, output_path):
    pdb = PDBFile(file)
    forcefield = ForceField('amber14-all.xml', 'amber14/protein.ff14SB.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = VerletIntegrator(0.001*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(maxIterations=100)
    print('Saving...')
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(output_path, 'w'))
    
    
   
def backbone_dihedrals(mol):
    N=[0]
    substructure = Chem.MolFromSmarts('CNC(=O)')
    for item in mol.GetSubstructMatches(substructure):
        N.append(item[1])
    N = list(set(N))
    N.sort()
    dihedral = []
    for i in range(len(N)-1):
        dihedral.append([N[i], N[i]+1, N[i]+2, N[i+1]])
        dihedral.append([N[i]+1, N[i]+2, N[i+1], N[i+1]+1])
        dihedral.append([N[i]+2, N[i+1], N[i+1]+1, N[i+1]+2])
    
    return dihedral
