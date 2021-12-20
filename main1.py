import sys
from mol2pdb import *
import os
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
import pandas as pd
import glob

if sys.argv[1]==3:
    seq = interprete(sys.argv[2])
else:
    seq = sys.argv[2]
    
mol = rdmolfiles.MolFromPDBFile(sys.argv[4])#converting the amino acid sequence to mol format

dihedral = get_dihedrals(mol)

AllChem.EmbedMolecule(mol, randomSeed=1) #converting the 2d molecule to 3d
if os.path.isdir('{}'.format(seq))==False:
    os.mkdir(seq)

rdmolfiles.MolToPDBFile(mol, '{}/{}.pdb'.format(seq, seq)) #writing mol to a pdb file

conf=mol.GetConformer(0) #get a conformation from the mol file


# for m,a in enumerate(dihedral):
#     for i in range(int(sys.argv[3])):
#         angle = Chem.rdMolTransforms.GetDihedralDeg(conf, a[0], a[1], a[2], a[3])
#         print('the {} dihedral angle between bonds {}'.format(a, angle))
#         newangle = np.random.uniform(-180, 180)
#         Chem.rdMolTransforms.SetDihedralDeg(conf, a[0], a[1], a[2], a[3], newangle)
#         print('New Angle: {}'.format(Chem.rdMolTransforms.GetDihedralDeg(conf, a[0], a[1], a[2], a[3])))
#         file = '{}/{}_{}_{}.pdb'.format(seq, 'dih', m, i)
#         Chem.MolToPDBFile(mol, file)
#         pdb = PDBFile(file)
#         forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
#         modeller = Modeller(pdb.topology, pdb.positions)
#         with open((file), 'w') as f:
#             PDBFile.writeFile(modeller.topology, modeller.positions, f)
#         file_path = clean_pdb(file, '{}/{}_{}_{}_Hs.pdb'.format(seq, 'dih', m, i))
#         if os.path.isdir('{}/minimized'.format(seq))==False:
#             os.mkdir('{}/minimized'.format(seq))
#         ffminimization(file_path, '{}/minimized/{}_{}_{}_Hs_output.pdb'.format(seq, seq, m, i))
# print('done')
df1 = pd.read_csv('dunbrack_peptide_lib.csv')
a = np.arange(-180, 180, 10)
for i in range(0, len(dihedral), 3):
    for j in range(int(sys.argv[3])):
        if 'P' in seq:
            idx = seq.index('P')
            if i==(idx-1)*3:
                psi_ang = float(np.random.choice(a))
                psi_atm = dihedral[i]
                Chem.rdMolTransforms.SetDihedralDeg(conf, psi_atm[0], psi_atm[1], psi_atm[2], psi_atm[3], psi_ang)
        else:
            phi_atm = dihedral[i+2]
            psi_atm = dihedral[i]
            omega_atm = dihedral[i+1]
            phi_ang = float(np.random.choice(a))
            psi_ang = float(np.random.choice(a))
            omega = df1.loc[(df1['Phi(+1)']==-phi_ang) & (df1['Psi(0)']==psi_ang) & (df1['ResTypeGroup']=='All')][['mW(+1)', 'sW(+1)']]
            omega_ang = np.random.uniform(float(omega['mW(+1)']) - float(omega['sW(+1)']),float(omega['mW(+1)']) + float(omega['sW(+1)']))
            print(Chem.rdMolTransforms.GetDihedralDeg(conf, phi_atm[0], phi_atm[1], phi_atm[2], phi_atm[3]), phi_ang)
            Chem.rdMolTransforms.SetDihedralDeg(conf, phi_atm[0], phi_atm[1], phi_atm[2], phi_atm[3], phi_ang)
            print(Chem.rdMolTransforms.GetDihedralDeg(conf, psi_atm[0], psi_atm[1], psi_atm[2], psi_atm[3]), psi_ang)
            Chem.rdMolTransforms.SetDihedralDeg(conf, psi_atm[0], psi_atm[1], psi_atm[2], psi_atm[3], psi_ang)
            print(Chem.rdMolTransforms.GetDihedralDeg(conf, omega_atm[0], omega_atm[1], omega_atm[2], omega_atm[3]), omega_ang)
            Chem.rdMolTransforms.SetDihedralDeg(conf, omega_atm[0], omega_atm[1], omega_atm[2], omega_atm[3], omega_ang)
            
        
        file = '{}/{}_{}_{}.pdb'.format(seq,seq, i, j)
        Chem.MolToPDBFile(mol, file)
        pdb = PDBFile(file)
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller = Modeller(pdb.topology, pdb.positions)
        with open((file), 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)
        file_path = clean_pdb(file, '{}/{}_{}_{}_Hs.pdb'.format(seq,seq, i, j))
        if os.path.isdir('{}/minimized'.format(seq))==False:
            os.mkdir('{}/minimized'.format(seq))
        ffminimization(file_path, '{}/minimized/{}_{}_{}_Hs_output.pdb'.format(seq,seq, i, j))
print('done')

files = glob.glob('{}/minimized/*.pdb'.format(seq))
for i,file in enumerate(files):
    point_mutation(sys.argv[4], file, i)
