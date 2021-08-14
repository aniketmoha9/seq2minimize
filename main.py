import sys
from mol2pdb import *
import os
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem

if sys.argv[1]==3:
    seq = interprete(sys.argv[2])
else:
    seq = sys.argv[2]
    
mol = Chem.MolFromSequence(seq) #converting the amino acid sequence to mol format

dihedral = get_dihedrals(mol)

AllChem.EmbedMolecule(mol, randomSeed=1) #converting the 2d molecule to 3d
if os.path.isdir('{}'.format(seq))==False:
    os.mkdir(seq)

rdmolfiles.MolToPDBFile(mol, '{}/{}.pdb'.format(seq, seq)) #writing mol to a pdb file

conf=mol.GetConformer(0) #get a conformation from the mol file


for m,a in enumerate(dihedral):
    for i in range(int(sys.argv[3])):
        angle = Chem.rdMolTransforms.GetDihedralDeg(conf, a[0], a[1], a[2], a[3])
        print('the {} dihedral angle between bonds {}'.format(a, angle))
        newangle = np.random.uniform(-180, 180)
        Chem.rdMolTransforms.SetDihedralDeg(conf, a[0], a[1], a[2], a[3], newangle)
        print('New Angle: {}'.format(Chem.rdMolTransforms.GetDihedralDeg(conf, a[0], a[1], a[2], a[3])))
        file = '{}/{}_{}_{}.pdb'.format(seq, 'dih', m, i)
        Chem.MolToPDBFile(mol, file)
        pdb = PDBFile(file)
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller = Modeller(pdb.topology, pdb.positions)
        with open((file), 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)
        file_path = clean_pdb(file, '{}/{}_{}_{}_Hs.pdb'.format(seq, 'dih', m, i))
        if os.path.isdir('{}/minimized'.format(seq))==False:
            os.mkdir('{}/minimized'.format(seq))
        ffminimization(file_path, 'seq/minimized/{}_{}_{}_Hs_output.pdb'.format(seq, m, i))
print('done')

