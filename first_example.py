#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# My first example with AutoDock Vina in python
#

from vina import Vina
import sys

print(5)
v = Vina(sf_name='vina')

v.set_receptor('mod_prot.pdbqt')

v.set_ligand_from_file('mod_pep.pdbqt')
v.compute_vina_maps(center=[-2, 19, 3], box_size=[25, 25, 25])

# Score the current pose
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

# Minimized locally the current pose
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
v.write_pose('prot_ligand_minimized.pdbqt', overwrite=True)

# Dock the ligand
v.dock(exhaustiveness=32, n_poses=20)
v.write_poses('prot_ligand_vina_out.pdbqt', n_poses=5, overwrite=True)
