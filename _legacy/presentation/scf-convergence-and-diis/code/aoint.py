"""Psi4 code in python to generate AO integrals for H20.
"""
__author__ = "Sai Vijay Mocherla"

import psi4
import numpy as np 
import sys
# setting memory, output, num_threads for psi4
psi4.set_memory('1 Gb')
numpy_memory = 2
psi4.core.set_output_file("../data/aoint.dat", False)
psi4.core.set_num_threads(2)
# setting other options : basis, ints_tolerance etc.
basis = sys.argv[1]
options = {'ints_tolerance': 1e-12,
            'basis': basis,
            'scf_type': 'pk'}
psi4.set_options(options)
# Loading the molecule object
mol = psi4.geometry("""
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
""")
# geom_opt = psi4.optimize('mp2/def2-sv(p)')
# mol.update_geometry() 
# Let's create a wavefunction object and AOintegrals
psi4.set_options(options)
wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('basis'))
mints = psi4.core.MintsHelper(wfn.basisset())
# One-electron integrals
S = np.asarray(mints.ao_overlap())
T = np.asarray(mints.ao_kinetic())
V = np.asarray(mints.ao_potential())
np.savez('../data/h2o_oeints.npz', overlap=S, kinetic=T, potential=V)
del(S, T, V) # Deleting this objects from memory
# Two-electron integrals
# Memory check for ERI tensor
nbf = wfn.nmo()
I_size = (nbf**4) * 8.e-9
print('\nSize of the ERI tensor will be {:4.2f} GB.'.format(I_size))
if I_size > numpy_memory:
    psi4.core.clean()
    raise Exception("Estimated memory utilization (%4.2f GB) exceeds allotted memory \
                     limit of %4.2f GB." % (I_size, numpy_memory))
# Building electron repulsion integral tensor
I = np.asarray(mints.ao_eri())
np.savez("../data/h2o_erints.npz", erints=I)
del(I)
psi4.core.clean()