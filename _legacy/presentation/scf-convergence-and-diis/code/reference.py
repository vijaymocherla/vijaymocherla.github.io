import psi4
import sys
import numpy as np 
import pandas as pd 
# setting memory, output, num_threads for psi4
psi4.set_memory('256 Gb')
numpy_memory = 2
psi4.core.set_output_file("e_reference.dat", False)
psi4.core.set_num_threads(11)
# setting other options : basis, ints_tolerance etc.
options = {'ints_tolerance': 1e-12,
            'scf_type': 'pk',
            'e_convergence': 1e-12,
            'd_convergence': 1e-8}
psi4.set_options(options)
# Loading the molecule object
mol = psi4.geometry("""
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
""")
# energy benchmarks
energies = {'nbf':[],
            'scf':[], 'mp2': [],  
            'ccsd':[], 'ccsd(t)':[],
            'pbe0':[], 'b3lyp':[]
            }
basissets = ['sto-3g', 'sto-6g', '3-21g', 
         '6-31g', '6-31g*', '6-31+g*',
         'cc-pvdz', 'cc-pvtz', 'cc-pvqz',
         'aug-cc-pvdz', 'aug-cc-pvtz', 'aug-cc-pvqz',
         'def2-sv(p)', 'def2-svp']

#geom_opt = psi4.optimize('mp2/def2-sv(p)')
#mol.update_geometry() 

for basis in basissets:    
    wfn = psi4.core.Wavefunction.build(mol, basis)
    energies['nbf'].append(wfn.basisset().nbf())
    e_scf, wfn_scf = psi4.energy('scf/'+basis, return_wfn=True)
    e_mp2, wfn_mp2 = psi4.energy('mp2/'+basis, return_wfn=True)
    e_ccsd, wfn_ccsd = psi4.energy('ccsd/'+basis, return_wfn=True)
    e_ccsd_t, wfn_ccsd_t = psi4.energy('ccsd(t)/'+basis, return_wfn=True)
    e_pbe0, wfn_pbe0 = psi4.energy('pbe0/'+basis, return_wfn=True)
    e_b3lyp, wfn_b3lyp = psi4.energy('b3lyp/'+basis, return_wfn=True)
    energies['scf'].append(e_scf - mol.nuclear_repulsion_energy())
    energies['mp2'].append(e_mp2 - mol.nuclear_repulsion_energy())
    energies['ccsd'].append(e_ccsd - mol.nuclear_repulsion_energy())
    energies['ccsd(t)'].append(e_ccsd_t - mol.nuclear_repulsion_energy())
    energies['pbe0'].append(e_pbe0 - mol.nuclear_repulsion_energy())
    energies['b3lyp'].append(e_b3lyp - mol.nuclear_repulsion_energy())  

# e_scf, e_cisd, e_ccsd, e_ccsd_t, e_pbe0, e_b3lyp
df = pd.DataFrame(energies, index=basissets)
print(df.describe())
print(df)
df.to_csv('energies_benchmark.csv')

