"""A python program demonstrating simple SCF procedure.
"""

import numpy as np

# Loading one- and two- electron integrals
S = np.load('../data/h2o_oeints.npz')['overlap']
T = np.load('../data/h2o_oeints.npz')['kinetic']
V = np.load('../data/h2o_oeints.npz')['potential']
I = np.load('../data/h2o_erints.npz')['erints']
H_core = T + V

# ==> Inspecting S for AO orthonormality <==
hope = np.allclose(S, np.eye(S.shape[0]))
print('\nDo we have any hope that our AO basis is orthonormal? %s!' % (hope))

# Diagonalise the Overlap matrix
svals, svecs = np.linalg.eigh(S)
# Taking the inverse square root
X = np.linalg.inv(np.diag(np.sqrt(svals)))
# Using the eigen vector back transform from eigen basis to previous basis.
X = np.einsum('aI,IJ,Jb', svecs, X, svecs.T)

# Orthonormalising S
S_p = np.einsum('ij,jk,kl->il', X.T, S, X)
# Checking if orthonormalised
orthonormalised = np.allclose(S_p, np.eye(S.shape[0]), atol=1e-08)
if orthonormalised:
    print("There is a chance for diagonalisation")
else:
    print("There's something wrong with the transformation.")

# Transforming the Fock matrix
F_p = np.einsum('ij,jk,kl->il', X.T, H_core, X)
# Diagonalising the Fock matrix eigenvalues and eigenvectors
e, C_p = np.linalg.eigh(F_p)
# Transform C_p back into AO basis
C = np.dot(X, C_p)
# slicing out occupied orbitals
ndocc = 5
C_occ = C[:, :ndocc]
# Building density matrix from C_occ
D = np.einsum('pi,qi->pq', C_occ, C_occ, optimize=True)


# ==> SCF Iterations <==
# Maximum SCF iterations
MAXITER = 50
# Energy convergence criterion
E_conv = 1.0e-10
# Pre-iteration energy declarations
SCF_E = 0.0
E_old = 0.0
print('==> Starting SCF Iterations <==\n')
# Begin Iterations
for scf_iter in range(1, MAXITER + 1):
    # Build Fock matrix
    J = np.einsum('pqrs,rs->pq', I, D, optimize=True)
    K = np.einsum('prqs,rs->pq', I, D, optimize=True)
    F = H_core + 2*J - K
    # Compute RHF energy
    SCF_E = np.einsum('pq,pq->', (H_core + F), D, optimize=True)
    print('SCF Iteration %3d: Energy = %4.16f dE = % 1.5E' % (scf_iter, SCF_E, SCF_E - E_old))
    
    # SCF Converged?
    if (abs(SCF_E - E_old) < E_conv):
        break
    E_old = SCF_E
    
    # Compute new orbital guess
    F_p =  np.einsum('ij,jk,kl->il', X, F, X)
    e, C_p = np.linalg.eigh(F_p)
    C = np.dot(X, C_p)
    C_occ = C[:, :ndocc]
    D = np.einsum('pi,qi->pq', C_occ, C_occ, optimize=True)
    
    # MAXITER exceeded?
    if (scf_iter == MAXITER):
        raise Exception("Maximum number of SCF iterations exceeded.")
        
# Post iterations
print('\nSCF converged.')
print('Final RHF Energy: %.8f [Eh]' % (SCF_E))