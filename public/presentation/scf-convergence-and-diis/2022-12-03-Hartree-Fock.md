---
layout: post
title:  "Accelerating SCF with DIIS"
date:   2023-01-12 05:30:00 -0500
categories: "quantum-chemistry" 
author_link: "/about.html"
author: Sai Vijay Mocherla
---

In this, notebook we learn : - How to apply the **SCF**(*self-consistent
field*) method to obtain a ground state energy for a closed-shell system
using Restricted Hartree-Fock formalism. - How to acclerate **SCF**
using an extrapolation method called **DIIS**(*Direct Inversion of the
Iterative Subspace*).

<center>
<img src="scf-convergence-and-diis/slides/0001.jpg" alt="slide 1" width="800" height="600" />
<center/>


<center>
<img src="scf-convergence-and-diis/slides/0002.jpg" alt="slide 2" width="800" height="600" />
<center/>
<center>
<img src="scf-convergence-and-diis/slides/0003.jpg" alt="slide 3" width="800" height="600" />
<center/>
<center>
<img src="scf-convergence-and-diis/slides/0004.jpg" alt="slide 4" width="800" height="600" />
<center/>

``` python
# Getting the Paraphernelia
import numpy as np
```

``` python
# Generating AO integrals and storing them as .npz files
!python ../code/aoint.py sto-3g
```


      Memory set to 953.674 MiB by Python driver.

    Size of the ERI tensor will be 0.00 GB.

``` python
S = np.load('../data/h2o_oeints.npz')['overlap']
T = np.load('../data/h2o_oeints.npz')['kinetic']
V = np.load('../data/h2o_oeints.npz')['potential']
I = np.load('../data/h2o_erints.npz')['erints']
H_core = T + V
```

<center>
<img src="scf-convergence-and-diis/slides/0005.jpg" alt="slide 5" width="800" height="600" />
<center/>

``` python
# ==> Inspecting S for AO orthonormality <==
hope = np.allclose(S, np.eye(S.shape[0]))
print('\nDo we have any hope that our AO basis is orthonormal? %s!' % (hope))
```


    Do we have any hope that our AO basis is orthonormal? False!

``` python
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
```

    There is a chance for diagonalisation

<center>
<img src="scf-convergence-and-diis/slides/0006.jpg" alt="slide 6" width="800" height="600" />
<center/>

``` python
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
```

``` python
# ==> SCF Iterations <==
# Maximum SCF iterations
MAXITER = 50
# Energy convergence criterion
E_conv = 1.0e-10
D_conv = 1.0e-8
# Pre-iteration energy declarations
SCF_E = 0.0
E_old = 0.0
```

``` python
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
    e, C_p = np.linalg.eig(F_p)
    C = np.dot(X, C_p)
    C_occ = C[:, :ndocc]
    D = np.einsum('pi,qi->pq', C_occ, C_occ, optimize=True)
    
    # MAXITER exceeded?
    if (scf_iter == MAXITER):
        raise Exception("Maximum number of SCF iterations exceeded.")
# Post iterations
print('\nSCF converged.')
print('Final RHF Energy: %.8f [Eh]' % (SCF_E))
```

    ==> Starting SCF Iterations <==

    SCF Iteration   1: Energy = -81.2881629233857979 dE = -8.12882E+01
    SCF Iteration   2: Energy = -80.7620614193841675 dE =  5.26102E-01
    SCF Iteration   3: Energy = -80.4876942697945026 dE =  2.74367E-01
    SCF Iteration   4: Energy = -80.4657917548717592 dE =  2.19025E-02
    SCF Iteration   5: Energy = -80.4616027349987775 dE =  4.18902E-03
    SCF Iteration   6: Energy = -80.4606338529288507 dE =  9.68882E-04
    SCF Iteration   7: Energy = -80.4604211214459184 dE =  2.12731E-04
    SCF Iteration   8: Energy = -80.4603728090642534 dE =  4.83124E-05
    SCF Iteration   9: Energy = -80.4603619905947767 dE =  1.08185E-05
    SCF Iteration  10: Energy = -80.4603595507151681 dE =  2.43988E-06
    SCF Iteration  11: Energy = -80.4603590022545205 dE =  5.48461E-07
    SCF Iteration  12: Energy = -80.4603588787711175 dE =  1.23483E-07
    SCF Iteration  13: Energy = -80.4603588509901186 dE =  2.77810E-08
    SCF Iteration  14: Energy = -80.4603588447377547 dE =  6.25236E-09
    SCF Iteration  15: Energy = -80.4603588433308659 dE =  1.40689E-09
    SCF Iteration  16: Energy = -80.4603588430142764 dE =  3.16589E-10
    SCF Iteration  17: Energy = -80.4603588429430090 dE =  7.12674E-11

    SCF converged.
    Final RHF Energy: -80.46035884 [Eh]

## Let’s try for a slightly larger basis

``` python
# Generating AO integrals for cc-pvtz
!python ../code/aoint.py cc-pvtz
```


      Memory set to 953.674 MiB by Python driver.

    Size of the ERI tensor will be 0.09 GB.

``` python
S = np.load('../data/h2o_oeints.npz')['overlap']
T = np.load('../data/h2o_oeints.npz')['kinetic']
V = np.load('../data/h2o_oeints.npz')['potential']
I = np.load('../data/h2o_erints.npz')['erints']
H_core = T + V
```

``` python
# Running scf code the same as above
!python ../code/scf.py
```


    Do we have any hope that our AO basis is orthonormal? False!
    There is a chance for diagonalisation
    ==> Starting SCF Iterations <==

    SCF Iteration   1: Energy = -69.1347968401195203 dE = -6.91348E+01
    SCF Iteration   2: Energy = -73.8555083348786212 dE = -4.72071E+00
    SCF Iteration   3: Energy = -78.9567919922574930 dE = -5.10128E+00
    SCF Iteration   4: Energy = -77.5705309079189647 dE =  1.38626E+00
    SCF Iteration   5: Energy = -80.4078732350268410 dE = -2.83734E+00
    SCF Iteration   6: Energy = -78.5636142545179439 dE =  1.84426E+00
    SCF Iteration   7: Energy = -80.8401891040232670 dE = -2.27657E+00
    SCF Iteration   8: Energy = -78.9968790717921650 dE =  1.84331E+00
    SCF Iteration   9: Energy = -81.0445277765533660 dE = -2.04765E+00
    SCF Iteration  10: Energy = -79.2502300566142992 dE =  1.79430E+00
    SCF Iteration  11: Energy = -81.1702581595921941 dE = -1.92003E+00
    SCF Iteration  12: Energy = -79.4207412202195684 dE =  1.74952E+00
    SCF Iteration  13: Energy = -81.2574417841521779 dE = -1.83670E+00
    SCF Iteration  14: Energy = -80.6382468146913993 dE =  6.19195E-01
    SCF Iteration  15: Energy = -81.8060592445228565 dE = -1.16781E+00
    SCF Iteration  16: Energy = -81.0375710019176410 dE =  7.68488E-01
    SCF Iteration  17: Energy = -82.0323551720707087 dE = -9.94784E-01
    SCF Iteration  18: Energy = -81.2104173916206662 dE =  8.21938E-01
    SCF Iteration  19: Energy = -82.1328210092766255 dE = -9.22404E-01
    SCF Iteration  20: Energy = -82.5674752866804624 dE = -4.34654E-01
    SCF Iteration  21: Energy = -83.1491071642107187 dE = -5.81632E-01
    SCF Iteration  22: Energy = -83.4617045573220651 dE = -3.12597E-01
    SCF Iteration  23: Energy = -83.7115925113282628 dE = -2.49888E-01
    SCF Iteration  24: Energy = -83.8462016480608270 dE = -1.34609E-01
    SCF Iteration  25: Energy = -83.9288413708520835 dE = -8.26397E-02
    SCF Iteration  26: Energy = -83.9716034037429750 dE = -4.27620E-02
    SCF Iteration  27: Energy = -83.9951748874826905 dE = -2.35715E-02
    SCF Iteration  28: Energy = -84.0072057347110501 dE = -1.20308E-02
    SCF Iteration  29: Energy = -84.0135682844026519 dE = -6.36255E-03
    SCF Iteration  30: Energy = -84.0168148188785722 dE = -3.24653E-03
    SCF Iteration  31: Energy = -84.0185045683967644 dE = -1.68975E-03
    SCF Iteration  32: Energy = -84.0193690014646108 dE = -8.64433E-04
    SCF Iteration  33: Energy = -84.0198159980625974 dE = -4.46997E-04
    SCF Iteration  34: Energy = -84.0200451724136883 dE = -2.29174E-04
    SCF Iteration  35: Energy = -84.0201633394359817 dE = -1.18167E-04
    SCF Iteration  36: Energy = -84.0202240070863695 dE = -6.06677E-05
    SCF Iteration  37: Energy = -84.0202552468108337 dE = -3.12397E-05
    SCF Iteration  38: Energy = -84.0202712979690887 dE = -1.60512E-05
    SCF Iteration  39: Energy = -84.0202795578515378 dE = -8.25988E-06
    SCF Iteration  40: Energy = -84.0202838036089048 dE = -4.24576E-06
    SCF Iteration  41: Energy = -84.0202859877546189 dE = -2.18415E-06
    SCF Iteration  42: Energy = -84.0202871107002238 dE = -1.12295E-06
    SCF Iteration  43: Energy = -84.0202876882823375 dE = -5.77582E-07
    SCF Iteration  44: Energy = -84.0202879852716933 dE = -2.96989E-07
    SCF Iteration  45: Energy = -84.0202881380138535 dE = -1.52742E-07
    SCF Iteration  46: Energy = -84.0202882165577876 dE = -7.85439E-08
    SCF Iteration  47: Energy = -84.0202882569512468 dE = -4.03935E-08
    SCF Iteration  48: Energy = -84.0202882777232389 dE = -2.07720E-08
    SCF Iteration  49: Energy = -84.0202882884055100 dE = -1.06823E-08
    SCF Iteration  50: Energy = -84.0202882938989717 dE = -5.49346E-09
    Traceback (most recent call last):
      File "../code/scf.py", line 79, in <module>
        raise Exception("Maximum number of SCF iterations exceeded.")
    Exception: Maximum number of SCF iterations exceeded.

### Can we accelerate the SCF cycles?

<center>
<img src="scf-convergence-and-diis/slides/0007.jpg" alt="slide 7" width="800" height="600" />
<center/>
<center>
<img src="scf-convergence-and-diis/slides/0008.jpg" alt="slide 8" width="800" height="600" />
<center/>
<center>
<img src="scf-convergence-and-diis/slides/0009.jpg" alt="slide 9" width="800" height="600" />
<center/>
<center>
<img src="scf-convergence-and-diis/slides/0010.jpg" alt="slide 10" width="800" height="600" />
<center/>
<center>
<img src="scf-convergence-and-diis/slides/0011.jpg" alt="slide 11" width="800" height="600" />
<center/>

``` python
# ==> SCF Iterations <==
# Maximum SCF iterations
MAXITER = 50
# Energy convergence criterion
E_conv = 1.0e-10
D_conv = 1.0e-8
# Pre-iteration energy declarations
SCF_E = 0.0
E_old = 0.0
```

``` python
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
```

    There is a chance for diagonalisation

``` python
# Start from fresh orbitals
F_p =  np.dot(np.dot(X,H_core),X)
e, C_p = np.linalg.eigh(F_p)
C = np.dot(X,C_p)
C_occ = C[:, :ndocc]
D = np.einsum('pi,qi->pq', C_occ, C_occ, optimize=True)

# Trial & Residual Vector Lists
diis_iter = MAXITER # int(N) if you want to switch of DIIS after 'N' scf cycles 
F_list = []
DIIS_RESID = []
```

``` python
# ==> SCF Iterations w/ DIIS <==
print('==> Starting SCF Iterations <==\n')
# Begin Iterations
for scf_iter in range(1, MAXITER + 1):
    # Build Fock matrix
    J = np.einsum('pqrs,rs->pq', I, D, optimize=True)
    K = np.einsum('prqs,rs->pq', I, D, optimize=True)
    F = H_core + 2*J - K
    if scf_iter <= diis_iter:
        #print('DIIS is ON \n')
        diis_bool = True
        # Build DIIS Residual
        diis_r = X.dot(F.dot(D).dot(S) - S.dot(D).dot(F)).dot(X)
        # Append trial & residual vectors to lists
        F_list.append(F)
        DIIS_RESID.append(diis_r)
    else:
        #print('DIIS is OFF \n')
        diis_bool = False
    # Compute RHF energy
    SCF_E = np.einsum('pq,pq->', (H_core + F), D, optimize=True)
    dE = SCF_E - E_old
    dRMS = np.mean(diis_r**2)**0.5
    if diis_bool:
        diis_str = '   DIIS is ON'
    else:
        diis_str = '   DIIS is OFF'    
    print(('SCF Iteration %3d: Energy = %4.16f dE = % 1.5E dRMS = %1.5E' % (scf_iter, SCF_E, dE, dRMS)) + diis_str )
    
    # SCF Converged?
    if (abs(dE) < E_conv):
        break
    E_old = SCF_E
    
    if scf_iter >= 2:
        # Build B matrix
        B_dim = len(F_list) + 1
        B = np.empty((B_dim, B_dim))
        B[-1, :] = -1
        B[:, -1] = -1
        B[-1, -1] = 0
        for i in range(len(F_list)):
            for j in range(len(F_list)):
                B[i, j] = np.einsum('ij,ij->', DIIS_RESID[i], DIIS_RESID[j], optimize=True)

        # Build RHS of Pulay equation 
        rhs = np.zeros((B_dim))
        rhs[-1] = -1
        
        # Solve Pulay equation for c_i's with NumPy
        coeff = np.linalg.solve(B, rhs)
        
        # Build DIIS Fock matrix
        F = np.zeros_like(F)
        for x in range(coeff.shape[0] - 1):
            F += coeff[x] * F_list[x]
    
    # Compute new orbital guess with DIIS Fock matrix
    F_p =  X.dot(F).dot(X)
    e, C_p = np.linalg.eigh(F_p)
    C = X.dot(C_p)
    C_occ = C[:, :ndocc]
    D = np.einsum('pi,qi->pq', C_occ, C_occ, optimize=True)
    
    # MAXITER exceeded?
    if (scf_iter == MAXITER):
        raise Exception("Maximum number of SCF iterations exceeded.")

# Post iterations
print('\nSCF converged.')
print('Final RHF Energy: %.8f [Eh]' % SCF_E)
```

    ==> Starting SCF Iterations <==

    SCF Iteration   1: Energy = -69.1347968401194919 dE = -6.91348E+01 dRMS = 1.33446E-01   DIIS is ON
    SCF Iteration   2: Energy = -73.8555083348785502 dE = -4.72071E+00 dRMS = 5.03754E-02   DIIS is ON
    SCF Iteration   3: Energy = -80.9851657019508906 dE = -7.12966E+00 dRMS = 5.41654E-02   DIIS is ON
    SCF Iteration   4: Energy = -83.3476769740025816 dE = -2.36251E+00 dRMS = 2.29219E-02   DIIS is ON
    SCF Iteration   5: Energy = -84.0071321201823480 dE = -6.59455E-01 dRMS = 2.42631E-03   DIIS is ON
    SCF Iteration   6: Energy = -84.0191685957183410 dE = -1.20365E-02 dRMS = 9.18560E-04   DIIS is ON
    SCF Iteration   7: Energy = -84.0201944120005066 dE = -1.02582E-03 dRMS = 1.80133E-04   DIIS is ON
    SCF Iteration   8: Energy = -84.0202825456275946 dE = -8.81336E-05 dRMS = 3.73431E-05   DIIS is ON
    SCF Iteration   9: Energy = -84.0202880247193207 dE = -5.47909E-06 dRMS = 7.08676E-06   DIIS is ON
    SCF Iteration  10: Energy = -84.0202882912433040 dE = -2.66524E-07 dRMS = 1.35952E-06   DIIS is ON
    SCF Iteration  11: Energy = -84.0202882993160500 dE = -8.07275E-09 dRMS = 4.20633E-07   DIIS is ON
    SCF Iteration  12: Energy = -84.0202882997008658 dE = -3.84816E-10 dRMS = 6.97036E-08   DIIS is ON
    SCF Iteration  13: Energy = -84.0202882997147356 dE = -1.38698E-11 dRMS = 1.02648E-08   DIIS is ON

    SCF converged.
    Final RHF Energy: -84.02028830 [Eh]

``` python
# Dataset of energies for H2O molecules computed using psi4 with different basis.
import pandas as pd
df = pd.read_csv('../data/energies_benchmark.csv', index_col='basis')
print(df)
```

                 nbf        scf        mp2       ccsd    ccsd(t)       pbe0  \
    basis                                                                     
    sto-3g         7 -82.944446 -82.993582 -83.015126 -83.015226 -83.243909   
    sto-6g         7 -83.659154 -83.708632 -83.730426 -83.730527 -83.963042   
    3-21g         13 -83.563679 -83.696978 -83.705981 -83.708087 -83.880180   
    6-31g         13 -83.954896 -84.097007 -84.104308 -84.105907 -84.283238   
    6-31g*        19 -83.977115 -84.177677 -84.187182 -84.189913 -84.302479   
    6-31+g*       23 -83.984260 -84.191978 -84.199963 -84.203457 -84.314272   
    cc-pvdz       24 -83.992162 -84.206491 -84.216072 -84.219958 -84.315199   
    cc-pvtz       58 -84.020288 -84.305513 -84.310393 -84.319489 -84.348190   
    cc-pvqz      115 -84.027569 -84.351165 -84.353690 -84.364479 -84.357060   
    aug-cc-pvdz   41 -84.005721 -84.240482 -84.247321 -84.253774 -84.335114   
    aug-cc-pvtz   92 -84.024000 -84.318366 -84.321684 -84.332038 -84.354125   
    aug-cc-pvqz  172 -84.028840 -84.356530 -84.357987 -84.369280 -84.359510   
    def2-sv(p)    18 -83.907364 -84.108361 -84.119167 -84.121896 -84.235178   
    def2-svp      24 -83.925270 -84.138846 -84.148504 -84.152310 -84.251446   

                     b3lyp  
    basis                   
    sto-3g      -83.314658  
    sto-6g      -84.037073  
    3-21g       -83.963402  
    6-31g       -84.369483  
    6-31g*      -84.388000  
    6-31+g*     -84.401838  
    cc-pvdz     -84.399149  
    cc-pvtz     -84.435662  
    cc-pvqz     -84.445071  
    aug-cc-pvdz -84.421785  
    aug-cc-pvtz -84.442275  
    aug-cc-pvqz -84.447817  
    def2-sv(p)  -84.320335  
    def2-svp    -84.335787  

**References** 1. A. Szabo and N. S. Ostlund, *Modern Quantum
Chemistry*, Introduction to Advanced Electronic Structure Theory.
Courier Corporation, 1996. 2. I. N. Levine, *Quantum Chemistry*.
Prentice-Hall, New Jersey, 5th edition, 2000. 3. P. Pulay. *Chem. Phys.
Lett.* **73**, 393-398 (1980) 4. C. David Sherrill. *“[Some comments on
accellerating convergence of iterative sequences using direct inversion
of the iterative subspace
(DIIS)](vergil.chemistry.gatech.edu/notes/diis/diis.pdf)”.* (1998) 5.
Smith, Daniel GA, et al. “Psi4NumPy: An interactive quantum chemistry
programming environment for reference implementations and rapid
development.” Journal of chemical theory and computation 14.7 (2018):
3504-3511.

<center>
<img src="scf-convergence-and-diis/slides/0012.jpg" alt="slide 12" width="800" height="600" />
<center/>
