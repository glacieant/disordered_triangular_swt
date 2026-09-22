# disordered_triangular_swt #

The C program to compute non-linear corrections to the distorted
moment created by a bond defect in the semi-classical 120-degree
ground state of the triangular lattice Heisenberg antiferromagnet.
These computations were inspired by the site-defect work by Wollny
et al. [1] and contributed to the larger bond-disorder program
studied by us [2].

Key features:

1. A real-space spin-wave theory computation using the generalized
   bosonic Bogoliubov transformation [3,4].

2. A saddle-point treatment of the non-linear spin-wave
   correction resulting from the defect-induced translation
   symmetry breaking; the correction enters through a cubic
   term and is therefore exact to leading order.

3. Calculation of various physical observables and the spin structure
   factor in reciprocal space.


## References ##

1. A. Wollny, L. Fritz, and M. Vojta, Fractional impurity moments in
two-dimensional noncollinear magnets, Phys. Rev. Lett. 107, 137204
(2011), https://doi.org/10.1103/PhysRevLett.107.137204

2. S. Dey, E. C. Andrade, and M. Vojta, Destruction of long-range order
in noncollinear two-dimensional antiferromagnets by random-bond
disorder, Phys. Rev. B 101, 020411(R) (2020),
https://doi.org/10.1103/PhysRevB.101.020411

3. J. H. P. Colpa, Diagonalization of the quadratic boson Hamiltonian,
Physica A 93, 327 (1978),
https://doi.org/10.1016/0378-4371(78)90160-7

4. S. Wessel and I. Milat, Quantum fluctuations and excitations in
antiferromagnetic quasicrystals, Phys. Rev. B 71, 104427 (2005),
https://doi.org/10.1103/PhysRevB.71.104427
