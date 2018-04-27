This is a more generalized sublattice method to make POSCAR supercells for input to VASP.
In the example an orthogonal supercell of magnesium oxide is replicated to make a slab with vacuum spacing.
One may account for effects of non-centrosymmetry in a crystal such as Zinc Oxide by adjusting the basis and extent of replication for each sublattice appropriately.
That is, this method requires only modification of the scr.py file for any configuration of any number of atoms.
For example, an adsorbate atom can be placed without any replication by giving replication indices of (1,1,1) and a lattice vector set of 0s, and giving the basis as the desired position.
Any number of species, including redundant species as the MgO example demonstrates, may be input. The supercell of calculation can differ from the supercell (or unit or primitive cell) lattice vectors provided, as is done here to make a slab.
In premise, this sublattice method is sufficient to make any geometry, because if one knows the basis adjustment and extent of replication needed, one may include the number of atoms required.
In practice, it is difficult to predict the required basis adjustment and extent of replication (see the scripts of Efficient creation and convergence of surface slabs), and the orthogonal method is simple to implement, though as the paper notes, it is more inefficient. 
