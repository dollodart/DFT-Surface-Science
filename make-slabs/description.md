This code creates supercells for slabs with varying thickness and a vacuum spacing. The python file generates the coordinates and the bash executive file runs the job and assembles the coordinates and information into the POSCAR.vasp format. This slab generation could be easily achieved for ZnO because, though it is is non-centrosymmetric and requires the deletion of physically meaningless surface terminations, the surfaces of interest are normal to the crystallographic replication vectors. It may be of some interest to make a slab generator for any Cartesian surface unit cell, in which the coordinates of a surface unit cell with perpendicular replications is made from a general surface unit cell. It may also be of interest to modify the code to allow for the inclusion of adsorbates. These additions may come in the future.

It is now possible to place adsorbates at the bottom and top surface.
Run the adsorbate python script before running the bash script.
The advantage is not necessarily in making the exact POSCAR file with the adsorbates positions.
One can manually edit positions to obtain lower coverage by removing entries, for example (remember to update the header containing atom numbers). 

Note that the non-polar case has spurious odd numbered surfaces from borrowing the polar script. This may be edited in the future.

Desired future features: library of adsorbate positions for hole, edge, and site adsorption.
