============
weak-lensing
============

Implement KSB methods for determining cosmic shear and confirm the emergence of bias. Explore Bayesian methods of inferring shear.



toy.py
======


uncorrelated galaxies for a given shear
---------------------------------------
Do not use the -p option 
E.g.
python toy.py -n <number of galaxies=3> -1 <shear component 1=-0.01> -2 <shear component 2=0.02> -e <sigma_e=0.05> -s <sigma_pr=0.3> -t <integrator tolerance=1.49e-06> -i <file index=1> 
Add -d if you want to display average time per galaxy

This will save a file "g<TIME>i<index>.csv" to "data/"
where <TIME> is a time stamp and <index> is the index you specified above. Make sure you provide unique indices to each cluster job so that the jobs don't write over each other and make a mess.

The CSV files are comma separated with the format P,Q1,Q2,R11,R12,R22.



correlated pairs
----------------

Use the -p option but don't specify -1 or -2
E.g.
python toy.py -p -n <number of galaxies=3> -e <sigma_e=0.05> -s <sigma_pr=0.3> -t <integrator tolerance=1.49e-06> -i <file index=1> 

This will draw n pairs (g,h) from generate_pairs.py.
It will calculate P,Q,R for these pairs.

And it will save files "g<TIME>i<index>.csv" and "h<TIME>i<index>.csv" to "data/"
where <TIME> is a time stamp and <index> is the index you specified above. Make sure you provide unique indices to each cluster job so that the jobs don't write over each other and make a mess.

The CSV files are comma separated with the format P,Q1,Q2,R11,R12,R22.


postprocess.py (for uncorrelated galaxies of a given shear)
===========================================================

This will process every CSV file in a specified directory assuming they contain P,Q,R for uncorrelated galaxies of a specified shear.
And it will infer that shear.
No command line options yet.


pairs_postprocess.py (for correlated pairs)
===========================================

This will process every CSV file in a specified directory assuming they contain P,Q,R for correlated pairs.
And it will try to infer the covmat.
No command line options yet.
