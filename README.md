# Protein Ion Coordination Survey

[![Join the chat at https://gitter.im/Becksteinlab/PDB_Ion_Survey](https://badges.gitter.im/Becksteinlab/PDB_Ion_Survey.svg)](https://gitter.im/Becksteinlab/PDB_Ion_Survey?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

As part of the [BioXFEL Ed Scientific Villages at ASU](https://www.bioxfel.org/education/high-school-learning-communities) program, [@kaceyreidy](https://github.com/kaceyreidy) is building a reusable data acquisition and analysis pipeline to characterize the interactions between ions and proteins. 

The input data are crystal structures from the [Protein Databank](http://www.pdb.org). [MDAnalysis](http://www.mdanalysis.org) is used to analyze the data. The initial focus is on a statistical analysis of the coordination numbers of ions, in particular oxygen-cation interactions.

# How to Use

* Use get_proteins to search the PDB for molecules of a certain type with a specified ion.
* Use get_pdb_file to download pdb files from the PDB.
* Use gee to create radial distributions for oxygen atoms around an ion.
* Use ofr to create a cumulative radial distribution function for oxygen atoms around an ion.
* Use gofr to use both gee and ofr.
* Talk to Martin Blech to use unparse or _emit.
