 .. _gromacs-tutorial-adding-ions-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

===========
Adding Ions
===========

We now have a solvated system that contains a charged protein. The output of pdb2gmx told us that the protein has a net charge of +8e (based on its amino acid composition). If you missed this information in the pdb2gmx output, look at the last line of your ``[ atoms ]`` directive in topol.top; it should read (in part) "qtot 8." Since life does not exist at a net charge, we must add ions to our system.

The tool for adding ions within GROMACS is called genion. What genion does is read through the topology and replace water molecules with the ions that the user specifies. The input is called a run input file, which has an extension of .tpr; this file is produced by the GROMACS grompp module (GROMACS pre-processor), which will also be used later when we run our first simulation. What grompp does is process the coordinate file and topology (which describes the molecules) to generate an atomic-level input (.tpr). The .tpr file contains all the parameters for all of the atoms in the system.

To produce a .tpr file with grompp, we will need an additional input file, with the extension .mdp (molecular dynamics parameter file); grompp will assemble the parameters specified in the .mdp file with the coordinates and topology information to generate a .tpr file.

An .mdp file is normally used to run energy minimization or an MD simulation, but in this case is simply used to generate an atomic description of the system. An example .mdp file (the one we will use) can be downloaded `here <http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp>`_. Download that file to your PC and upload it to binder :ref:`uploading-files-label`.

Alternatively you can call this command:

.. code-block:: bash

   $ wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp

And move the file to the correct ``~/gromacs_tutorial`` directory.

In reality, the .mdp file used at this step can contain any legitimate combination of parameters. I typically use an energy-minimization script, because they are very basic and do not involve any complicated parameter combinations.

.. note::

   The files provided with this tutorial are intended **only** for use with the OPLS-AA force field. Settings, particularly non-bonded interaction settings, will be different for other force fields.

Assemble your .tpr file with the following:

.. code-block:: bash

   $ gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr

Now we have an atomic-level description of our system in the binary file ions.tpr. We will pass this file to genion:

.. code-block:: bash

   $ gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

When prompted, choose group 13 "SOL" for embedding ions. You do not want to replace parts of your protein with ions.

In the genion command, we provide the structure/state file (-s) as input, generate a .gro file as output (-o), process the topology (-p) to reflect the removal of water molecules and addition of ions, define positive and negative ion names (-pname and -nname, respectively), and tell genion to add only the ions necessary to neutralize the net charge on the protein by adding the correct number of negative ions (-neutral, which in this case will add 8 Cl- ions to offset the +8 charge on the protein). You can also use genion to add a specified concentration of ions in addition to simply neutralizing the system by specifying the -neutral and -conc options in conjunction. Refer to the genion man page for information on how to use these options.

The names of the ions specified with -pname and -nname were force field-specific in previous versions of GROMACS, but were standardized in version 4.5. The specified ion names are always the elemental symbol in all capital letters, which is the ``[ moleculetype ]`` name that is then written to the topology. Residue or atom names may or may not append the sign of the charge (+/-), depending on the force field. **Do not use atom or residue names in the genion command, or you will encounter errors in subsequent steps.**

Your ``[ molecules ]`` directive should now look like:

.. code-block::

   [ molecules ]
   ; Compound      #mols
   Protein_A         1
   SOL           10636
   CL                8

Next step: :ref:`gromacs-energy-minimization-label`
