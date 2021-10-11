 .. _gromacs-tutorial-generate-topology-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

=================
Generate Topology
=================

We must download the protein structure file with which we will be working. For this tutorial, we will utilize hen egg white lysozyme (PDB code 1AKI). Open a Terminal on Binder and create a new directory:

.. code-block:: bash

   $ cd gromacs_tutorial

Currently the directory is empty, except the ``plot_xvg.py`` script, that will help us plot our results. Have a look at it with:

.. code-block:: bash

   $ ls
   $ python3 plot_xvg.py --help:

Download the 1AKI.pdb file from the `RCSB <http://www.rcsb.org/pdb/home/home.do>`_ via:

.. code-block:: bash

   $ wget https://files.rcsb.org/view/1AKI.pdb -O 1AKI.pdb

Have a look at the file with

.. code-block:: bash

   $ less 1AKI.pdb

Once you have downloaded the structure, you can visualize the structure using a viewing program such as VMD, Chimera, PyMOL, etc. As we are limited by binder, we will use NGLView to visualize the structure in a jupyter-notebook. For that select New and then Python 3 in binder. Paste this code into the notebook and press :boldcode:`[Shift]` + :boldcode:`[Enter]` to run the code. Alternatively, you can click the *run* button.

.. code-block:: python

   import nglview as ngl
   view = ngl.show_file('gromacs_tutorial/1AKI.pdb', gui=True)
   view.clear_representations()
   view.add_ball_and_stick()
   view

.. image:: ../_static/pics/gromacs_tutorial/lysozyme_nglview.png
   :target: ../_static/pics/gromacs_tutorial/lysozyme_nglview.png
   :alt: Ball and Stick and Cartoon representation of the Lysozyme protein.

Once you've had a look at the molecule, you are going to want to strip out the crystal waters. To delete the water molecules (residue "HOH" in the PDB file), either use a plain text editor like vi, emacs (Linux/Mac), or Notepad (Windows). Do not use word processing software! Alternatively, you can use grep to delete these lines very easily:

.. code-block:: bash

   $ grep -v HOH 1aki.pdb > 1AKI_clean.pdb

Note that such a procedure is **not universally appropriate** (e.g., the case of a tightly bound or otherwise functional active-site water molecule). For our intentions here, we do not need crystal water.

Always check your .pdb file for entries listed under the comment MISSING, as these entries indicate either atoms or whole residues that are not present in the crystal structure. Terminal regions may be absent, and may not present a problem for dynamics. Incomplete internal sequences or any amino acid residues that have missing atoms will cause pdb2gmx to fail. These missing atoms/residues must be modeled in using other software packages. Also note that pdb2gmx is not magic. It cannot generate topologies for arbitrary molecules, just the residues defined by the force field (in the *.rtp files - generally proteins, nucleic acids, and a **very** finite amount of cofactors, like NAD(H) and ATP).

Now that the crystal waters are gone and we have verified that all the necessary atoms are present, the PDB file should contain only protein atoms, and is ready to be input into the first GROMACS module, pdb2gmx. The purpose of pdb2gmx is to generate three files:

1. The topology for the molecule.
2. A position restraint file.
3. A post-processed structure file.

The topology (topol.top by default) contains all the information necessary to define the molecule within a simulation. This information includes nonbonded parameters (atom types and charges) as well as bonded parameters (bonds, angles, and dihedrals). We will take a more detailed look at the topology once it has been generated.

Execute pdb2gmx by issuing the following command:

.. code-block:: bash

   $ gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce

The structure will be processed by pdb2gmx, and you will be prompted to choose a force field:

.. code-block::

   Select the Force Field:
    From '/usr/local/gromacs/share/gromacs/top':
     1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)
     2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)
     3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)
     4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)
     5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)
     6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
     7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)
     8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)
     9: GROMOS96 43a1 force field
    10: GROMOS96 43a2 force field (improved alkane dihedrals)
    11: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)
    12: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)
    13: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)
    14: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)
    15: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)

The force field will contain the information that will be written to the topology. This is a very important choice! You should always read thoroughly about each force field and decide which is most applicable to your situation. For this tutorial, we will use the all-atom OPLS force field, so type 15 at the command prompt, followed by 'Enter'.

There are many other options that can be passed to pdb2gmx. Some commonly used ones are listed here:

* -ignh: Ignore H atoms in the PDB file; especially useful for NMR structures. Otherwise, if H atoms are present, they must be in the named exactly how the force fields in GROMACS expect them to be. Different conventions exist, so dealing with H atoms can occasionally be a headache! If you need to preserve the initial H coordinates, but renaming is required, then the Linux sed command is your friend.
* -ter: Interactively assign charge states for N- and C-termini.
* -inter: Interactively assign charge states for Glu, Asp, Lys, Arg, and His; choose which Cys are involved in disulfide bonds.

You have now generated three new files: 1AKI_processed.gro, topol.top, and posre.itp. 1AKI_processed.gro is a GROMACS-formatted structure file that contains all the atoms defined within the force field (i.e., H atoms have been added to the amino acids in the protein). The topol.top file is the system topology (more on this in a minute). The posre.itp file contains information used to restrain the positions of heavy atoms (more on this later).

One final note: many users assume that a .gro file is mandatory. **This is not true**. GROMACS can handle many different file formats, with .gro simply being the default for commands that write coordinate files. It is a very compact format, but it has limited precision. If you prefer to use, for instance, PDB format, all you need to do is to specify an appropriate file name with .pdb extension as your output. The purpose of pdb2gmx is to produce a force field-compliant topology; the output structure is largely a side effect of this purpose and is intended for user convenience. The format can be just about anything you like (see the GROMACS manual for different formats).

Next step: :ref:`gromacs-tutorial-examine-topology-label`
