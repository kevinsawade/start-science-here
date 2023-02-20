 .. _gromacs-analysis-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

========
Analysis
========

Now that we have simulated our protein, we should run some analysis on the system. What types of data are important? This is an important question to ask before running the simulation, so you should have some ideas about the types of data you will want to collect in your own systems. For this tutorial, a few basic tools will be introduced.

The first is trjconv, which is used as a post-processing tool to strip out coordinates, correct for periodicity, or manually alter the trajectory (time units, frame frequency, etc). For this exercise, we will use trjconv to account for any periodicity in the system. The protein will diffuse through the unit cell, and may appear "broken" or may "jump" across to the other side of the box. To account for such behaviors, issue the following:

.. code-block:: bash

   $ gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center

Select 1 ("Protein") as the group to be centered and 0 ("System") for output. We will conduct all our analyses on this "corrected" trajectory. Let's look at structural stability first. GROMACS has a built-in utility for RMSD calculations called rms. To use rms, issue this command:

.. code-block:: bash

   $ gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns

Choose 4 ("Backbone") for both the least-squares fit and the group for RMSD calculation. The -tu flag will output the results in terms of ns, even though the trajectory was written in ps. This is done for clarity of the output (especially if you have a long simulation - 1e+05 ps does not look as nice as 100 ns). The output plot will show the RMSD relative to the structure present in the minimized, equilibrated system:

If we wish to calculate RMSD relative to the crystal structure, we could issue the following:

.. code-block:: bash

   $ gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns

Use:

.. code-block::

   $ python3 plot_xvg.py rmsd.xvg rmsd_xtal.xvg -lbl ref:eq ref:crystal

To plot both xvg files and manually set the labels

.. sshsolution::
   :class: dropdown

   .. image:: ../_static/pics/gromacs_tutorial/plotted_analysis.png
      :target: ../_static/pics/gromacs_tutorial/plotted_analysis.png
      :alt: Difference between reference equilibrated and reference crystal is that RMSD distance to crystal is generally greater.

Both time series show the RMSD levels off to ~0.1 nm (1 Å), indicating that the structure is very stable. Subtle differences between the plots indicate that the structure at t = 0 ns is slightly different from this crystal structure. This is to be expected, since it has been energy-minimized, and because the position restraints are not 100% perfect, as discussed previously.
