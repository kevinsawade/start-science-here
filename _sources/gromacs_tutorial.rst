 .. _gromacs-tutorial-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

================
Gromacs Tutorial
================

Credit
======

Thanks to Justin A. Lemkul, Ph.D. (Virginia Tech Department of Biochemistry) for his GROMACS tutorials, I adjusted them a bit to include ``wget`` and use a plotting tool for the command line.

Some Gromacs Basics
===================

With the release of version 5.0 of GROMACS, all of the tools are essentially modules of a binary named "gmx". This is a departure from previous versions, wherein each of the tools was invoked as its own command. In 5.0, this still works via symlinks, but they will go away in future versions, so it is best to get accustomed to the new way of doing things. To get help information about any GROMACS module, you can invoke either of the following commands:

.. code-block:: bash

   $ gmx help (module)

or

.. code-block:: bash

   $ gmx (module) -h

where (module) is replaced by the actual name of the command you're trying to issue. Information will be printed to the terminal, including available algorithms, options, required file formats, known bugs and limitations, etc. For new users of GROMACS, invoking the help information for common commands is a great way to learn about what each command can do. Now, on to the fun stuff!


Lysozyme Tutorial
=================

Start the tutorial: :ref:`gromacs-tutorial-generate-topology-label`

.. toctree::
   :maxdepth: 1

   gromacs_tutorial_generate_topology
   gromacs_tutorial_examine_topology
   gromacs_tutorial_solvation
   gromacs_tutorial_adding_ions
   gromacs_tutorial_energy_minimization
   gromacs_tutorial_equilibration
   gromacs_tutorial_equilibration_2
   gromacs_tutorial_production
   gromacs_tutorial_analysis
   gromacs_tutorial_more_analysis
