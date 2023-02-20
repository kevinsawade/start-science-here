 .. _gromacs-energy-minimization-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

===================
Energy Minimization
===================

The solvated, electroneutral system is now assembled. Before we can begin dynamics, we must ensure that the system has no steric clashes or inappropriate geometry. The structure is relaxed through a process called energy minimization (EM).

The process for EM is much like the addition of ions. We are once again going to use grompp to assemble the structure, topology, and simulation parameters into a binary input file (.tpr), but this time, instead of passing the .tpr to genion, we will run the energy minimization through the GROMACS MD engine, mdrun.

Assemble the binary input using grompp using `this <http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp>`_ input parameter file, which you can get with ``wget``.

.. code-block:: bash

   $ wget http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp

.. code-block:: bash

   $ gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr

Make sure you have been updating your topol.top file when running genbox and genion, or else you will get lots of nasty error messages ("number of coordinates in coordinate file does not match topology," etc).

We are now ready to invoke mdrun to carry out the EM:

.. code-block:: bash

   $ gmx mdrun -v -deffnm em

The -v flag is for the impatient: it makes mdrun verbose, such that it prints its progress to the screen at every step. The -deffnm flag will define the file names of the input and output. So, if you did not name your grompp output "em.tpr," you will have to explicitly specify its name with the mdrun -s flag. In our case, we will get the following files:

* em.log: ASCII-text log file of the EM process
* em.edr: Binary energy file
* em.trr: Binary full-precision trajectory
* em.gro: Energy-minimized structure

There are two very important factors to evaluate to determine if EM was successful. The first is the potential energy (printed at the end of the EM process, even without -v). E\ :sub:`pot` should be negative, and (for a simple protein in water) on the order of 10\ :sup:`5`-10\ :sup:`6`, depending on the system size and number of water molecules. The second important feature is the maximum force, F\ :sub:`max`, the target for which was set in minim.mdp - "emtol = 1000.0" - indicating a target F\ :sub:`max` of no greater than 1000 kJ mol\ :sup:`-1` nm\ :sup:`-1`. It is possible to arrive at a reasonable Epot with F\ :sub:`max` > emtol. If this happens, your system may not be stable enough for simulation. Evaluate why it may be happening, and perhaps change your minimization parameters (integrator, emstep, etc).

Let's do a bit of analysis. The em.edr file contains all of the energy terms that GROMACS collects during EM. You can analyze any .edr file using the GROMACS energy module:

.. code-block:: bash

   $ gmx energy -f em.edr -o potential.xvg

At the prompt, type "10 0" to select Potential (10); zero (0) terminates input. You will be shown the average of Epot, and a file called "potential.xvg" will be written. To plot this data, you will need to call one of the provided scripts like so:

.. code-block:: bash

   $ python3 plot_xvg.py potential.xvg

.. sshsolution::
   :class: dropdown

   Your data will now be plotted to your terminal.

   .. image:: ../_static/pics/gromacs_tutorial/plotted_to_terminal.png
      :target: ../_static/pics/gromacs_tutorial/plotted_to_terminal.png
      :alt: Decay of potential energy during minimization.

Now that our system is at an energy minimum, we can begin real dynamics.
