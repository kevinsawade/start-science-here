 .. _gromacs-production-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

=============
Production MD
=============

Upon completion of the two equilibration phases, the system is now well-equilibrated at the desired temperature and pressure. We are now ready to release the position restraints and run production MD for data collection. The process is just like we have seen before, as we will make use of the checkpoint file (which in this case now contains preserve pressure coupling information) to grompp. We will run a 1-ns MD simulation, the script for which can be found `here <http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp>`_.

.. code-block:: bash

   $ wget http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp

Intermission
============

The production MD will probably not run in binder. We have to find a way to run the productioon MD on an actual computer and not on the cloud via binder. Here are the steps to achieve this:

1. Download a .zip archive with our current simulation.
2. Get to a linux terminal via using a computer with linux, or the Windows subsystem for linux: :ref:`wsl-label`.
3. Install gromacs on that computer. Continue the run.
4. Upload files to WSL.

Download a .zip archive with the current simulation
---------------------------------------------------

Go to the directory ``~/home/jovyan`` in binder and use ``ls`` to make sure, that the directory you've ran the sims in (``gromacs_tutorial/`` or similar) is present. Otherwise browse the filesystem to find that directory. Use:

.. code-block:: bash

   $ zip -r gromacs_tutorial.zip gromacs_tutorial/

The ``-r`` flag tells zip to zip recursively (all sub-directories in ``gromacs_tutorial/``). Check the filesize of the zip archive with:

.. code-block:: bash

   $ du -h gromacs_tutorial.zip

Download the .zip archive on binder: :ref:`downloading-files-label`.

Get to a linux terminal
-----------------------

If your computer runs any flavor of linux, you can skip this step, if you're running Windows, please use the guide over at :ref:`wsl-label` to install the Windows subsystem for linux.

Install Gromacs
---------------

Once you are greeted by the bash prompt: ``username@hostname:~$`` make sure, that you are running ubuntu, via

.. code-block ::

   $ cat /etc/*release

If you are running Ubuntu, you can use ``sudo apt`` to install gromacs.

.. code-block::

   $ sudo apt install gromacs

If you are running any other flavor of linux, please refer to the the note below:

.. note::
   :class: dropdown

   Please refer to the man pages of your package manager. Depending on what package manager you use you might want to use one of these commands:

   .. code-block::

      $ brew install gromacs # Mac OS X
      $ yum install gromacs # RHEL
      $ pacman install gromacs # Arch


Upload files to WSL
-------------------

Because the WSL is somewhat buried inside the Windows system, we have to take some extra steps to put our .zip archive in our home directory. Inside WSL type

.. code-block:: bash

   $ explorer.exe .

A window of the Windows explorer will open and you can put your .zip archive there. Last step is to unzip it via:

.. code-block:: bash

   $ unzip gromacs_turorial.zip


Continuation
============

.. code-block:: bash

   $ gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

grompp will print an estimate for PME load, which will dictate how many processors should be dedicated to the PME calculation, and how many for the PP calculations. Refer to the GROMACS 4 `publication <http://dx.doi.org/10.1021/ct700301q>`_ and the manual for details.

.. code-block::

   Estimate for the relative computational load of the PME mesh part: 0.22

For a cubic box, the optimal setup will have a PME load of 0.25 (3:1 PP:PME - we're very close to optimal!); for a dodecahedral box, the optimal PME load is 0.33 (2:1 PP:PME). When executing mdrun, the program should automatically determine the best number of processors to assign for the PP and PME calculations. Thus, make sure you indicate an appropriate number of threads/cores for your calculation (the value of -nt X), so that you can get the best performance.

Now, execute mdrun:

.. code-block:: bash

   $ gmx mdrun -deffnm md_0_1

In GROMACS 2018, the PME calculations can be offloaded to graphical processing units (GPU), which speeds up the simulation substantially. Using a Titan Xp GPU, this system can be simulated at an astounding 295 ns/day!

Running GROMACS on GPU
======================

As of version 4.6, GROMACS supports the use of GPU accelerators for running MD simulations. With the release of version 2018, the nonbonded interactions and PME are calculated on the GPU, with only bonded forces calculated on the CPU cores. When building GROMACS (see www.gromacs.org for installation instructions), GPU hardware will automatically be detected, if present. The minimum requirements for using GPU acceleration are the CUDA libraries and SDK, and a GPU with a compute capability of >= 2.0. A nice list of some of the more common GPUs and their specifications can be found `here <https://developer.nvidia.com/cuda-gpus>`_.

Assuming you have one GPU available, the mdrun command to make use of it is as simple as:

.. code-block:: bash

   $ gmx mdrun -deffnm md_0_1 -nb gpu

If you have more than one GPU available, or require customization of how the work is divided up via the hybrid parallelization scheme available in GROMACS, please consult the GROMACS manual and webpage. Such technical details are beyond the scope of this tutorial.