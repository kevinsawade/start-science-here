.. start-science-here documentation master file, created by
   sphinx-quickstart on Sun Jul 25 08:57:53 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================
Start Science Here!
===================

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/kevinsawade/start-science-here/HEAD?urlpath=%2Ftree%2F
   :alt: Binder

.. image:: https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/kevinsawade/bcd9d09bc682b4743b84fc6e967478ac/raw/endpoint.json
   :target: https://www.chemie.uni-konstanz.de/ag-peter/
   :alt: MadeWithLove

**An introduction to scientific programming, Python, molecular dynamics, and more**

Welcome to start science here! This project is meant to introduce you to the world of scientific programming and teach you the most important concepts. This project is subdivided into different main topics. You can dive right into the Python tutorial by visiting: https://mybinder.org/v2/gh/kevinsawade/start-science-here/HEAD?urlpath=%2Ftree%2F. If you want to take the full course, take the recommended route:

1. Linux tutorial
2. Windows Subsystem for Linux (Windows users only)
3. Installing Docker
4. Gromacs
5. Python
6. Setting up your own Python
7. Git
8. Your first project

There is no set up required to follow along with this tutorial. It can take three to four days to work through everything (the Python tutorial is very long, but can be put on hold after the beginner tutorials).

.. prereq::

   * A working internet connection
   * An up-to-date browser (Chrome, Firefox, Safari, Opera, Edge)

0. Proposing Changes and Issues with start science here
=======================================================

Before you leave this page and start the tutorials check out the pages :ref:`raising-issues-label` and :ref:`propose-changes-label`.

We are trying to continually improve start science here and are more than happy if users bring us feedback. These two pages allow you to give us such feedback. Be it spelling mistakes, factual mistakes or bugs. You can also raise issues on GitHub if you run into any problems that you can't fix yourself.

1. Linux tutorial
=================

Start the recommended route with the linux tutorial here: :ref:`linux-tutorial-label`.

2. Windows Subsystem for Linux (Windows users only)
===================================================

You now have a pretty solid understanding about unix-like operating systems and know how to work with a command-line interface. Binder makes it very easy to learn because all you need is a browser and a working internet connection. However, we can't use binder forever, at some point we will want to switch to our own computers. The caveat is that there are three prevalent OSes. Mac OS X and our favorite Linux distribution (Ubuntu) are very similar, but Windows 10 is very different. Luckily, there is the Windows subsystem for Linux, which gives you a working Ubuntu machine for your Windows 10 PC. Here's a summary of how to install the latest version of Ubuntu: :ref:`wsl-label`.

3. Installing Docker
====================

Docker allows to bring you pre-configured ubuntu installations. Instead of writing pages upon pages of installation instructions for all kinds of programs, we will install docker and can provide you different environments by calling:

.. code-block:: bash

   $ docker exec -it kevinsawade:gromacs bash

For a working gromacs installation, or

.. code-block:: bash

   $ docker exec -it kevinsawade:conflict-resolution

Although some gromacs exercises can be done on binder, docker is much easier and faster for larger computations. Follow these instructions to install it: :ref:`install-docker-label`.

4. Gromacs
==========

Gromacs is a simulation package for MD simulations of proteins, DNA, and lipids. It is not only fast, utilizing CPUs and GPUs, but also flexible and extensible. We use it for all-atom and coarse-grained simulations of proteins and DNA, but also for advanced sampling techniques like replica exchange and path sampling. There are great tutorials for gromacs available at: http://www.mdtutorials.com/gmx/

I adjusted the first tutorial "Lysozyme in water" to work with binder. So follow along at: :ref:`gromacs-tutorial-label`.

5. Python
=========

The longest part of Start Science Here! is the Python tutorial. If you feel confident with Python you might just glimpse over the quick reference at :ref:`python-tutorial-label`. Otherwise you can open this Repo in binder and start making your way through the notebooks.

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   linux_tutorial/linux_tutorial
   windows_susbsystem_for_linux
   install_docker
   gromacs_tutorial/gromacs_tutorial
   python_tutorial/python_tutorial
   setting_up_python
   git_tutorial
   first_project



.. toctree::
   :maxdepth: 1
   :caption: Extras

   raising_issues
   propose_changes
   downloading_files
   uploading_files

.. toctree::
   :maxdepth: 1
   :caption: Having Problems?

   troubleshooting

6. Setting up your own Python
=============================

7. Git
======

8. Your first project
=====================

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
