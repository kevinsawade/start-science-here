.. start-science-here documentation master file, created by
   sphinx-quickstart on Sun Jul 25 08:57:53 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================
Start Science Here!
===================

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/kevinsawade/start-science-here/HEAD
   :alt: Binder


.. image:: https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/kevinsawade/bcd9d09bc682b4743b84fc6e967478ac/raw/endpoint.json
   :target: https://www.chemie.uni-konstanz.de/ag-peter/
   :alt: MadeWithLove

**An introduction to scientific programming, python, molecular dynamics and more**

Welcome to start science here! This project is meant to introduce you to the world of scientific programming and teach you the most important concepts. This project is subdivided into different main topics. You can dive right into the python tutorial by visiting: https://mybinder.org/v2/gh/kevinsawade/start-science-here/HEAD. If you want to take the full course, take the recommended route:

1. Working with Linux from the command line
2. Installing the Windows subsystem for Linux (Windows users only)
3. Installing Docker
4. Gromacs
5. Python
6. Setting up your own python
7. Git
8. Your first project

There is no set up required to follow along with this tutorial. It can take three to four days to work through everything (the python tutorial is very long, but can be put on hold after the beginner tutorials).

.. prereq::

   * A working internet connection
   * An up-to-date browser (Chrome, Firefox, Safari, Opera, Edge)

Linux Tutorial
==============

The recommended route starts with the linux tutorial: :ref:`linux-tutorial-label`.

Windows subsystem for Linux
===========================

You now have a pretty solid understanding about unix-like operating systems and how to work with a command-line interface. Binder makes it very easy to learn because all you need is a working browser and internet connection. However, we can't use binder forever, at some point we want to switch to our onw compters. The caveat is, that there are three prevalent OSes. Mac OS X and our favorite Linux distribution (ubuntu) are very similar, but Windows 10 is very different. Luckily, there is the Windows subsystem for linux, which gives you a working ubuntu machine for your Windows 10 PC. Here's a summary of how to install the latest version of ubuntu: :ref:`wsl-label`.

Docker
======

Docker allows us to bring you pre-configured ubuntu installations. Instead of writing pages upon pages of installation instructions for all kinds of pograms, we will install docker and can bring you different environments by calling:

.. code-block:: bash

   $ docker exec -it kevinsawade:gromacs bash

For a working gromacs installation, or

.. code-block:: bash

   $ docker exec -it kevinsawade:merge_conflict

Although some gromacs exercises can be done on binder, docker is much easier and faster for larger computations. Follow these instructions to install it: :ref:`install-docker-label`.

.. toctree::
   :maxdepth: 1
   :hidden:

   linux_tutorial
   windows_susbsystem_for_linux
   install_docker


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
