 .. _linux-tutorial-06-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

============
Tutorial six
============

UNIX Variables
==============

Variables are a way of passing information from the shell to programs when you run them. Programs look "in the environment" for particular variables and if they are found will use the values stored. Some are set by the system, others by you, yet others by the shell, or any program that loads another program.

Standard UNIX variables are split into two categories, environment variables and shell variables. In broad terms, shell variables apply only to the current instance of the shell and are used to set short-term working conditions; environment variables have a farther reaching significance, and those set at login are valid for the duration of the session. By convention, environment variables have UPPER CASE and shell variables have lower case names.

Environment Variables
=====================

An example of an environment variable is the OSTYPE variable. The value of this is the current operating system you are using. Type

.. code-block:: bash

   $ echo $OSTYPE

More examples of environment variables are

* USER (your login name)
* HOME (the path name of your home directory)
* HOST (the name of the computer you are using)
* ARCH (the architecture of the computers processor)
* DISPLAY (the name of the computer screen to display X windows)
* PRINTER (the default printer to send print jobs)
* PATH (the directories the shell should search to find a command)

Finding out the current values of these variables
-------------------------------------------------

ENVIRONMENT variables are set using the ``setenv`` command, displayed using the ``printenv`` or ``env`` commands, and unset using the ``unsetenv`` command.

To show all values of these variables, type

.. code-block:: bash

   $ printenv | less


Shell variables
===============

An example of a shell variable is the history variable. The value of this is how many shell commands to save, allow the user to scroll back through all the commands they have previously entered. Type

.. code-block:: bash

   $ echo $history

More examples of shell variables are

* cwd (your current working directory)
* home (the path name of your home directory)
* path (the directories the shell should search to find a command)
* prompt (the text string used to prompt for interactive commands shell your login shell)

Finding out the current values of these variables
-------------------------------------------------

SHELL variables are both ``set`` and displayed using the set command. They can be unset by using the ``unset`` command.

To show all values of these variables, type

.. code-block:: bash

   $ set | less

So what is the difference between PATH and path ?
-------------------------------------------------

In general, environment and shell variables that have the same name (apart from the case) are distinct and independent, except for possibly having the same initial values. There are, however, exceptions.

Each time the shell variables home, user and term are changed, the corresponding environment variables HOME, USER and TERM receive the same values. However, altering the environment variables has no effect on the corresponding shell variables.

PATH and path specify directories to search for commands and programs. Both variables always represent the same directory list, and altering either automatically causes the other to be changed.

Using and setting variables
===========================

Each time you login to a UNIX host, the system looks in your home directory for initialisation files. Information in these files is used to set up your working environment. The bash shell (``$``)  uses a file called ``.bashrc`` (note that this file's name begins with a dot).

At login the bash shell reads the ``bashrc``. This file is used to prepare the user environment and set environment variables that other programs might access.

Finished
========

You've reached the end. Take a break and come back, when you feel like it.