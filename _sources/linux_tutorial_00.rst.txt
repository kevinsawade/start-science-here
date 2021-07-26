 .. _linux-tutorial-00-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

=============
Tutorial zero
=============

After accessing the terminal on binder you will be greeted by a black screen showing the prompt:

.. image:: _static/pics/linux_tutorial/example_terminal_binder.png
   :target: _static/pics/linux_tutorial/example_terminal_binder.png
   :alt: Terminal Session on Binder showing blank.

We will now make our way through the first commands:

Listing files and directories
=============================

:boldcode:`ls (list)`

When you first login, your current working directory is your home directory. Your home directory has the same name as your user-name (for us, that username is ``jovyan`` and the home directory is ``/home/jovyan``). Normally your personal files can be found here, but for us, the content of the whole "start science here" project is here. To find out, what exactly can be found here, tye:

.. code-block:: bash

   $ ls
   binder  docs  ideas.txt  LICENSE  pics  python_tutorial  README.md  requirements.txt

The ``ls`` command lists the contents of your current working directory.

ls does not, in fact, cause all the files in your home directory to be listed, but only those ones whose name does not begin with a dot (.) Files beginning with a dot (.) are known as hidden files and usually contain important program configuration information. They are hidden because you should not change them unless you are very familiar with UNIX!

To list all files in your home directory including those whose names begin with a dot, type:

.. code-block:: bash

   $ ls -a

Now you can see the hidden ``.git`` directory and the ``.gitignore`` file. This command also lists two files, ``.`` and ``..``. The first is the file of the current working cirectory and the second is the file of the parent directory (in our case that would be ``/home``).

The ``-a`` part of the ``ls`` command is reffered to as an option (sometimes also called flag). Options change the behaviour of the command. There are online manual pages that tell you which options a particular command can take, and how each option modifies the behaviour of the command. (See later in this tutorial)

Changing into a different directory
===================================

:boldcode:`cd (change directory)`

The command ``cd`` :italiccode:`directory` means change the current working directory to 'directory'. The current working directory may be thought of as the directory you are in, i.e. your current position in the file-system tree.

To change to the directory you have just made, type:

.. code-block:: bash

   $ cd linux_tutorial

You can have a look around that directory with the ``ls`` or ``ls -a`` command. Now the two files ``.`` and ``..`` come into play. Typing

.. code-block:: bash

   $ cd .

Means: Stay in the same directory (nothing changes), but typing:

.. code-block:: bash

   $ cd ..

Moves you back up one directory, back to ``/home/jovjan`` (enter ``pwd`` to check your current working directory).

Now we will move two directories at once. So type:

.. code-block:: bash

   $ cd linux_tutorial/tutorial_00

And have a look around with ``ls``. To swith between your two most recent directories, you can use ``cd -``. Type these commands in succession:

.. code-block:: bash

   $ pwd
   /home/jovyan/linux_tutorial/tutorial_00
   $ cd -
   $ pwd
   /home/jovyan
   $ cd -
   $ pwd
   /home/jovyan

The ``cd`` command without anything will bring you back to your home directoy.


Making directories
==================

:boldcode:`mkdir (make directory)`

We will now make a subdirectory in the ``~/linux_tutorial/tutorial_00`` directory to hold the files you will be creating and using in the course of this tutorial.

.. note::

   All changes and files that you upload to binder will be deleted. Please refer to :ref:`downloading-files-label` for how to download files that you might want to keep.


To make a subdirectory called ``unixstuff`` in your current working directory type:

.. code-block:: bash

   $ mkdir unixstuff

Verify that your directory creation was successful by calling the ``ls`` command.

.. important::

   Your exercise: Make another directory inside the ``unixstuff`` directory called ``backups``.


Pathnames
=========

:boldcode:`pwd (print working directory)`

We've already used the ``pwd`` command extensively. But let's talk about the *filesystem*. The filesystem controls how data is stored. We can traverse it with the above commands. The filesystem imposes limits on our PC. The disk can be full or the filesystem can have a largest possible filesize. In Linux, the root filesystem is denoted as ``/`` and you can change to it with:

.. code-block:: bash

   $ cd /

.. note::

   You can always get back to your home directory by calling ``cd``.

Here you will find somewhat imposing directories. You can read up about them here: https://www.howtogeek.com/117435/htg-explains-the-linux-directory-structure-explained/

You can also skip this part and return to ``~/linux_tutorial/tutorial_00``

Understanding pathnames
=======================

Type:

.. code-block:: bash

   $ ls unixstuff

Now type

.. code-block:: bash

   $ ls backups

You will get a message like this -

.. code-block:: bash

   backups: No such file or directory

The reason is, ``backups`` is not in your current working directory. To use a command on a file (or directory) not in the current working directory (the directory you are currently in), you must either cd to the correct directory, or specify its full pathname. To list the contents of your backups directory, you must type

.. code-block:: bash

   $ ls unixstuff/backups

~ (your home directory)
=======================

Home directories can also be referred to by the tilde ``~`` character. It can be used to specify paths starting at your home directory. These tow commands result in you landing in the same directory:

.. code-block:: bash

   $ cd /home/jovyan/linux_tutorial/tutorial_00
   $ cd ~/linux_tutorial/tutorial_00

The command:


.. code-block:: bash

   $ ls ~/unixstuff

will list the contents of your unixstuff directory, no matter where you currently are in the file system.


.. important::

   Exercises:

What do you think

.. code-block:: bash

   $ ls ~

would list?

What do you think

.. code-block:: bash

   $ ls ~/..

would list?

Calling the check.py scripts
============================

In every ``linux_tutorial/tutorial_XX`` directory (where XX can be 00, 01, 02, etc.) there is a file called ``check.py`` which will check, whether your Exercises have been executed correctly. You can execute them with this command:

.. code-block:: bash

   $ ./check.py

So type these commands in succession:

.. code-block:: bash

   $ cd /home/jovyan/linux_tutorial/tutorial_00
   $ ./check.py

And you will either see:

.. code-block::

   Success! `unixstuff` directory is present!

Or:

.. code-block::

   The `unixstuff` directory is not present. Try again and make sure, that you are in `/home/jovyan/linux_tutorial/tutorial_00`, when you make the directory.

Continue
========

Continue to the next exercise: 