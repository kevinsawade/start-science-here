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

boldcode:`ls (list)`

When you first login, your current working directory is your home directory. Your home directory has the same name as your user-name (for us, that username is ``jovyan`` and the home directory is ``/home/jovyan``). Normally your personal files can be found here, but for us, the content of the whole "start science here" project is here. To find out, what exactly can be found here, tye:

.. code-block:: bash

   $ ls

Sample text.

Heading 1
=========

Sample text2

Subheading 1
------------

Sample Text again

Heading 2
=========

More text

Heading 3
=========

Even more