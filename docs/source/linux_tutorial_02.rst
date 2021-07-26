 .. _linux-tutorial-02-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

============
Tutorial two
============

Redirection
===========

Most processes initiated by UNIX commands write to the standard output (that is, they write to the terminal screen), and many take their input from the standard input (that is, they read it from the keyboard). There is also the standard error, where processes write their error messages, by default, to the terminal screen.

We have already seen one use of the ``cat`` command to write the contents of a file to the screen.

Now type ``cat`` without specifing a file to read

.. code-block:: bash

   $ cat

Then type a few words on the keyboard and press the :boldcode:`[Return]` key.

Finally hold the :boldcode:`[Ctrl]` key down and press :boldcode:`[d]` (written as ``^D`` for short) to end the input.

What has happened?

If you run the ``cat`` command without specifing a file to read, it reads the standard input (the keyboard), and on receiving the'end of file' (^D), copies it to the standard output (the screen).

In UNIX, we can redirect both the input and the output of commands.

Redirecting the Output
======================

We use the ``>`` symbol to redirect the output of a command. For example, to create a file called **list1** containing a list of fruit, type  

.. code-block:: bash

   $ cat > list1

Then type in the names of some fruit. Press [Return] after each one.

.. code-block::

   pear
   banana
   apple
   ^D (Control D to stop)

What happens is the ``cat`` command reads the standard input (the keyboard) and the ``>`` redirects the output, which normally goes to the screen, into a file called **list1**.

To read the contents of the file, type

.. code-block:: bash

   $ cat list1

.. important::

   Exercise 2a: Using the above method, create another file called list2 containing the following fruit: orange, plum, mango, grapefruit. Read the contents of list2

   The form >> appends standard output to a file. So to add more items to the file list1, type

   % cat >> list1

   Then type in the names of more fruit

   peach
   grape
   orange
   ^D (Control D to stop)

   To read the contents of the file, type

   % cat list1

   You should now have two files. One contains six fruit, the other contains four fruit. We will now use the cat command to join (concatenate) list1 and list2 into a new file called biglist. Type

   % cat list1 list2 > biglist

   What this is doing is reading the contents of list1 and list2 in turn, then outputing the text to the file biglist

   To read the contents of the new file, type

   % cat biglist