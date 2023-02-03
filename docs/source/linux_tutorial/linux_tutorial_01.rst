 .. _linux-tutorial-01-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

============
Tutorial one
============

Copying files
=============

:boldcode:`cp (copy)`
---------------------

`cp` :italiccode:`file1` :italiccode:`file2` is the command which makes a copy of **file1** in the current working directory and calls it **file2**.

What we are going to do now, is to take a file stored in ``/home/jovyan/linux_tutorial`` and use the ``cp`` comamnd to copy it to ``/home/jovyan/linux_tutorial/tutorial_01``.

First, cd to your unixstuff directory.

.. code-block:: bash

   $ cd ~/linux_tutorial/tutorial_01/unixstuff


Then at the UNIX prompt, type,

.. code-block:: bash

   $ cp /home/jovyan/linux_tutorial/science.txt .

(Note: Don't forget the dot (``.``) at the end. Remember, in UNIX, the dot means the current directory.)

The above command means copy the file **science.txt** to the current directory, keeping the name the same.

.. sshexercise:: Exercise 1a

   Create a backup of your **science.txt** file by copying it to a file called **science.bak**.

Moving files
============

:boldcode:`mv (move)`
---------------------

`mv` :italiccode:`file1` :italiccode:`file2` moves (or renames) **file1** to **file2**.

To move a file from one place to another, use the ``mv`` command. This has the effect of moving rather than copying the file, so you end up with only one file rather than two.

It can also be used to rename a file, by moving the file to the same directory, but giving it a different name.

We are now going to move the file science.bak to your backup directory.

First, change directories to your unixstuff directory (can you remember how?). Then, inside the unixstuff directory, type

.. code-block:: bash

   $ mv science.bak backups/.


Type ``ls`` and ``ls backups`` to see if it has worked.

Removing files and directories
==============================

:boldcode:`rm (remove), rmdir (remove directory)`
-------------------------------------------------

To delete (remove) a file, use the ``rm`` command. As an example, we are going to create a copy of the **science.txt** file then delete it.

Inside your unixstuff directory, type

.. code-block:: bash

   $ cp science.txt tempfile.txt
   $ ls # (to check if it has created the file)
   $ rm tempfile.txt
   $ ls # (to check if it has deleted the file)

You can use the ``rmdir`` command to remove a directory (make sure it is empty first). Try to remove the backups directory. You will not be able to since UNIX will not let you remove a non-empty directory.

.. sshexercise::

   Exercise 1b: Create a directory called ``tempstuff`` using ``mkdir``, then remove it using the ``rmdir`` command.


Displaying the contents of a file on the screen
===============================================

:boldcode:`clear (clear screen)`
--------------------------------

Before you start the next section, you may like to clear the terminal window of the previous commands so the output of the following commands can be clearly understood.

At the prompt, type

.. code-block:: bash

   $ clear

This will clear all text and leave you with the $ prompt at the top of the window.

:boldcode:`cat (concatenate)`
-----------------------------

The command ``cat`` can be used to display the contents of a file on the screen. Type:

.. code-block:: bash

   $ cat science.txt

As you can see, the file is longer than than the size of the window, so it scrolls past making it unreadable.

I am obliged to also tell you that ``cat``'s original intention is to combine multiple files into one (that's what concatenation is). So normally one would call:

.. code-block:: bash

   $ cat file1.txt file2.txt >> out_file.txt

But ``cat`` is easy to enter and works just fine, if you like to get a glimpse at a file.

:boldcode:`less`
----------------

The command ``less`` writes the contents of a file onto the screen a page at a time. Type

.. code-block:: bash

   $ less science.txt

Press the :boldcode:`[space-bar]` if you want to see another page, type :boldcode:`[q]` if you want to quit reading. As you can see, less is used in preference to cat for long files.

:boldcode:`head`
----------------

The ``head`` command writes the first ten lines of a file to the screen.

First clear the screen then type

.. code-block:: bash

   $ head science.txt

Then type

.. code-block:: bash

   $ head -5 science.txt

What difference did the ``-5`` do to the head command?

:boldcode:`tail`
----------------

The ``tail`` command writes the last ten lines of a file to the screen.

Clear the screen and type

.. code-block:: bash

   $ tail science.txt

How can you view the last 15 lines of the file?

.. code-block:: bash

   $ tail -15 science.txt

Searching the contents of a file
================================

Simple searching using less
---------------------------

Using ``less``, you can search though a text file for a keyword (pattern). For example, to search through science.txt for the word 'science', type

.. code-block:: bash

   $ less science.txt

then, still in less (i.e. don't press :boldcode:`[q]` to quit), type a forward slash :boldcode:`[/]:boldcode:` followed by the word to search


.. code-block::

   /science

As you can see, less finds and highlights the keyword. Type :boldcode:`[n]` to search for the next occurrence of the word.

grep (don't ask why it is called grep)
--------------------------------------

``grep`` is one of many standard UNIX utilities. It searches files for specified words or patterns. First clear the screen, then type

.. code-block:: bash

   $ grep science science.txt

As you can see, grep has printed out each line containg the word science.

Or has it????

Try typing

.. code-block:: bash

   $ grep Science science.txt

The grep command is case sensitive; it distinguishes between Science and science.

To ignore upper/lower case distinctions, use the ``-i`` option, i.e. type

.. code-block:: bash

   $ grep -i science science.txt

To search for a phrase or pattern, you must enclose it in single quotes (the apostrophe symbol). For example to search for spinning top, type

.. code-block:: bash

   $ grep -i 'spinning top' science.txt

Some of the other options of grep are:

``-v`` display those lines that do NOT match
``-n`` precede each maching line with the line number
``-c`` print only the total count of matched lines
Try some of them and see the different results. Don't forget, you can use more than one option at a time, for example, the number of lines without the words science or Science is

.. code-block:: bash

   $ grep -ivc science science.txt

There is also the option to print the following or preceding lines of a line matching the grep pattern. That can come in handy, when you search for some code or other text you wrote long ago and want to have the context of the matching line and not just the line itself. ``-A`` will print lines **A**fter the match, ``-B`` will print lines **B**efore the match.

.. code-block:: bash

   $ grep -A 5 # 5 lines after match
   $ grep -A 5 -B 5 # 5 lines before and after match

:boldcode:`wc (word count)`
---------------------------

A handy little utility is the ``wc`` command, short for word count. To do a word count on science.txt, type

.. code-block:: bash

   $ wc -w science.txt

To find out how many lines the file has, type

.. code-block:: bash

   $ wc -l science.txt


Calling the check.py scripts
============================

To check whether your exercises succeeded type:

.. code-block:: bash

   $ cd /home/jovyan/linux_tutorial/tutorial_01
   $ python3 check.py


Summary
=======

+--------------------------+-------------------------------------------------+
| cp file1 file2           | copy file1 and call it file2                    |
+==========================+=================================================+
| ``mv file1 file2``       | move or rename file1 to file2                   |
+--------------------------+-------------------------------------------------+
| ``rm file``              | remove a file                                   |
+--------------------------+-------------------------------------------------+
| ``rmdir directory``      | remove a directory                              |
+--------------------------+-------------------------------------------------+
| ``cat file``             | display a file                                  |
+--------------------------+-------------------------------------------------+
| ``more file``            | display a file a page at a time                 |
+--------------------------+-------------------------------------------------+
| ``head file``            | display the first few lines of a file           |
+--------------------------+-------------------------------------------------+
| ``tail file``            | display the last few lines of a file            |
+--------------------------+-------------------------------------------------+
| ``grep 'keyword' file``  | search a file for keywords                      |
+--------------------------+-------------------------------------------------+
| ``wc file``              | count number of lines/words/characters in file  |
+--------------------------+-------------------------------------------------+


Continue
========

Continue to the next exercise: :ref:`linux-tutorial-02-label`
