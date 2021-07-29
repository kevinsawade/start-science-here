 .. _linux-tutorial-04-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

=============
Tutorial four
=============

File system security (access rights)
------------------------------------

In your ``unixstuff`` directory, type

.. code-block:: bash

   $ ls -l # (l for long listing!)

You will see that you now get lots of details about the contents of your directory, similar to this:

.. code-block::

   -rw-r--r-- 1 jovyan users 700 Jul 27  2021 file1
   drwxr-xr-x 1 jovyan users  78 Mar 21  2021 directory

**1st column:**

These associated access rights (each file and directory has them). The first symbol (``-``, or ``d``) means file or directory. The next 9 symbols are 3 blocks of three symbols. They are:

* ``r``: read (file), list (directory)
* ``w``: write (file), change files in dir (directory)
* ``x``: execute (file), access files (directory), given that these files have read permission for you.

They are always in this order and a ``-`` means, that this permission is not given.

The first block gives the rights of the owner (``jovyan``), the second block of the group (``ls -lg`` adds the group to the ``ls -l`` command and ``groups`` prints your current groups), the last 3 characters give the permissions for everyone.

+---------+-----------------------------------+
| Column  | Content                           |
+=========+===================================+
| 2       | The number of links to that file  |
+---------+-----------------------------------+
| 3       | The owner of the file             |
+---------+-----------------------------------+
| 4       | The group that owns the file      |
+---------+-----------------------------------+
| 5       | Size of file                      |
+---------+-----------------------------------+
| 6       | Date of creation                  |
+---------+-----------------------------------+
| 7       | Name of file                      |
+---------+-----------------------------------+


**Some examples:**

+-------------+------------------------------------------------------------------------------------------------------------------------------------------+
| -rwxrwxrwx  | a file that everyone can read, write and execute (and delete).                                                                           |
+=============+==========================================================================================================================================+
| -rw-------  | a file that only the owner can read and write - no-one else can read or write and no-one has execution rights (e.g. your mailbox file).  |
+-------------+------------------------------------------------------------------------------------------------------------------------------------------+

Changing access rights
======================

:boldcode:`chmod (change mode)`

Only the owner of a file can use ``chmod`` to change the permissions of a file. The options of chmod are as follows

+---------+-----------------------+
| Symbol  | Meaning               |
+=========+=======================+
| u       | user                  |
+---------+-----------------------+
| g       | group                 |
+---------+-----------------------+
| o       | other                 |
+---------+-----------------------+
| a       | all                   |
+---------+-----------------------+
| r       | read                  |
+---------+-----------------------+
| w       | write                 |
+---------+-----------------------+
| x       | execute               |
+---------+-----------------------+
| \+      | add permission        |
+---------+-----------------------+
| \-      | take away permission  |
+---------+-----------------------+

For example, to remove read write and execute permissions on the file biglist for the group and others, type

.. code-block:: bash

   $ chmod go-rwx biglist

This will leave the other permissions unaffected.

To give read and write permissions on the file biglist to all,

.. code-block:: bash

   $ chmod a+rw biglist

.. exercise:: Exercise 4a

   Try changing access permissions on the file science.txt and on the directory backups

   Use ls -l to check that the permissions have changed.


Processes and Jobs
==================

A process is an executing program identified by a unique PID (process identifier). To see information about your processes, with their associated PID and status, type

.. code-block:: bash

   $ ps

A process may be in the foreground, in the background, or be suspended. In general the shell does not return the UNIX prompt until the current process has finished executing.

Some processes take a long time to run and hold up the terminal. Backgrounding a long process has the effect that the UNIX prompt is returned immediately, and other tasks can be carried out while the original process continues executing.

Running background processes
----------------------------

To background a process, type an ``&`` at the end of the command line. For example, the command ``sleep`` waits a given number of seconds before continuing. Type

.. code-block:: bash

   $ sleep 10

This will wait 10 seconds before returning the command prompt ``$``. Until the command prompt is returned, you can do nothing except wait.

To run sleep in the background, type

.. code-block:: bash

   $ sleep 10 &
   6259

The ``&`` runs the job in the background and returns the prompt straight away, allowing you do run other programs while waiting for that one to finish.

The first line in the above example is typed in by the user; the next line, indicating job number and PID, is returned by the machine. The user is be notified of a job number (numbered from 1) enclosed in square brackets, together with a PID and is notified when a background process is finished. Backgrounding is useful for jobs which will take a long time to complete.

Backgrounding a current foreground process
------------------------------------------

At the prompt, type

.. code-block:: bash

   $ sleep 100

You can suspend the process running in the foreground by holding down the :boldcode:`[Control]` key and typing :boldcode:`[z]` (written as ``^Z``) Then to put it in the background, type

.. code-block:: bash

   $ bg

.. note::

   do not background programs that require user interaction e.g. pine

Listing suspended and background processes
==========================================

When a process is running, backgrounded or suspended, it will be entered onto a list along with a job number. To examine this list, type

.. code-block:: bash

   $ jobs

An example of a job list could be

.. code-block:: bash

   [1] Suspended sleep 100
   [2] Running netscape
   [3] Running nedit

To restart (foreground) a suspended processes, type

.. code-block:: bash

   $ fg %jobnumber

For example, to restart sleep 100, type

.. code-block:: bash

   $ fg %1

Typing ``fg`` with no job number foregrounds the last suspended process.

Killing a process
=================

:boldcode:`kill (terminate or signal a process)`

It is sometimes necessary to kill a process (for example, when an executing program is in an infinite loop)

To kill a job running in the foreground, type ``^C`` (:boldcode:`[Control]` +  :boldcode:`c`). For example, run

.. code-block:: bash

   $ sleep 100
   $ ^C

To kill a suspended or background process, type

.. code-block:: bash

   $ kill %jobnumber

For example, run

.. code-block:: bash

   $ sleep 100 &
   $ jobs

If it is job number 4, type

.. code-block:: bash

   $ kill %4

To check whether this has worked, examine the job list again to see if the process has been removed.

:boldcode:`ps (process status)`

Alternatively, processes can be killed by finding their process numbers (PIDs) and using kill PID_number

.. code-block:: bash

   $ sleep 100 &
   $ ps

   PID TT S TIME COMMAND
   20077 pts/5 S 0:05 sleep 100
   21563 pts/5 T 0:00 netscape
   21873 pts/5 S 0:25 nedit

To kill off the process sleep 100, type

.. code-block:: bash

   $ kill 20077

and then type ``ps`` again to see if it has been removed from the list.

If a process refuses to be killed, uses the ``-9`` option, i.e. type

.. code-block:: bash

   $ kill -9 20077

.. note::

   It is not possible to kill off other users' processes.

Summary
=======

+---------------------------+--------------------------------------------+
| ``ls -lag``               | list access rights for all files           |
+===========================+============================================+
| ``chmod [options] file``  | change access rights for named file        |
+---------------------------+--------------------------------------------+
| ``command &``             | run command in background                  |
+---------------------------+--------------------------------------------+
| ``^C``                    | kill the job running in the foreground     |
+---------------------------+--------------------------------------------+
| ``^Z``                    | suspend the job running in the foreground  |
+---------------------------+--------------------------------------------+
| ``bg``                    | background the suspended job               |
+---------------------------+--------------------------------------------+
| ``jobs``                  | list current jobs                          |
+---------------------------+--------------------------------------------+
| ``fg %1``                 | foreground job number 1                    |
+---------------------------+--------------------------------------------+
| ``kill %1``               | kill job number 1                          |
+---------------------------+--------------------------------------------+
| ``ps``                    | list current processes                     |
+---------------------------+--------------------------------------------+
| ``kill 26152``            | kill process number 26152                  |
+---------------------------+--------------------------------------------+


Continue
========

Continue to the next exercise: :ref:`linux-tutorial-05-label`