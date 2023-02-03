 .. _linux-tutorial-05-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

=============
Tutorial five
=============

Other useful UNIX commands
--------------------------

:boldcode:`df`
The df command reports on the space left on the file system. For example, to find out how much space is left on the fileserver, type

.. code-block:: bash

   $ df .

:boldcode:`du`
The du command outputs the number of kilobyes used by each subdirectory. Useful if you have gone over quota and you want to find out which directory has the most files. In your home-directory, type

.. code-block:: bash

   $ du

:boldcode:`compress`
This reduces the size of a file, thus freeing valuable disk space. For example, type

.. code-block:: bash

   $ ls -l science.txt

and note the size of the file. Then to compress science.txt, type

.. code-block:: bash

   $ compress science.txt

This will compress the file and place it in a file called science.txt.Z

To see the change in size, type ls -l again.

To uncomress the file, use the uncompress command.

.. code-block:: bash

   $ uncompress science.txt.Z

:boldcode:`gzip`
This also compresses a file, and is more efficient than compress. For example, to zip science.txt, type

.. code-block:: bash

   $ gzip science.txt

This will zip the file and place it in a file called science.txt.gz

To unzip the file, use the gunzip command.

.. code-block:: bash

   $ gunzip science.txt.gz

:boldcode:`file`
file classifies the named files according to the type of data they contain, for example ascii (text), pictures, compressed data, etc.. To report on all files in your home directory, type

.. code-block:: bash

   $ file *

:boldcode:`history`
The C shell keeps an ordered list of all the commands that you have entered. Each command is given a number according to the order it was entered.

.. code-block:: bash

   $ history (show command history list)

If you are using the C shell, you can use the exclamation character (!) to recall commands easily.

.. code-block:: bash

   $ !! # (recall last command)
   $ !-3 # (recall third most recent command)
   $ !5 # (recall 5th command in list)
   $ !grep # (recall last command starting with grep)

You can increase the size of the history buffer by typing

.. code-block:: bash

   $ set history=100

Continue
========

Continue to the next exercise: :ref:`linux-tutorial-06-label`
