#!/usr/bin/python3

import os

if not os.path.isfile('unixstuff/science.txt'):
    print("The file `unixstuff/science.txt` is not present. Make sure to copy `/home/jovyan/linux_tutorial/science.txt` to `/home/jovyan/linux_tutorial/tutorial_01/unixstuff/science.txt`.")
else:
    print("Success! The file `unixstuff/science.txt` is where it shoudl be!")

if not os.path.isfile('unixstuff/backups/science.bak'):
    print("The file `unixstuff/backups/science.bak` is not present. Make sure to copy `/home/jovyan/linux_tutorial/tutorial_01/unixstuff/science.txt` to `/home/jovyan/linux_tutorial/tutorial_01/unixstuff/backups/science.bak`.")
else:
    print("Success! The file `unixstuff/backups/science.bak` is where it shoudl be!")

if not os.path.getsize('unixstuff/backups/science.bak') > 0:
    print("The file `unixstuff/backups/science.bak` is empty. Try again to copy it.")

if os.path.isfile('unixstuff/tempfile.txt'):
    print("You created the file `unixstuff/tempfile.txt`, but you didn't delete it. Try again.")
else:
    print("Success! The file `unixstuff/tempfile.txt` was deleted.")

if os.path.isdir('unixstuff/tempstuff'):
    print("There is still a directory called `unixstuff/tempstuff`. Try to remove it again.")
else:
    print("Success! The directory `unixstuff/tempstuff` does not exist.")
