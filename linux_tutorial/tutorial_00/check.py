#!/usr/bin/python3

import os

if os.path.isdir('unixstuff'):
    print("Success! `unixstuff` directory is present!")
else:
    print("The `unixstuff` directory is not present. Try again and make sure, that you are in `/home/jovyan/linux_tutorial/tutorial_00`, when you make the directory.")

if os.path.isdir('unixstuff/backups'):
    print("Success! `unixstuff/backups` directory is present!")
else:
    print("The `unixstuff/backups` directory is not present. Try again and make sure, that you are in `/home/jovyan/linux_tutorial/tutorial_00`, when you make the directory.")
