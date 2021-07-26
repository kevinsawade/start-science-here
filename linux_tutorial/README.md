# Linux tutorial

Linux is an open-source operating system, first released in 1991. Linux is what's running on virtually all supercomputers and cluster all over the world. Linux is very modular and comes often as the underlying software of so-called distributions. Such distributions are the popular Ubuntu/Debian Family, of the commercial Red Hat Enterprise Linux distribution. The most significant difference between Windows and Linux for new users is the usage text-only commands. Instead of opening applications or settings pages and working your way through it with a mouse, you just open a terminal and type away.


## Working with the terminal

First, let's introduce you to the terminal. The terminal allows you to control your computer, make changes to files, move files, and start calculations using hundreds of processors on a remote supercomputer. Let's open a terminal by pressing ```Ctrl``` + ```Alt``` + ```T``` simultaneously. A new window will open and you are greeted by the *prompt*. The *prompt* looks like this:

```bash
username@computername:~$
```
The *prompt* contains useful information, such as your username, the name of the computer and your current directory. The tilde symbol (~) shows you, that you are in your current come directory. Try and change the directory to your Downloads directory with the ```cd``` command. Commands which are to be entered into the console are often prepended with a dollar sign ($):

```bash
$ cd Downloads
```

The *prompt* will change and display the following:

```bash
username@computername:~/Downloads$
```

To change back to your home directory you can "climb up" the directory tree once with the ```cd``` command. Enter:

```bash
$ cd ..
```

You can do more awesome stuff in the terminal. If you want to learn more you can visit the linux tutorial **AFTER** completing the step "working with jupyter notebooks"

## Check if git is installed

The files are hosted on github. Github is a service by Microsoft which allows to collaborate on through the internet on software development. We use the program *git* to track changes, control versions of programs, and upload these changes to github.com. To download projects from github.com the program *git* needs to be installed on your computer. Let us check if the program is installed. Open a terminal and enter the following command:

```bash
$ whereis git
```

If this returns something like this:

```bash
git: /usr/bin/git /usr/share/man/man1/git.1.gz
```

You are all set and can continue. If this command returns something like this:

```bash
git:
```

*git* is not installed on your computer. Please ask someone to help you installing *git*

## Download the files

To avoid clutter it is advised to download this tutorial into a directory that is not the home directory (~). Let's create a new directory four all our *git* downloads.

```bash
$ mkdir git
$ cd git
```

In this directory we can donwload the tutorial files with:

```bash
$ git clone THIS NEEDS TO BE CHANGED WHEN THE TUTORIAL HAS BEEN PUSHED TO GITHUB.COM
```

## Working with jupyter notebooks

The next part of this tutorial can be found in a html file. To open it enter:

```bash
$ cd ~/git/tutorial_new_members
$ firefox 00_introduction_notebooks.html &
```