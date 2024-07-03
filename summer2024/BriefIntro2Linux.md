# A Brief Intro to Linux

Summer Course in Population Genetics

University of Copenhagen


# Introduction

For this course you are not required to know linux beforehand. BUT a lot of the programs we going to use are only available in linux (or easier to use in linux), so below is a brief how to for linux. It should be enough to get by for the entire week. However, if you get stuck in an exercise because there is a linux command you don't understand or can't get to work then PLEASE just ask :)

# Open a terminal and ssh
More or less everything we are going to do in linux will be through a 'terminal'. You will get a login and password to a server.

**Server Name**

You will be working on the server called emily.popgen.dk.

**Windows:**
For Windows users we recommend using MobaXterm.

Download the program [MobaXTerm](http://mobaxterm.mobatek.net/download-home-edition.html), open and press start a new session. Enter server, login and password.

**Mac:**
You should have Xquartz or similar installed so that you can view graphics from the terminal. Else ask Ida when you are here if you have problems.

```
ssh -Y USERNAME@emily.popgen.dk
change USERNAME to your user name
```

**Linux:**
If you use linux then you hopefully know how to open a terminal!

```
ssh -X USERNAME@emily.popgen.dk
change USERNAME to your user name
```

# How to

**Copy-paste**

Copy: To copy paste something mark the text you want to copy by holding your left mouse button down while moving the mouse from one end of the text.

Paste: just point your cursor where you want the text pasted and click once one your middle mouse button.

**Save a file from a browser to the file system**

Right click on the link an press "save-link-as". Then click your way to the directory you want to store the file in and finally click on "save".

**Open a terminal**

More or less everything we are going to do in linux will be through a 'terminal'. To open a terminal click on the button that looks like a small computer screen in the buttom left part of your screen. A terminal is a simple window where you can enter a command and make the computer execute it by pressing enter. Here is an example of typing a command to the commandline:

```
ls
```

**Find your way around in the file system in the terminal**

Here are some very useful commands to help you move around in the file system in the terminal:

ls: shows the content of the directory you are in. Example:

```
ls
subdir1  mytextfile1.txt  mytextfile2.txt  mytextfile3.txt
```

cd: will make you mv to another directory. You can either move down into a subdirectory or move a level up in directories. Examples:

```
cd subdir
cd ..
```

**Change the file system through a terminal**

Here are some commands for copying and moving files and for making new directories:

cp: copies a file to another directory or to another file in the same dirctory. Examples:

```
cp file1.txt subdir1
cp mytextfile1.txt copyofmytextfile1.txt
```

mv: moves a file to another directory or to another file in the same directory (rename a file). Example:
``` 
mv mytextfile1.txt subdir1
mv mytextfile1.txt mytextfile1withanothername.txt
```

mkdir: creates a new subdirectory. Example:

```
mkdir subdir2
```

rm: remove a file, Example:

```
rm mytextfile1.txt
```

**View the content of different types of files from the terminal**

less: shows the content of a text file. In less mode you can use arrows to navigate up and down in the text. To get out of less mode type q. Example:

```
less file1.txt
```

evince: shows the content of pdf, eps, ps and png files. Example:

```
evince myplot.png &
```

**Edit files from the terminal**

gedit: a simple editor that is similar to msome of the simple Windows text editors. Example:

```
gedit mytextfile1.txt &
```

**Make your life easier...**

Use "history": you can go through your previous commands by using the up and down arrows.

Use "tab completion": you don't always have to type the whole filename. Instead try to use the tab button once you have written part of the filename. The computer will then finish the file name if it is the only one filenames that starts with what you wrote or finish the part that it can unambiguously.

Use "regular expressions": if you want to do something with all the files that start or end with something specific you can use a * to denote the vaying part. Example:

```
mv *txt subdir1
```
