The Bash Command Line
=====================

This tutorial makes you familiar with **bash**, the Linux command line.
You will learn to:

-   navigate directories
-   manipulate files
-   execute programs

If you have no previous experience with Unix-like systems or know a few
commands but would like to know more, this tutorial is for you.

------------------------------------------------------------------------

Goal
----

In this tutorial, you will be looking for two words with 22 characters
in total:

All characters are hidden in the exercises below.

This tutorial is run **in a terminal**. Log in to the server, then type the
commands shown in the grey boxes one at a time and read what comes back.


------------------------------------------------------------------------

## 1. Directories and files
------------------------

```bash
mkdir -p IntroToBash
cp /course/data/popgenmsc26_exercises/linux/Exercises.zip IntroToBash/
echo "Practice workspace created at: $(pwd)/IntroToBash"
```

You created your exercises folder `IntroToBash` and moved the Exercises file into it.

Now we need to decompress it so you can find the hidden characters

```bash
cd IntroToBash
unzip Exercises.zip
```

### 1.1. Navigating directories

The **first character** is hidden in a file somewhere in the `exercise_1` directory tree. Use the commands

cd <directory_name>

```bash
cd unix_tutorial
ls
```

This is how you move from one directory to the next. Look through subdirectories until you find one with the name `solution_1.1` and list its contents. If you went to a wrong directory, you can go back one level by typing:

`cd ..`

```bash
cd exercises
ls  
cd exercise_1
ls
```

Look within the directories in `exercise_1`

```bash
cd directoryA
ls  
cd maybe_here
ls
```

It doesn't seem to have the answer. We need to go back and try another directory

`pwd` - will tell you where you are (our path)

```bash
cd ..
pwd
ls  
cd ..
pwd
ls
```

We are back to the first exercise directory. We will look into the second one `directoryB`

```bash
cd directoryB
ls
```

We found `solution_1.1`!! We need to access it and we have our first letter

```bash
cd solution_1.1
ls
```

### 1.2. Show a hidden file

Some files are not visible immediately. To see them, you need the command

`ls -a`

The **second character**, is in the same directory as the first one `solution_1.1`, but in a hidden file.

```bash
ls -a
```

### 1.3. Execute a program

Use cd .. to go back to the directory `exercise_1/directoryB/`. When listing its contents, you should see a **shell script file** `program.sh`.

To find the **third character**, you need to execute the program. On bash, this is done by typing source and the name of the program:

```bash
cd ..
ls  
source program.sh
```

### 1.4. Find out how big a file is

Go to the folder `exercise_1/directoryC/`. To find **the fourth character**, you need to find out how big the text file in the directory is. This is done with the command

`ls -l`

```bash
cd ..
cd directoryC
ls -l
```


In the table the command produces, you will find the file size in bytes,
the file’s owner, permissions to read and modify it, and the date/time
of the last modification.

To obtain the fourth character look up the file size in the [Table of
printable ASCII
characters](https://en.wikipedia.org/wiki/ASCII#Printable_characters)
(use the decimal column):

<img src="https://alpharithms.s3.amazonaws.com/assets/img/ascii-chart/ascii-table-alpharithms-scaled.jpg"
     alt="Printable ASCII character table"
     width="700">

*ASCII Table, Public Domain*

When typing names of directories or files, press `[TAB]` after the first
few characters. Unix tries to guess what you are typing.


------------------------------------------------------------------------

## 2. Edit text files
------------------

For the second section, we will change to the directory `exercise_2`.

```bash
cd ../..
cd exercise_2
ls
```

### 2.1. See what is in a text file

In the directory *exercise\_2/*, you will find a text file *solution\_2.1.txt*. The **fifth character** is inside that file. To see its contents, use the command

`less <filename>`
or you can print the content into the terminal
`cat solution_2.1.txt`

Type ‘q’ on your keyboard to exit the screen.

### 2.2. Edit text files

To get **character number six**, you need to create a text file in the
`exercise_2` directory. Use the editor `nano`:

`nano <filename>`

**To exit nano, type Ctrl-X.**

Along the bottom of the nano screen you will see its list of shortcuts. The
**sixth character** is the one you need to press to save/write a file in `nano`.


------------------------------------------------------------------------

## 3. Copy and remove files
------------------------

Please go to the directory `exercise_3`.

```bash
cd ..
cd exercise_3
ls
```

### 3.1. Create a directory and copy a file to it.

To find **characters seven and eight**, you need to create a subdirectory named *solution* in `exercise_3/` and copy the files from the `code1/` and `code2/` folders into it.

For creating directories, use the command:

`mkdir <directory name>`

For copying, you can use the command

`cp <filename from> <filename to>`

Type `ls -l solution/*` afterwards to see the solution.

```bash
mkdir solution
ls  
cp code1/* solution/
cp code2/* solution/
ls -l solution/*
```

### 3.2. Removing files

In the `data` directory, all files with an `Y` need to be deleted. To do so, use the command:

`rm <filename>`

Also, there are more files to be deleted in the *data* directory. To remove more than one file at once, you can use `*` as a wildcard, i.e. `rm ju*` will delete all of `junk.txt, juniper.txt` and `june.docx`.

To get **characters nine and ten**, look at the files that remain after deleting all that contain a `Y`.

To remove an empty directory, you can use

`rmdir <directory name>`

The command

`rm -r <directory name>`

deletes a directory and everything in it.

**Be careful! On Unix, it is not possible to undelete files!**

This makes removing files with the `*` symbol **very** dangerous, because you could wipe out everything with a single command (e.g. if you type the wrong directory by accident). Backups become an even better idea after learning this command.

```bash
ls data  
rm data/*Y*
ls data
```


------------------------------------------------------------------------

## 4. Process text data
--------------------

Please go to the directory exercise\_4.

```bash
cd ..
cd exercise_4
ls
```

### 4.1. comparing two files

There are two different versions of a quote, `ai.txt`, and `artificial_intelligence.txt`. To find out, how they differ, Unix provides the command

`diff <filename1> <filename2>`

The **11th character** of the solution is the single character in which the two files differ.

```bash
diff ai.txt artificial_intelligence.txt
```

### 4.2. Sorting a text file

Unix has a small program to sort text files alphabetically. It is called by

`sort <filename>`

The **12th character** of the solution is the first character of the last word in the alphabetically sorted file elephant.txt.

To store the sorted lines in a new file, you can add an output file, like

`sort <filename> > result.txt`

```bash
echo "<<this is how elephant.txt looks like>>"
head elephant.txt  
echo "<<this is after sort>>"
sort elephant.txt
```

```bash
sort elephant.txt > result.txt
head result.txt
```

### 4.3. Finding words in a text file

To look for specific words in a text file, use the command

`grep <word> <filename>`

It produces all lines from the given file that contain the given word. The `grep` command is very powerful and can handle Regular Expressions.

To find the **13th character**, search for the word **fire** in the file `datascience.txt` and take the **first** character of the output.

You can search through many files at once by including a \* in the filename.

```bash
grep "fire" datascience.txt
```


--------------------------------------------------------------------------

## 5. Unzip files
--------------

Please go to the directory exercise\_5.

```bash
cd ..
cd exercise_5
ls
```

### 5.1. unzipping archives

Unzipping compressed files is a very basic and important task. On Unix, you often encounter WinZip archives, .tar archives and .gz compressed files. For unpacking Win zip files, use

`unzip <filename>`

for .tar and .tar.gz files

`tar -xf <filename>`

and for .gz files,

`gunzip <filename>`

The **14th and 15th character** of the solution are in a multiply wrapped archive in the exercise\_5 directory.

To pack a directory and everything within, you can use the command

`tar -cf backup.tar <directory>`

To subsequently compress it, use

`gzip backup.tar`

```bash
unzip wrapped.zip
ls  
echo "The file `wrapped.zip` has been unzipped. The contents are now available in the current directory."
tar -xf wrapped.tar
ls
```

You decompressed a directory!

```bash
cd solution
ls  
echo "There is still one file to decompress"

gunzip solution.txt.gz  
echo "Finally!"
head solution.txt
```


------------------------------------------------------------------------

## 6. Command-line tools
---------------------

Please go to the directory `exercise_6`.

```bash
cd ../../exercise_6
ls -l
```

### 6.1. Changing file access rights

Each file on Unix has separate permissions for reading ‘r’, writing ‘w’, and executing ‘x’. Displaying them with:

`ls -l`

There is one triplet of permissions for the owner of the file owner, one triplet for a group of users, and one for all others. The `chmod`
command allows to change these permissions, e.g.

`chmod a+x <filename>`

grants all users the permission to execute a file, while chmod u-w forbids the current user (oneself) to write to the file (thereby protecting it from being deleted accidentally).

To see **characters 16 and 17** of the solution, make the program `permissions.sh` executable. Then execute it with:

`./permissions.sh`

You can grant permissions for a whole directory tree using

`chmod -R a+x <directory>`

```bash
chmod a+x permissions.sh
ls -l
```

To run the script:

```bash
./permissions.sh
```

### 6.2. How much disk space have I left?

To find out, how much disk space you have left, you can use the command

`df`

`df` lists all hard drive partitions, CD-ROMs, pendrives and some logical partitions Unix uses. All numbers are given in kilobyte (1000 byte or one 1000000th GB).

To obtain the **18th character**, check out the version of the `df` program. Find out how to do that with:

`df --version`

The solution is the last character of the first authors’ first name.

```bash
df --version
```

### 6.3. Set an environment variable

To install some programs, it is necessary to set so-called environment variables. These can be set using the command

`export <variable-name>=<value>`

You can see all variables by the command

`env`

To obtain the **19th character**, you need to use `export` to set the variable *GIVEME* to the value **SOLUTION** (i.e. assign the value ‘SOLUTION’ to the variable name ‘GIVEME’).

`echo $GIVEME`

Find out the **character’s position in the alphabet** with:

`echo $GIVEME | wc -c`

The symbol ‘\|’ is called a pipe and is often used to connect Unix programs to each other.

By default, changes to environment variables only affect the current terminal.

If you want to set environment variables for each console window, write the export command to the file `.bashrc` in your home directory (it is a hidden file).

```bash
export GIVEME=SOLUTION
echo $GIVEME
echo $GIVEME | wc -c
```

### 6.4. Check whether you have internet

The easiest way to check from the Unix command line whether the internet connection works, is to send a request to a known server (e.g. www.google.com) using the command

`ping <web address>`

The command reports, how long a message takes back and forth to the given server.

The **20th character** is the `ping` option that sets the maximum number of requests sent. In other words, look for the option that will stop sending packets after sending a user-defined number.

Check the documentation with:

`man ping`

```bash
man ping
```

### 6.5. Managing processes

To see what programs are running on your machine, type

`top`

It displays you a list of all currently active programs. *Shift+P* sorts them by the CPU time they are using, *Shift+M* by the amount of memory they are using (if you don’t see any program consuming lots of memory, start a web browser). Quit `top` by pressing *q*.

The **last two characters** of the solution are the first two characters of the second word in the line containing the column labels.

If you want to get rid of one of the programs you started (e.g. because it crashed), you can do so by typing

`kill -s 9 <pid>`

You find the pid number in the first column of the *top* output. Of course, you may only interrupt your own programs, not those owned by *root*, the system administrator.

**Congratulations! You’ve completed the tutorial :)**


------------------------------------------------------------------------

### License

**© 2010 Dr. Kristian Rother**

This tutorial is published under the Creative Commons Attribution
Share-alike License 4.0

You can find the full sources on
[\<https://github.com/krother/bash\_tutorial\>](https://github.com/krother/bash_tutorial).

### Acknowledgements

I thank Janusz M. Bujnicki, Allegra Via, Pedro Fernandes and Joachim
Jacob for their help with testing and reviewing the material. Further
thanks go to the German Academic Exchange Service (DAAD) for financial
support.

### Contact

Dr. Kristian Rother

[krother@academis.eu](krother@academis.eu)

[www.academis.eu](www.academis.eu)
