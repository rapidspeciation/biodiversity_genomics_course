
# Introduction to command line

## Log in to server

```bash
ssh genomics@toko.uncu.edu.ar

```

Where are you?
Use print working directory (pwd)
```bash
pwd
```

Which file are there?
```bash
ls ./

```

## Make your working directory
You have seen a list representing the directories and files in this directory.

We can change to another directory by using the command 'change directory' cd
Go into the users directory

```bash
cd users

```
Here is where you make your directory for the course. 
This is your own directory 
Make a directory with your first name and the first letter of your surname, common letters connected with an underscore _

mkdir karin_n

Remember no spaces! No special characters!

```bash
# first make sure you are in the right place 
pwd

# you should be in /home/genomics/scratch/

mkdir firstname_x # change this to your name!!
```

Does everyone have a directory with the correct name?


Go to that directory
```bash
cd karin_n

```

This is your working place for the next week!

## Navigating in the filesystem

Create a file

```bash
touch test_file.txt

```

```
ls ./
```

```bash
touch test_file2.txt

```


ls have many options, check the manual of ls

```bash
man ls

```
Exit the manual by typing q

```bash

q

```

Check out some of the list options
Make a long list (more information on each file), human readable, time sorted

```bash
ls -l ./
ls -h ./
ls -t ./

```

And combine the options
```bash
ls -lht ./

```

Create a directory

```bash
mkdir test_dir1 test_dir2

ls ./

```

Change diretory

```bash
cd test_dir1

ls ./
```


```bash
touch file3
touch file4

ls ./
```

Go back to your first directory
```bash
cd ../

pwd
ls ./
```


Relative paths

```bash

cd test_dir2

pwd

ls ./

ls /home/genomics/scratch/yourname_x/test_dir2

ls ../

ls ../test_dir1/


```

Use the relative path to go to the other directory

```bash
cd ../test_dir1/
pwd

```


# Viewing files

Printing the content to screen with the command cat


```bash

touch new_file.txt
cat new_file.txt 

```

Well it is empty, so we need some content.
nano is a text editor, but there are others vi, vim, emacs, etc

Open the file with nano and write something.
Save and exit with Ctrl + X

```bash

nano new_file.txt
cat new_file.txt 

```
If we do not want the whole file printed to screen

```bash
less new_file.txt

```
You can search in less with /text
n next match
N previous match

Quit with q


Take a look at the first lines of the file
```bash

head file
```

Take a look at the first lines of the file
```bash

tail file
```


# Moving and copying

Command (cp or mv) source target

Source is the file (or directory) you want to copy or move
Target is the file (or directory) you want to copy or move it to

```bash

mv new_file.txt ../

```
Mv can also be used to rename files

```bash
touch new_file2.txt
mv new_file2.txt new_name.txt

```


Copy
```bash

cp new_name.txt new_name_copy.txt

ls ./
```

Copy the file to a new directory.

```bash
cp new_name.txt ../

```

Be careful, you could overwrite if there already is another file with the same name as the file you want to move or copy.


# Cleaning up

```bash

rm -i new_name_copy.txt

ls ./
```
The flag -i returns a question before removing, are you sure you want to remove teh file?

rm -i test 
rm: remove regular file 'test'? y


Remove a directory

```bash
cd ../
ls ./

ls test_dir1/
rm -i -r test_dir1/

ls ./
```

No undo button so be sure where you are and what you are deleting.
Good practice is to first list (ls) what you think of removing so you know which directories and files the command will delete.



## Using wild cards and pattern matching

```bash
mkdir my_files.txt
cd my_files.txt
touch abc.txt abc.jpg xyz.txt xyz.jpg cat.txt car.txt
```

First, list all text files.

```bash
ls *.txt
```

Then we’ll show all files with the name xyz.

```bash
ls xyz*
```

What about if we want to list all files except those with xyz in the name?
```bash
ls -Ixyz*
```

This example requires the -I flag to ls - i.e. ignore. This is one way to that, but you could also use more formal pattern matching which is more flexible and more powerful as it can be used with other commands such as mv.

```bash
ls [^x]*
```

Here we are essentially saying ‘show me everything except things that start with x’.


We can easily extend this to make it exclude objects that do not start with x or a. Like so:

```bash
ls [^xa]*
```

Or only files that start with ‘c’?

```bash
ls c*
```

Or all files where the name contains ‘c’?

```bash
ls *c*
```
Finally here’s a little example of how to use something like this with copy – i.e. the cp command. We want to copy all .txt files from dir1 to dir2. As we learned previously, cp acts much like mv, except it only copies files. You can still accidentally overwrite things though so beware!

```bash
cp *.txt ../test_dir2
ls ../test_dir2
```













