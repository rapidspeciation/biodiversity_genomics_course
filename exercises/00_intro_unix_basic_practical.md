
# Introduction to command line

## Log in to server

```bash
ssh <your_username>@168.176.34.122
ssh perseus

cd /scratchsan/C_computacion/<your_username>

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

## Make an intro directory
You have seen a list representing the directories and files in this directory. It is quite empty now.

Make a directory for this intro session with the command

mkdir 

Remember no spaces! No special characters!

```bash
# first make sure you are in the right place 
pwd

# you should be in /scratchsan/C_computacion/<your_username>

mkdir intro 
```

Does everyone have a directory with the correct name?

```
ls .

```

Go to that directory

```bash
cd intro

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

Check out some of the list options.

Make a -l long list (more information on each file), -h human readable, -t time sorted

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

Change directory

```bash
cd test_dir1

ls ./
```


```bash
touch file3.txt
touch file4.txt

ls ./
```

Go back to your first directory
```bash
cd ../

pwd
ls ./
```


Relative paths

Tell the path relative to the directory you are in.

./ is the current directory

../ is the directory one level above

```bash

cd test_dir2

pwd

ls ./
# this is the same as the absolute path:
ls /scratchsan/C_computacion/<your_username>/intro/test_dir2

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

touch file5.txt
cat file5.txt 

```

Well it is empty, so we need some content.

Here I will use nano is one text editor, but there are others like vi, vim, emacs, etc.

Open the file with nano and write something.

Save and exit with Ctrl + X (same as ^X)

```bash

nano file5.txt
cat file5.txt 

```
If we do not want the whole file printed to screen we can use the command less.

```bash
less file5.txt

```
You can search in less with /text

n next match

N previous match

Quit with q


To take a look at the first lines of a file, we can use head. The default is 10 lines but you can decode the number of lines with the flag -n
```bash

head file5.txt

head -n3 file5.txt
```

Take a look at the last lines of the file with tail.
```bash

tail file

tail -n3 file
```


# Moving and copying

Command (cp or mv) source target

Source is the file (or directory) you want to copy or move.

Target is the file (or directory) you want to copy or move it to.

```bash

mv file5.txt ../

```
mv can also be used to rename files

```bash
touch test_file2.txt
mv test_file2.txt file6.txt

```
Be careful you can overwrite if there is a file with the same name where you are moving the file to.

We can also move directories:

```bash
#make directories
mkdir test1 test2

#move directories
mv test1/ test2/
ls ./
ls test2/

```



Copy
```bash

cp file6.txt file6_copy.txt

ls ./

```

Copy the file to a new directory.

```bash
cp file.txt ../

```

Copy a directory, but using the option -r (recursive).

```bash
cp -r test2/ ../

```

Be careful, you could overwrite if there already is another file or directory with the same name as the file you want to copy.


# Cleaning up

Removing files.

```bash

rm -i new_name_copy.txt

ls ./
```
The flag -i returns a question before removing, are you sure you want to remove the file?

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

To remove an empty directory:

```
mkdir to_remove

rmdir to_remove

```

There is no undo button so be sure where you are and what you are deleting.
Good practice is to first list (ls) what you think of removing so you know which directories and files the command will delete.


## Using wild cards and pattern matching

```bash
mkdir my_files
cd my_files
touch abc.txt abc.jpg xyz.txt xyz.jpg cat.txt car.txt
```

First, list all text files.
The asterisk * match any character and any number of characters.
For example, *.txt match anything that ends with .txt.

```bash
ls *.txt
```

Then we can list all files that starts with the letters xyz.

```bash
ls xyz*
```

What about if we want to list all files except those that starts with x?
This is one way to that, by using a formal pattern matching which is flexible and powerful as it can be used with other commands such as mv.

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













