# Awk tutorial

awk is a domain-specific language designed for text manipulation. The abbreviation awk could easily stand for awkward as the syntax is a bit unintuitive. However, it actually stands for the authors Alfred Aho, Peter Weinberger and Brian Kernighan. Awk works on a line-by-line base. It goes through each line in the file and does to the line whatever the user has specified. It is very versatile and can do very different text manipulations such as computing means, combining columns, adding new columns, selecting specific lines etc. The awk command is put in '' (single quotes). Below you can see some examples.

```shell

# First let's make a folder to work in
mkdir awk

# Let's get some files
cp /home/genomics/scratch/data/martin2019/martin2019.* ./

# Let's have a look at these file (note '| column -t' makes the columns nicely aligned)
head martin2019.info | column -t
head martin2019.imiss | column -t

```

In the simplest way of running awk, we can tell awk to print out all lines where a specific condition is true. The condition, e.g. line number is 1, is put between single quotes.

```shell

# Let's get the first line
# in awk, NR means number of line
awk 'NR==1' martin2019.info
# get the tenth line
awk 'NR==10' martin2019.info

# Get only every 10th line
awk 'NR%10==0' martin2019.info

# awk can be used like grep to select all lines with a specific word such as "Htim" which is found in multiple sample names.
awk '/Htim/' martin2019.info
# this gives the same output as
grep 'Htim' martin2019.info

# Let's use awk for something else that grep cannot do
# Get all samples with more than 30% missing data
# note: $5 means fifth column
awk '$5>0.3' martin2019.imiss

# you can can also combine the previous codes to select only lines with Htim with more than 20% missing data
awk '/Htim/ && $5>0.2' martin2019.imiss

```

By using curly brackets {} we can specify more complex awk commands.

```shell
# Of the samples with missing data proportion above 0.2, print the first column ($1) which contains the individual labels
awk '{if($5>0.2) print $1}' martin2019.imiss
# with if() we can define a specific condition that needs to be fulfilled (e.g. column 5 needs to be greater than 0.01)

# We could also specify now that the letters Htim need to be in the first column and the fifth column needs to be higher than 0.2.
awk '{if($1~/Htim/ && $5>0.2) print $1}' martin2019.imiss
# Note that && means that both conditions need to be true. You could specify lots of different conditions.

# if we now wanted to print the entire line fulfilling both conditions (containing Htim and >20% missing data), we have to specify print $0 ($0 stands for the entire line)
awk '{if($1~/Htim/ && $5>0.2) print $0}' martin2019.imiss

```
Everything until now was always performed for each line. If we wanted to do something before or after reading the lines, we can use BEGIN{} and END{}, before or after the main code {}, respectively.

```shell

# Let's use END{} to get the mean missing data proportion. The part in END{} is preformed after going through all lines. a is a variable that will sum up the column five entries (missing data proportions). In the end statement, we specify that it should calculate the sum of missing data proportions divided by the number of lines.
awk '{a+=$5}END{print a/NR}' martin2019.imiss

# Actually the header is still in here, so let's skip the first line
awk '{if(NR>1) a+=$5}END{print a/(NR-1)}' martin2019.imiss

# Now if we want to get the mean missing data proportion per group
# we need to first combine the two files to get the missing data proportion and group information in a single file. As they are ordered the same way, we can just combine them with paste
paste martin2019.info.tmp martin2019.imiss > martin2019.combined
sed -i 's/ /\t/g' martin2019.combined # as the info file has blanks instead of tabs as delimiters 

# Get the mean missing data proportion per group
awk '{if(NR>1) missing[$2]+=$7; count[$2]++}END{for(g in missing){print g,(missing[g]/count[g])}}' martin2019.combined | column -t



#### String manipulation

Now that we have learned to create variables, we can also explore how to manipulate them. This is not always straightforward in bash, but again it is really worth learning how to do this as it can make script writing much more straightforward. In most cases we will use string manipulation on filenames and paths, so we will use this as an example now.

First, let's declare a variable. We'll make a dummy filename in this instance:

```shell
FILE="$HOME/an_example_file.txt"
```
Let's echo this back to the screen:

```shell
echo $FILE
```
An important point to note here is that the `$HOME` variable has been interpreted so that we now have the entire file path.

Let's say we just want the actual filename, i.e. without the directory or path? We can use `basename`:

```shell
basename $FILE
```
Alternatively, we could remove the filename and keep only the directory or path:

```shell
dirname $FILE
```
For now though, we want to operate on the filename itself, so let's redeclare the variable so it is *only* the filename.

```shell
FILE=$(basename $FILE)
echo $FILE
```
Note that here we have to wrap the `basename $FILE` command in `$()` because it is an actual command.

OK so onto some proper string manipulation. Let's remove the `.txt` suffix.

```shell
echo ${FILE%.*}
```
What did we do here? First we have to wrap the entire variable name in curly brackets - this will not work without them. The `%` denotes that we want to delete everything after the next character, which in this case is `.*` - i.e. everything after the period. Note that the following would have also worked:

```shell
echo ${FILE%.txt}
```

We don't need to limit ourselves to the suffix. We could also delete everything after the last underscore. Like so:

```shell
echo ${FILE%_*}
```
We could also set it so that we delete everything after the *first* underscore:

```shell
echo ${FILE%%_*}
```
 We can also delete from infront of the characters in our string manipulation example. For example:

```shell
echo ${FILE#*.}
```
This deletes everything up to and including the period character. We could also do the same with the underscores:

```shell
echo ${FILE#*_}
echo ${FILE##*_}
```

Where again, a single `#` states we want to delete only after the last occurrence and a double `##` denotes we want to delete everything after the first occurrence.

You might be wondering, what exactly is the point of this? Well altering filenames is very important in most bioinformatics pipelines. So for example, with simple string manipulation you can change the suffix of a filename quickly and easily:

```shell
echo $FILE
echo ${FILE%.*}.jpg
```
One last point here; string manipulation in bash is not straightforward. It takes a lot of practice to get right and remember properly. We google [this excellent tutorial](https://www.tldp.org/LDP/abs/html/string-manipulation.html) nearly all the time!

#### Bash control flow

[Control flow](https://en.wikipedia.org/wiki/Control_flow) is an important part of many different programming languages. It is essentially a way of controlling how code is carried out.

Imagine you have to perform the same operation on many different files - do you want to type out a command for each and everyone of them? Of course not! This is why you might use control flow to repeat a command multiple times. There are many different types of control flow, but for now we will focus on the most common one - a `for` loop.

Let's have a look at a simple example:

```shell
for i in {1..10}
do
 echo "This is $i"
done
```

All this is is saying is that for each number between 1 and 10, echo a "This is 1", "This is 2" and so on to the screen. `do` and `done` initiate and stop the loop respectively.  Here, the variable `i` is used within the loop but this is completely arbitrary - you can use whatever variable you would like. Indeed, it is often much more convenient to use a variable that makes sense to you. For example:

```shell
for NUMBER in {1..10}
do
echo "This is $NUMBER"
done
```
It doesn't just have to be numbers either. You can use a loop to iterate across multiple strings too. For example:

```shell
for NAME in Mario Link Luigi Peach Zelda
do
echo "My name is $NAME"
done
```

Of course, this is a silly example, but you could easily substitute this with filenames - making it quite clear why control flow is an essential skill for effective bash programming in bioinformatics.

#### Declaring arrays

Imagine we want to run a `for` loop on some text files. We'll make five of them to demonstrate - and we can actually do this with a `for` loop too:

```shell
for i in {1..5}
do
touch file_${i}.txt
done
```
Use `ls` after running this code and you'll see five text files.

Now, what if we want to do something simple like go through all of them and print their names to the screen? We could do it like this:

```shell
for FILE in *.txt
do
echo $FILE
done
```
This works really well for this simple example. However as your code becomes more advanced, it is easy for something like this to become quite dangerous. Imagine for example that in our `for` loop, we create a new `.txt` file each time? We would be in danger of creating an infinite loop, that continually prints the names of the new files it creates.

For this reason, it is **best practice** in bash to use **arrays**. These are essentially predefined lists of variables. They are easy to make too. Let's try a simple example.

```shell
ARRAY=(Link Zelda Gannon)
```
Now we can try printing this to the screen:

```shell
echo $ARRAY
```
This only prints the first value of our array. Actually, arrays have indexes, so we can print any value we specify like so:

```shell
echo ${ARRAY[0]}
echo ${ARRAY[1]}
echo ${ARRAY[2]}
```
Notice that like python (and unlike R) everything in bash is zero-indexed - i.e. the first variable is zero and so on.

What if we want to print everything in the array?

```shell
echo ${ARRAY[@]}
echo ${ARRAY[*]}
```
Either of these will work fine.

We can also loop through the array, like so:

```shell
for CHARACTER in ${ARRAY[@]}
do
echo $CHARACTER
done
```

The purpose of the array here is that it ensures the **scope** of our loop is limited and that it doesn't get carried away, operating on things it shouldn't do.

One last point about arrays - it is often quite cumbersome to define them by hand. Imagine if you wanted to make an array for hundreds of files? Luckily you can also *declare* them from for loops too:

```shell
ARRAY2=($(for i in *.txt
do
echo $i
done))
```

This doesn't look very neat though - you can actually write a for loop like this on a single line - i.e.:

```shell
ARRAY2=($(for i in *.txt ;do echo $i; done))
```
Where `;` indicates a separate line (as seen in the above example).

The convenience of arrays, like much of this tutorial will become much more apparent as you become more experienced in using Unix for bioinformatics.

### Writing a bash script

So far we have learned a lot about bash as a programming language - but can we use it to write a program? Well actually... yes! This is extremely easy and it is exactly what we set out to do when we write a bash script.

Let's start with a really basic example. Type `nano` into the command line in order to open the `nano` text editor.

Then we can write a simple bash script, like so:

```shell
#!/bin/sh

# a simple bash script
echo "Hello world"

exit
```

Save it as `my_first_script.sh`. You can take another look at the output with `cat` or `less` if you want to check you saved it properly.

Let's breakdown some of the script. First of all there is this `#!/bin/sh` line. You don't need to worry too much about that - it's just good practice to ensure the script is run in the bash language. We also have another line starting with `#` - this is just a comment. Here it explains something about the script. Comments are really important actually and you should fill your script with them - they are a good way of letting yourself know what you have done. Again, they can be invaluable when you come back to your scripts after some time away...

Now we can actually run the script. We do that like so:

```shell
sh my_first_script.sh
```

You just wrote your first program! [What it said is also quite relevant](https://en.wikipedia.org/wiki/%22Hello,_World!%22_program).

You can also write a script that is interactive. Let's write another script to take an input from the command line.

Open `nano` and create a script with the following:

```shell
#!/bin/sh

# a simple bash script with input
echo "My name is ${1} and my best friend is ${2}"

exit
```

Save it as `name_script.sh`. In this script, the `${1}` and `${2}` variables are just specifying that this script will take the first and second arguments to the script from the command line. Let's see it in action.

```shell
sh name_script.sh Mark Milo
```

Feel free to add whatever combination you want in here. You can actually try running this without the arguments too and see that it still works, just the output doesn't make much sense.

### A more serious scripting example

Now we have learned a little about how to script, let's write one that will do something for us. We'll create five files (again using a `for` loop) and then convert them all from `.txt` to `.jpg`.

Firstly, let's make those files:

```shell
for i in {1..5}
do
touch file_${i}.txt
done
```

Now, we can open up `nano` and write out script. We would do this like so:

```shell
#!/bin/sh

# a script to rename files

# declare an array
ARRAY=($(for i in file*.txt; do echo $i; done))

# loop over array
for FILE in ${ARRAY[@]}
do
echo "Creating ${FILE%.*}.jpg"
mv $FILE ${FILE%.*}.jpg
done
```
Now write this script out as a `file_renamer.sh`.

If you run this as `sh file_renamer.sh` - it will print the name of each file it converts to the screen. You can then use `ls` to see that it has indeed converted all the `.txt` files to `.jpg`.

### A script writing challenge

With all the skills we have learned in this tutorial, it is now time for you to put them to the test. Return to the `unix_exercises` directory you created when you used `git clone` and write a short script to do the following:

* make an array of the three files
* loop through the array and count the number of lines in the file
* print the number of lines and the name of the file to the output

<details><summary>Click here to see a possible solution.</summary>
<p>


```shell
#!/bin/sh

# a possible solution script

# declare an array
ARRAY=($(for i in *.t*; do echo $i; done))

# loop over array
for FILE in ${ARRAY[@]}
do
echo "${FILE}"
wc -l ${FILE}
done
```

</p>
</details>
