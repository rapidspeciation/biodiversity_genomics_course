# Some more UNIX

## Gleaning files

### head

head will display the ‘head’ of a file - i.e. the first 10 lines by default. This command is essential if you are going to be working with bioinformatic data as often your files are millions of lines long and you just want to have a quick peek at what it contains.


For example, we can use head on the ilMelLudo1.reduced.fa file to see the first 10 lines

```bash
head ilMelLudo1.reduced.fa
```

We can also specify exactly how many lines we want to display. For example, if we want to see 20 lines:

```bash
head -n20 ilMelLudo1.reduced.fa
```

### tail

tail is much the same as head, except it operates at the other end - i.e. it shows you the last ten lines of a file. It is also useful for skipping the start of a file.

First of all, let’s look at the last 10 lines of the udhr.txt file.

```bash
tail ilMelLudo1.reduced.fa
```

Or the last 20?

```bash
tail -n20 ilMelLudo1.reduced.fa

```

If we use tail with + flag, we can skip lines from the start of the file. For example:

```bash
tail -n+3 ilMelLudo1.reduced.fa
```

The -n+3 argument skips the first three lines of the file. This is very useful for removing lines you are not interested in.

## Redirecting output

Let’s first just extract the first 10 lines of the declaration using head and redirect the output to a file


```bash
head ilMelLudo1.reduced.fa > my_file.txt
```

We could add the last 10 line to the file
Appending them to file with >>

```bash
tail ilMelLudo1.reduced.fa >> my_file.txt
```





## Counting 

wc you can count the number of words, number of characters and number of lines 

```bash
wc -c ilMelLudo1.reduced.fa

wc -l ilMelLudo1.reduced.fa

wc -c ilMelLudo1.reduced.fa
```




## Searching, subsetting and making changes to files

### grep
```bash

grep ">" ilMelLudo1.reduced.fa > chromosome_file.txt
```

Have a look at the file we created. 

### sed 
Stream editor

```bash
sed -n 3p chromosome_file.txt
```

In this command, 3p is just telling the -n flag we want to see the third line. We could also extract lines 3-5 like so:

```bash
sed -n 3,5p chromosome_file.txt
```

But sed can actually do much more than this. For example, it can replace text. Let’s replace all instances of “SUPER” with another name, like “MelLudo1”:


```bash
sed 's/SUPER/MelLudo1/g' chromosome_file.txt
```

This is just a small demonstration of what it is possible to do with sed. It is a very useful tool, especially for file conversion and well worth getting more familiar with.




