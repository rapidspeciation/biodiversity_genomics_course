#! /usr/bin/env python

# python version 3
# by Joana, script to convert vcf to phylip
# Some functions are recycled from the RAD python script written by Sam Wittwer

import sys, getpass, re, argparse, gzip

# Define the neccessary functions:

# Function to extract header info: extractHeaderInfo
# takes a vcf(.gz) file from line 1, goes through header, returns:
# 0[string]             header
# 1[list[string]]       individual IDs
# 2[int]                number of individuals in vcf file
# 3[int]                length of header
def extractHeaderInfo(input):
        linecounter = 0
        header = ""
        for line in input:
                linecounter += 1
                if re.match("^\#\#", line):
                        header+=line
                else:
                        header+=line
                        n_individuals = len(str(line.strip('\n')).split('\t'))-9
                        id_individuals = str(line.strip('\n')).split('\t')[9:]
                        break
        return header, id_individuals, n_individuals, linecounter

# Function to make all sample names of equal length (could be simplified)
# FillUp: accepts list, fills up list with optional character (default is " ") until they are of equal length
def fillUp(list, fill = " "):           # Checks length of entries in a list, fills up with specified fill string.
        returnlist = []
        for entry in list:
                if len(entry) < len(max(list, key=len)):
                        returnlist.append(entry+(len(max(list, key=len))-len(entry))*fill)
                else:
                        returnlist.append(entry)
        return returnlist

        
# Function to write the the lines in phylip format        
def writePhylipSequences(samplenames, sequences, outputdestination, writeref):
        if writeref:
                beginning = 0
                nsamples = str(len(samplenames))
        else:
                beginning = 1
                nsamples = str(len(samplenames)-1)
        nbases = str(len(sequences[0]))
        outputdestination.write(nsamples +"\t"+nbases+"\n")
        outstring = ""
        for i in range(beginning, len(samplenames)):
                outstring += samplenames[i]+"".join(sequences[i])
                outstring +="\n"
        outputdestination.write(outstring.strip("\n"))
        

# Parse the arguments provided
parser = argparse.ArgumentParser(description='Convert vcf file to phylip file format')

parser.add_argument('-i', '--input', dest='i', help="input file in vcf(.gz) format [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="output file [required]", required=True)
parser.add_argument('-r', '--ref', action='store_true', help="if -r is specified, the reference sequence will be included in the phylip file)", default=False)
parser.add_argument('-f', '--fill', action='store_true', help="if -f is specified, all sites not in the vcf file will be printed as missing (N)", default=False)
parser.add_argument('-e', '--exclIndels', action='store_true', help="if -e is specified, indels are not printed (else replaced by N)", default=False)
parser.add_argument('-m', '--mtDNA', action='store_true', help="if -m is specified, haploid genotype calls are expected in the vcf", default=False)

args = parser.parse_args()
        
# Set the default values:
if args.i.endswith('.gz'):
        input = gzip.open(args.i,'rt')
else:
        input = open(args.i,'r')
output = open(args.o,'w')
writeref = args.ref
fill = args.fill        
noIndels = args.exclIndels
haploid = args.mtDNA
prev=100000000000000000000000000000000000

# How to convert vcf-style genotypes to single letters
AmbiguityMatrix = [["A","M","R","W","N"],["M","C","S","Y","N"],["R","S","G","K","N"],["W","Y","K","T","N"],["N","N","N","N","N","N"]]
CoordinatesDictionary = { "A":0 , "C":1, "G":2, "T":3, ".":4 }
GetGenotype = lambda individual, altlist: AmbiguityMatrix[CoordinatesDictionary[altlist[int(individual[0])]]][CoordinatesDictionary[altlist[int(individual[1])]]]

# If -f and -e are specified 
if noIndels and fill:
        print("-f and -e are incompatible! Please decide if you want all sites or not")
        sys.exit(2)

# Get the header info
headerinfo = extractHeaderInfo(input)        #skips header, retains IDs in headerinfo

# Get the sample labels
IDs = []        #IDs holds all sample names without equal spacing 
IDs.append("reference ")
resultsequences = []        #resultsequences holds all complete sequences to be written
resultsequences.append([])                #resultsequences[0] holds reference
for entry in headerinfo[1]:
        IDs.append(entry.replace("-",".")+" ")
        resultsequences.append([])
samplenames = fillUp(IDs)
        
linecounter = 0
print("\ngenerating phylip file with ",len(samplenames)-1," individuals")

# Go through the lines to get the genotypes
for line in input:
        site = line.strip('\n').split('\t')
        indel = False
        pos = int(site[1])
        
        # If missing positions should be filled up with Ns (-f specified, e.g. sites of low quality that were filtered out)
        if pos > (prev + 1) and fill:
                addLine=pos-(prev+1)
                individualcounter = 1
                for individual in site[9:]:
                        resultsequences[individualcounter]+= "N" * addLine
                        individualcounter += 1
                linecounter += addLine
                resultsequences[0] += "N" * addLine  # reference

        # site contains a deletion, replace by missing data
        if len(site[3])>1:
                indel=True

        else:
                alleles = site[4].split(",")  # if there are more than 1 alternative alleles, they are separated by commas
                for alt in alleles:
                        if len(alt)>1 or '*' in alt:  # in case of an insertion
                                indel=True
                                break

        alternativeslist = []
        alternativeslist.append(site[3])
        if site[4] != ".":
                for entry in site[4].split(","):
                        alternativeslist.append(entry)
        individualcounter = 1

        # If the site is not an indel, i.e. if it is a SNP or monomorphic site
        if not indel:
                for individual in site[9:]:
                        indGeno = individual.split(":")
                        if haploid:
                                if '.' in indGeno[0]:
                                        resultsequences[individualcounter]+="N"
                                else:
                                        resultsequences[individualcounter] += alternativeslist[int(indGeno[0])]
                        else:
                                if '.' in indGeno[0]:
                                        resultsequences[individualcounter]+="N"
                                else:
                                        resultsequences[individualcounter] += GetGenotype(re.split('/|\|',indGeno[0]), alternativeslist)
                        individualcounter += 1
                linecounter += 1
                resultsequences[0] += site[3]


        # If the site is an indel and noIndels is not specified, print as missing data (else not printed)
        elif not noIndels:
                for individual in site[9:]:
                        resultsequences[individualcounter]+= "N"
                        individualcounter += 1
                linecounter += 1
                resultsequences[0] += site[3][:1]  # reference

        prev=int(site[1])

input.close()
        
writePhylipSequences(samplenames, resultsequences, output, writeref)
output.write("\n")
output.close()


