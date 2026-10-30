#!/usr/bin/env python

import argparse
import re

# This script assumes single end reads and a sorted SAM file.
# Here is the bash command needed to sort you SAM file: [samtools sort -o <sorted.sam> <input.sam>]

def get_args():
    parser = argparse.ArgumentParser(
                                    description="Program to deduplicate PCR samples in a sorted SAM file with single end reads",
                                    epilog=(
                                            "Assumptions:\n"
                                            " - Input SAM file is sorted by chromosome AND position\n"
                                            " - Sorted SAM file reads are single end reads\n"
                                            " - UMIs are in the QNAME (last field and colon seperated) and match the UMI list provided"),
                                            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--infile", help="Path to SORTED SAM file", type=str)
    parser.add_argument("-o", "--outfile", help="Path to SORTED SAM Output File", type=str)
    parser.add_argument("-u", "--umi", help="File containing the UMI list", type=str)
    return parser.parse_args()
args = get_args()

#using this to call in my which strand function
def is_revcomp(flag: int) -> bool:
    '''when given a bitwise flag this function will check the sixteenth
     bit to see if it is set'''
    return (flag & 16) == 16

def which_strand(flag: str)-> str:
    '''Strand determination based on rev_comp determination'''
    if is_revcomp(flag):
        return "-"
    else:
        return "+"
    
def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    '''when given a cigar string function will utilize regex and output the cigar
    separated by all numbers (length) before a single upper case letter (operation) as a pair
    ex) input: 40M25N5M output: [('5', 'S'), ('25', 'M')]'''
    pattern = r"(\d+)([MIDNS=X])" #reg ex to grab CIGAR pairs as mentioned above
    parsed_cigar = [] #make a list of length and letter tuples built IN ORDER they're seen
    for number, letter in re.findall(pattern, cigar): #findall returns digits and letter
        length = int(number) #convert the number from a str to an int to it can be used for math
        parsed_cigar.append((length, letter)) #store the pair  keeping the order
    return parsed_cigar 

def true_start(pos: int, cigar: str, strand: str ) -> int:
    '''Determine the true start site of a strand by checking strandedness, then grabbing pos in col 4 of the read
    and looking at CIGAR string for any soft clipping. 
    PLUS: subtract that number from col 4 pos to output the true 5' start site. 
    MINUS: add the reference consuming lengths back and right soft clip to pos, then subtract 1'''
    options = parse_cigar(cigar) #list of tuples holding the CIGAR string info
    first_len, first_op = options[0] #leftmost CIGAR string info 5' end on PLUS
    last_len, last_op = options[-1] #rightmost CIGAR string info 5' end on MINUS
    if strand == "+": #plus strand: then 5' is on the leftmost side of alignment
        if first_op == "S": #if there is soft clipping present at the 5' end
            return pos - first_len #shift start to the left where it would have begun preclipping
        return pos #if not clipped true start alr there so return either way
#THIS IS FOR THE MINUS STRAND: 
    end = pos #start from the leftmost pos
    for length, op in options:
        if op in "MDN=X": # if there is anything consuming the reference
            end += length #move the end by adding that val
    if last_op == "S": #if theres soft clipping at the rightmost (5' end) of the minus read
        end += last_len #add that length to the end to get the true start
    return end - 1 #end is ONE past the last base so sub 1

def get_umi(qname: str) -> str:
    '''determine the UMI from the QNAME'''
    divide_qname = qname.split(":") #split the QNAME by the colons
    umi = divide_qname[7] #grab the UMI from the list
    return umi

def is_mapped(flag: int) -> bool:
    '''when given a bitwise flag this function will check the fourth
    bit to see if it is set and return if it is mapped(!=4)'''
    return (flag & 4) != 4


## Initalizing sets and counters:

#sets:
seen = set()
umi_set = set()
#counters:
headerline_count = 0
unmapped_count = 0
bad_umi_count = 0
pcr_dupe_count = 0
unique_reads_count = 0

## Populating my umi set with umi file!

with open(args.umi, "r") as umi_fh:
    for line in umi_fh:
        fields = line.strip("\n").split("\t")
        umi_set.add(fields[0])

## Starting the deduping process!
with open(args.infile, "r") as fh, open(args.outfile, "w") as ofh: #opening input and output files defined in args
    for line in fh: #for each line in my input file
        if line.startswith("@"): #if it starts with an @ 
            header_info = line #IT IS A HEADER LINE
            headerline_count += 1 # so count it 
            ofh.write(header_info) #and write the header lines out to the outfile
            continue #MOVE ON! once we have reads
        read = line.split("\t") #split the reads to make them a list
        chromo = read[2] #chromosome col
        pos = int(read[3]) #pos col had to wrap in int to do math later
        flag = int(read[1]) #bitwise flag col had to wrap in int to call in a function as an int
        cigar = read[5] #CIGAR string col
        qname = read [0] #need this for get UMI function
        if not is_mapped(flag): #skips unmapped reads
            unmapped_count += 1 #add to my unmapped counter
            continue #after counting it move ON!
        umi = get_umi(qname) #calling my get umi func on the qname
        if umi not in umi_set: #if there is an unknown umi
            bad_umi_count += 1 #add to the bad umi count and move on
            continue #MOVEEE ON NO BAD UMIS
        strand = which_strand(flag) #check which strand you are on
        five_p_start = true_start(pos, cigar, strand) #find the true 5' start by doing some math with my function (:)
        key = (umi, chromo, five_p_start, strand) #make allll the criteria for determining a PCR dupe a list and store in mem
        if key in seen: #check if the combo of all those deets has been seen before and if it has then == PCR DUPE! GET OUTTA HERE
            pcr_dupe_count += 1 #add the the counter
            continue #then move onnnnnnnn
        seen.add(key) #if not seen add to my set please!
        ofh.write(line) #write out any unique mapped reads to my outfile
        unique_reads_count += 1 #and update my counter for unique reads

## Print Statements for final counts!
print(f"The number of header lines is: {headerline_count}\n\
        The number of unique reads: {unique_reads_count}\n\
        The number of PCR dupes seen: {pcr_dupe_count}\n\
        The number of unmapped reads: {unmapped_count}\n\
        The number of reads with bad umis: {bad_umi_count}\n")



#testing print statements ignore
    # print(f"This is my chr:{chromo}\n\
    #         This is my pos:{pos}\n\
    #         This is my bitflag:{flag}\n\
    #         This is my cigar:{cigar}\n\
    #         This is my qname:{qname}\n\
    #         This is my UMI:{umi}\n\
    #         This is my strand:{strand}\n\
    #         This is my cigar string parsed:{parse_cigar(cigar)}\n\
    #         This is my true 5' start pos:{five_p_start}\n")
