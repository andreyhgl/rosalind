#!/usr/bin/env python

import sys

def read_file(file):
  with open(file, "r") as f:
    content = f.read().strip()
  return content

def count_character(content):
  # extract the unique characters from the string, keep in alphabetic order
  chars = "".join(sorted(set(content)))

  # assign the counted chars to output
  output = ""

  # loop over each char and count, save as string w/ whitespace
  for char in chars:
    output += str(content.count(char)) + " "
  
  print(output.strip())

# get first argument
file = sys.argv[1]

# read the sequence
sequence = read_file(file)

# count the nucleotides
count_character(sequence)