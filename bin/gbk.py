#!/usr/bin/env python

from datetime import datetime
from Bio import Entrez
import sys

def read_file(file):
  with open(file, "r") as f:
    content = f.read().strip()
  return content

def entrez_search(organism, start_date, end_date):
  # pad dates with zero if needed
  # strptime creates a datetime object
  start_date = datetime.strptime(start_date, "%Y/%m/%d")
  #print(start_date)
  # strftime creates a string
  start_date = start_date.strftime("%Y/%m/%d")
  #print(start_date)

  # parse the query
  term = '"' + organism + '"' + "[Organism]" + " AND " + '"' + start_date + '"' + "[Publication Date]" + " : " + '"' + end_date + '"' + "[Publication Date]"
  #print(term)
  Entrez.email = "dummy@aces.su.se"
  handle = Entrez.esearch(db="nucleotide", term=term)
  record = Entrez.read(handle)
  print(record["Count"])
  #out = pd.Series(record)
  #print(out)

# get arguments
file = sys.argv[1]

# parse file content
content = read_file(file).split("\n")
organism, start_date, end_date = content

entrez_search(organism, start_date, end_date)