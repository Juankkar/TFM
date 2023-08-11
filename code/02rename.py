#!/usr/bin/env python

#################
##  LIBRERIES  ##
#################
import pandas as pd
import os

################
##  VARIBLES  ##
################

## We will use some of the fields of report.tsv in metadata
report = pd.read_csv("metadata/report.tsv", sep="\t")
path_origbam = "data/original_bam/" 

submitted_ftp_list = list(report.submitted_ftp)
old_filename_list = []
for name in range(len(submitted_ftp_list)):
    names = path_origbam + submitted_ftp_list[name].split("/")[-1]
    old_filename_list.append(names)

run_accession_list = list(report.run_accession)
new_filename_list = []
for name in range(len(run_accession_list)):
    names = path_origbam + run_accession_list[name] + ".bam"
    new_filename_list.append(names)

#################
##  EXECUTION  ##
#################

for name in range(len(old_filename_list)):
    old_filename = old_filename_list[name]
    new_filename = new_filename_list[name]

    # Rename the file
    os.rename(old_filename, new_filename)