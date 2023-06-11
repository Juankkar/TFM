#!/usr/bin/env python

#################
##  LIBRERIES  ##
#################

import os

################
##  VARIBLES  ##
################

old_filename_list = ["data/original_bam/C2466_GCCAAT.bwa.sorted.rmdup.recal.realigned.fixed.bam",
                     "data/original_bam/C2467_CGATGT.bwa.sorted.rmdup.recal.realigned.fixed.bam",
                     "data/original_bam/C2468_TTAGGC.bwa.sorted.rmdup.recal.realigned.fixed.bam",
                     "data/original_bam/C3481_GCCAAT.bwa.sorted.rmdup.recal.realigned.fixed.bam",
                     "data/original_bam/C3852_CTTGTA.bwa.sorted.rmdup.recal.realigned.fixed.bam",
                     "data/original_bam/C3890_TGACCA.bwa.sorted.rmdup.recal.realigned.fixed.bam",
                     "data/original_bam/NG-5502_C2463.bwa.sorted.rmdup.recal.realigned.sorted.fixed.bam",
                     "data/original_bam/NG-5502_C2461.bwa.sorted.rmdup.recal.realigned.sorted.fixed.bam",
                     "data/original_bam/C-3861_GCCAAT_L007.bwa.sorted.rmdup.recal.realigned.fixed.bam",
                     "data/original_bam/C2213_CGATGT.bwa.sorted.rmdup.recal.realigned.fixed.bam",
                     "data/original_bam/C896_TGACCA.bwa.sorted.rmdup.recal.realigned.fixed.bam",
                     "data/original_bam/C3988_ACAGTG.bwa.sorted.rmdup.recal.realigned.fixed.bam"]

new_filename_list = ["data/original_bam/ERR696683.bam",
                     "data/original_bam/ERR753368.bam",
                     "data/original_bam/ERR753369.bam",
                     "data/original_bam/ERR753370.bam",
                     "data/original_bam/ERR753371.bam",
                     "data/original_bam/ERR753372.bam",
                     "data/original_bam/ERR753373.bam",
                     "data/original_bam/ERR753374.bam",
                     "data/original_bam/ERR753375.bam",
                     "data/original_bam/ERR753376.bam",
                     "data/original_bam/ERR753377.bam",
                     "data/original_bam/ERR753378.bam"]

#################
##  EXECUTION  ##
#################

## rename the files
if len(old_filename_list) != len(new_filename_list):
    print("Error: The number of BAM filenames does not match the number of new filenames.")
    exit(1)

for i in range(len(old_filename_list)):
    old_filename = old_filename_list[i]
    new_filename = new_filename_list[i]

    # Rename the file
    os.rename(old_filename, new_filename)
