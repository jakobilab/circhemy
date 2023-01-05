#!/bin/bash

# step one: remove circatlas names for this step
# sed will mess up the columns because we have identical names and locations in case of non-gene circRNAs

# cut -f2 hg38_hg19_v2.0.txt > /tmp/foo

awk 'OFS="\t" {if($4 != "-") print $1,$2}' $1 > /tmp/foo

# hg19
awk 'OFS="\t" {($1="");($2="");($3="hg19"); if($4 != "-") print $0}' $1 | sed 's/:/\t/g' | sed 's/|/\t/g' | awk 'OFS="\t" {$3=$3-1; print $0}' > /tmp/boo


paste /tmp/foo /tmp/boo > rn5.v1
awk 'OFS="\t" {print "homo_sapiens",$0}' rn5.v1 | tail -n +2 > rn5.v2


awk 'OFS="\t" {if($3 != "-") print $1,$2}' $1 > /tmp/foo


# hg38
awk 'OFS="\t" {($1="");($2="");($4="hg38"); if($3 != "-") print $0}' $1 | sed 's/:/\t/g' | sed 's/|/\t/g' | awk 'OFS="\t" {$3=$3-1; print $4,$1,$2,$3,$5,$6,$7,$8}' > /tmp/boo


paste /tmp/foo /tmp/boo > rn6.v1
awk 'OFS="\t" {print "homo_sapiens",$0}' rn6.v1 | tail -n +2 > rn6.v2
