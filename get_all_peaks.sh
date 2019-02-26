#!/bin/bash

#grep chr `ls *loco.txt` |awk '{if ($12<0.00000251188) {print $1, "\t", $2,  "\t", $12}}' |sed "s/_\|:/\t/g" |cut -f 1,2,3,7,8,9
 grep chr `ls *loco.txt` |awk '{if ($15<0.00000251188*5) {print $1, "\t", $2,  "\t", $15}}' |sed "s/_\|:/\t/g" |cut -f 1,2,6,7,8

