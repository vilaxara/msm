#!/bin/bash

gawk 'NR==FNR{a[$1][$0];next} $0 in a {for (i in a[$0]) print i}' cc.txt S_1.txt > temp

sort -nr -k2 temp > s_1.txt

gawk 'NR==FNR{a[$1][$0];next} $0 in a {for (i in a[$0]) print i}' cc.txt S_2.txt > temp

sort -nr -k2 temp > s_2.txt

gawk 'NR==FNR{a[$1][$0];next} $0 in a {for (i in a[$0]) print i}' cc.txt S_3.txt > temp

sort -nr -k2 temp > s_3.txt

gawk 'NR==FNR{a[$1][$0];next} $0 in a {for (i in a[$0]) print i}' cc.txt S_4.txt > temp

sort -nr -k2 temp > s_4.txt

gawk 'NR==FNR{a[$1][$0];next} $0 in a {for (i in a[$0]) print i}' cc.txt S_5.txt > temp

sort -nr -k2 temp > s_5.txt

gawk 'NR==FNR{a[$1][$0];next} $0 in a {for (i in a[$0]) print i}' cc.txt S_6.txt > temp

sort -nr -k2 temp > s_6.txt
