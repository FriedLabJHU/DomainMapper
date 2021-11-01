#!/bin/bash

pip install .
cd data 
dommap -f UP000000625-Escherichia_coli-K12.hmm.out -o test.map.out > out.log
cd ..
