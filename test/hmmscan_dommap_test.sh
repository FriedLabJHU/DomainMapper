#!/bin/bash -i

#####################
# DomainMapper Test #
#####################

# NOTICE:
# Run this test after installing HMMER3 and DomainMapper

# Setting Variables
ecoli_fasta="UP000000625_83333.fasta"
ecoli_hmm_out="Ecoli.hmm.out"
ecoli_dommap_out="Ecoli.dommap.tsv"

# Downloading ECOD HMM Profiles
ecod_hmm_link="http://prodata.swmed.edu/ecod/distributions/ecodf.hmm.tar.gz"
wget $ecod_hmm_link

# Untar & Unzip
ecod_hmm_tar="ecodf.hmm.tar.gz"
tar -xzvf $ecod_hmm_tar

# Pressing ECOD HMM Profiles with HMMER3
ecod_hmm="ecodf.hmm"
hmmpress $ecod_hmm

# Scanning E. coli proteome
# The number of CPUs can be increased to improve speed
hmmscan --cpu 4 -o ${ecoli_hmm_out} ${ecod_hmm} ${ecoli_fasta}

# Running DomainMapper
dommap -f ${ecoli_hmm_out} -o ${ecoli_dommap_out}