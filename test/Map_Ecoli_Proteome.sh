#!/bin/bash

# This example starts before the creation of an hmmscan full output
# For an example of parsing your own hmmscan full output skip to part 6
ecoli_proteome="uniprot-proteome_UP000000625.fasta"
ecoli_proteome_tar="../data/Uniprot\ Proteomes\ \(canonical\ and\ isoforms\)/uniprot-proteome_UP000000625.fasta.gz"

tar -xvzf ${ecoli_proteome_tar} -C .

# Download the ECOD hmm profiles
wget http://prodata.swmed.edu/ecod/distributions/ecodf.hmm.tar.gz

tar -xvzf ecodf.hmm.tar.gz -C .

# Create HMMR3 binaraies
hmmpress ecodf.hmm

# Crearte hmmscan full output
hmmscan -o Ecoli.hmm.out ecodf.hmm ${ecoli_proteome}

# Mapping Protein Domains
python src/DomainMapper.py -f Ecoli.hmm.out -o Ecoli.mapped.out

# Compare the output of this mapping to the data file for E. coli in data/hmmscan Full Output/UP000000625-Escherichia_coli-K12.hmm.out.tar.gz