# DomainMapper
DomainMapper is an HMM parser designed to annotate proteins' domain structures from sequence alone using the ECOD database.  
DomainMapper annotates non-contiguous, insertional, and circularly permuted domains as well.

To begin using DomainMapper, clone the repository to a destination on your local computer

0) git clone https://github.com/FriedLabJHU/DomainMapper

To use DomainMapper, the full output from an HMMR3 hmmscan is required as an input.  
Here is how to generate this file:
1) Get a file containing protein sequences in .fasta format (for instance, a proteome file from uniprot), save this in ./DomainMapper
2) Obtain the ECOD HMM profile database at http://prodata.swmed.edu/ecod/distributions/ecodf.hmm.tar.gz
3) Place this file in the same directory as DomainMapper, and unzip (tar -xzvf ecodf.hmm.tar.gz)

At this point, you will need to run an hmmscan.  The user needs to have separately installed HMMER3.
HMMER3 is pretty easy to install. On Mac, brew install hmmer. On Ubuntu, apt install hmmer. 
Visit http://hmmer.org/documentation.html for details.
4) hmmpress ecodf.hmm
5) hmmscan -o DomainMapper/your_hmmscan_output.hmm.out ecodf.hmm DomainMapper/your_fasta_file.fasta

6) python DomainMapper/src/DomainMapper.py -f your_hmmscan_output.hmm.out -o your_hmmscan_output.mapped.out

To read the various additional options available, read the docstring with python DomainMapper.py -h.

*If you do not have wget:
