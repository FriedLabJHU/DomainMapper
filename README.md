# DomainMapper
DomainMapper is an HMM parser designed to annotate proteins' domain structures from sequence alone using the ECOD database.  
DomainMapper annotates non-contiguous, insertional, and circularly permuted domains as well.  
DomainMapper works on Linux / MacOS / Windows  

## Dependencies

In order to use DomainMapper you will need ```numpy```, ```scipy``` and  ```bio```
These can be installed using ```pip```

##### Installing Dependencies
```
# installing numpy
pip install numpy

# installing scipy
pip install scipy

# installing BioPython
pip install bio
```

## Installing DomainMapper
```
cd /DomainMapper
pip install .

# or

pip install ./DomainMapper
```

## Running DomainMapper

``` dommap -f raw_hmmscan_output.hmm.out -o mapped_protein_domains.mapped.out```

## Documentation

```
usage: dommap [-h] [-f F] [-o O] [--ecod_domains ECOD_DOMAINS] [--intra_gap INTRA_GAP] [--inter_gap INTER_GAP] [--overlap OVERLAP] [--eval_cutoff EVAL_CUTOFF] [--update]

arguments:
  -h, --help            show this help message and exit
  -f F                  Input path to file from 'hmmscan'
  -o O                  Output path for mapped domains
  --ecod_domains ECOD_DOMAINS
                        Path to ECOD 'Latest Domains' text file (default = file is automatically downloaded [165 MB
                        Free Space Required (deleted after parsing)] [2 MB File Saved])
  --intra_gap INTRA_GAP, --intra_domain_gap_tolerance INTRA_GAP
                        Optional gap size between HMM sequence and query sequence for non-contiguous alignment within
                        a domain (default = 30)
  --inter_gap INTER_GAP, --inter_domain_gap_tolerance INTER_GAP
                        Optional gap size between two domains sequences for non-contiguous merging (default = 30)
  --overlap OVERLAP, --domain_overlap_tolerance OVERLAP
                        Optional overlap between HMM domain sequence and fasta aligment in consecutive or split
                        domains (default = 40)
  --eval_cutoff EVAL_CUTOFF
                        Optional upper bound tolerance of the E-value (default = 1e-5)
  --update              Update ECOD 'Latest Domains'
```

## HMMR3 & DomainMapper Tutorial

To begin using DomainMapper, clone the repository to a destination on your local computer

0) git clone [https://github.com/FriedLabJHU/DomainMapper](https://github.com/FriedLabJHU/DomainMapper)

To use DomainMapper, the full output from an HMMR3 hmmscan is required as an input.  
Here is how to generate this file:

1) Get a file containing protein sequences in .fasta format (for instance, a proteome file from uniprot)

2) Obtain the [ECOD HMM profile database](http://prodata.swmed.edu/ecod/distributions/ecodf.hmm.tar.gz)

3) Place this file on your computer, untar and unzip ```tar -xzvf ecodf.hmm.tar.gz```

At this point, you will need to run hmmscan.  The user needs to have separately installed HMMER3.

##### Installing HMMR3
```
# On MacOS (Intel Only)
brew install hmmer

# On Ubuntu (x86/x64 Processors Only)
apt install hmmer
```
*Currently there are no versions of HMMER for Windows but HMMER3 on WSL could be a good alternative*

*Additionally, there are currently no fully-supported distributions of HMMER on ARM-bases processors*

Visit [http://hmmer.org/documentation.html](http://hmmer.org/documentation.html) for details.

4) ```hmmpress ecodf.hmm```

5) ```hmmscan -o your_hmmscan_output.hmm.out /path/to/ecodf.hmm /path/to/your_fasta_file.fasta```

6) ```dommap -f /path/to/your_hmmscan_output.hmm.out -o /path/to/your_hmmscan_output.mapped.out```

To read the various additional options available, read the help docstring  
```dommap -h```

## Citation

Pending 

## Funding

* NSF Division of Molecular and Cellular Biology (2045844)

* NIH Training Grant - Program in Molecular Biophysics (T32GM135131)
