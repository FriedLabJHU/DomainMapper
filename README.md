# DomainMapper
DomainMapper is an HMM parser designed to annotate proteins' domain structures from sequence alone using the ECOD database.  
DomainMapper annotates non-contiguous, insertional, and circularly permuted domains as well.

## Dependencies

In order to use DomainMapper you will need ```wget``` and  ```BioPython```
These can be installed using ```pip```

##### Installing Dependencies
```
# installing wget
pip install wget

# installing BioPython
pip install Bio
```

## Running DomainMapper

``` python DomainMapper.py -o raw_hmmscan_output.hmm.out -f mapped_protein_domains.mapped.out```

## Documentation

  -h, --help            show this help message and exit
  -f F                  Input path to file from 'hmmscan'
  -o O                  Output path for mapped domains
  --ecod_domains ECOD_DOMAINS
                        Path to ECOD 'Latest Domains' text file (default = file is automatically downloaded [165 MB
                        Free Space Required (deleted after parsing)] [2 MB File Saved])
  --gap GAP, --domain_gap_tolerance GAP
                        Optional gap size between HMM domain sequence and fasta aligment (default = 40)
  --overlap OVERLAP, --domain_overlap_tolerance OVERLAP
                        Optional overlap between HMM domain sequence and fasta aligment in consecutive or split
                        domains (default = 15)
  --eval_cutoff EVAL_CUTOFF
                        Optional upper bound tolerance of the E-value (default = 1e-5)
  -v                    Verbosity

## HMMR3 & DomainMapper Tutorial

To begin using DomainMapper, clone the repository to a destination on your local computer

0) git clone [https://github.com/FriedLabJHU/DomainMapper](https://github.com/FriedLabJHU/DomainMapper)

To use DomainMapper, the full output from an HMMR3 hmmscan is required as an input.  
Here is how to generate this file:

1) Get a file containing protein sequences in .fasta format (for instance, a proteome file from uniprot), save this in ```./DomainMapper```

2) Obtain the [ECOD HMM profile database](http://prodata.swmed.edu/ecod/distributions/ecodf.hmm.tar.gz)

3) Place this file in the DomainMapper directory, untar and unzip ```tar -xzvf ecodf.hmm.tar.gz```

At this point, you will need to run an hmmscan.  The user needs to have separately installed HMMER3.

##### Installing HMMR3
```
# On MacOS
brew install hmmer
# n Ubuntu
apt install hmmer
```

Visit [http://hmmer.org/documentation.html](http://hmmer.org/documentation.html) for details.

4) ```hmmpress ecodf.hmm```

5) ```hmmscan -o your_hmmscan_output.hmm.out ecodf.hmm your_fasta_file.fasta```

6) ```python src/DomainMapper.py -f your_hmmscan_output.hmm.out -o your_hmmscan_output.mapped.out```

To read the various additional options available, read the help docstring
```python DomainMapper.py -h```

## Citation

Coming Soon!

## Funding

* NSF Division of Molecular and Cellular Biology (2045844)
* NIH Training Grant - Program in Molecular Biophysics (T32GM008403)