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
usage: dommap [-h] [-f F] [-o O] [--ecod_domains ECOD_DOMAINS] [--intra_gap INTRA_GAP] [--inter_gap INTER_GAP] [--overlap OVERLAP] [--frac_overlap FRAC_OVERLAP] [--eval_cutoff EVAL_CUTOFF] [--update]

arguments:
  -h, --help            show this help message and exit
  -f F                  Input path to file from 'hmmscan'
  -o O                  Output path for mapped domains
  --dom_def DOM_DEF     Path to ECOD 'Latest Domains' text file (default = file is automatically downloaded [165 MB Free Space Required (deleted
                        after parsing)] [2 MB File Saved])
  --intra_gap INTRA_GAP, --intra_domain_gap_tolerance INTRA_GAP
                        Optional minimum gap size within a high-scoring pair for those residues to be carved out, generating a non-contiguous hit
                        (default = 30)
  --inter_gap INTER_GAP, --inter_domain_gap_tolerance INTER_GAP
                        Optional minimum gap size between two high-scoring pairs for the residues inbetween to be left out, generating a non-
                        contiguous hit (default = 30)
  --overlap OVERLAP, --domain_overlap_tolerance OVERLAP
                        Optional overlap between high-scoring pairs to mandate an elimination (default = 40)
  --frac_overlap FRAC_OVERLAP, --fractional_domain_overlap_tolerance FRAC_OVERLAP
                        Optional fractional overlap between high-scoring pairs to mandate an elimination (0.0 - 1.0) (default = 0.7)
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

## Reading DomainMapper Output

The first 31 lines consist of a header containing input and output file names, and parameter values. This is followed by summary statistics for the domains identified in the FASTA file submitted.  

Specifically:

* Total Domains - Self explanatory
* NC - Number of non-contiguous domains
* IS - Number of insertional domains
* CP - Number of circularly permutated domains

The second part of the output contains the list of identified domains. Each of which is assigned values in 9 columns (for ECOD searches).
* Column 1 - Protein indentifier from FASTA
* Column 2 - Conditional E-value
* Column 3 - Domain Range
* Column 4 - Domain Property (NC, IS, CP)
* Column 5 - ECOD Architecture (Will be NA if ECOD HMMs are not used)
* Column 6 - ECOD X-Group (Will be NA if ECOD HMMs are not used)
* Column 7 - ECOD T-Group (Will be NA if ECOD HMMs is not used)
* Column 8 - ECOD F-Group
* Column 9 - ECOD F-ID (Will be NA if ECOD HMMs is not used)

Domains with complex topologie (i.e. nesting, weaving, etc.) can be identified by the "NC IS" property flags in column 4.

## Citation

Manriquez-Sandoval, E, Fried, SD. DomainMapper: Accurate domain structure annotation including those with non-contiguous topologies. Protein Science. 2022; 31( 11):e4465. https://doi.org/10.1002/pro.4465

## Funding

* NSF Division of Molecular and Cellular Biology (2045844)

* NIH Training Grant - Program in Molecular Biophysics (T32GM135131)
