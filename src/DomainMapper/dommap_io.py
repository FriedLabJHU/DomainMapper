import sys

def file_header(time, in_file, out_file, intra_gap_tol, inter_gap_tol, overlap_tol, eval_tol, Tot_prot_cnt, Tot_cnt, NC_cnt, CP_cnt, IS_cnt):
    fileHeader = """#===========================================================================================
#  DOMAIN MAPPER v3.0.2
#  Johns Hopkins Univeristy - September 22nd, 2022
#  Edgar Manriquez-Sandoval, M.S. - Dept. of Biophysics
#  emanriq1@jhu.edu
#  & 
#  Stephen D. Fried, Ph.D. - Dept. of Chemistry
#  sdfried@jhu.edu
#===========================================================================================
#  Excecuted on:
#               {}
#  Input HMM: 
#               {}
#  Output:
#               {}
#  Options:
#               Intra domain gap = {:2d}
#               Inter domain gap = {:2d}
#               overlap = {:2d}
#               E-value cutoff = {:1.2e}
#  Domain Counts:
#               Total Proteins: {:6d}         Total Domains:  {:6d}
#                                                        NC : {:3d} ({:.2%})
#                                                        CP : {:3d} ({:.2%})
#                                                        IS : {:3d} ({:.2%})
#  Property Definitions:
#               CP = Circular Permutant Domain
#               NC = Non-Contiguous Domain
#               IS = InSertional Domain
#===========================================================================================
# Accession\tE-Value\tResidue Range\tProperty\tArchitecture\tX-group\tT-group\tF-group\tF-id
""".format(time, in_file, out_file, intra_gap_tol, inter_gap_tol, overlap_tol, eval_tol, Tot_prot_cnt, Tot_cnt, NC_cnt, float(NC_cnt)/float(Tot_cnt), CP_cnt, float(CP_cnt)/float(Tot_cnt), IS_cnt, float(IS_cnt)/float(Tot_cnt),)
    return fileHeader

# This was stolen from: Greenstick @ https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console?page=1&tab=votes#tab-top
# Headless and fast
def progress_bar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

def error_msg(msg):
    ErrMsg = 'ERROR: ' + msg + '\n' + 'System Exiting...\n'
    sys.stderr.write(ErrMsg)
    sys.stderr.flush()
    return sys.exit()

def warning_msg(msg):
    WrnMsg = 'WARNING: ' + msg + '\n'
    sys.stderr.write(WrnMsg)
    sys.stderr.flush()
    return None

def notice_msg(msg):
    NtcMsg = 'NOTICE: ' + msg + '\n'
    sys.stderr.write(NtcMsg)
    sys.stderr.flush()
    return None

descriptionText = \
"""
DomainMapper is a HMMER3 output parser designed to annotate protein domains using the ECOD database.\n
\n
`https://github.com/FriedLabJHU/DomainMapper`

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
"""