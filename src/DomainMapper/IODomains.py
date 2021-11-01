import sys

def FileHeader(time, in_file, out_file, gap_tol, overlap_tol, eval_tol, Tot_cnt, NC_cnt, CP_cnt, IS_cnt):
    fileHeader = \
"""#==========================================================================
#  DOMAIN MAPPER v1.3.0
#  Johns Hopkins Univeristy - November 1st, 2021
#  Edgar Manriquez-Sandoval, B.S. - Dept. of Biophysics
#  emanriq1@jhu.edu
#  & 
#  Stephen D. Fried, Ph.D. - Dept. of Chemistry
#  sdfried@jhu.edu
#==========================================================================
#  Excecuted on:
#               {}
#  Input HMM: 
#               {}
#  Output:
#               {}
#  Options:
#               gap = {:2d}
#               overlap = {:2d}
#               E-value cutoff = {:1.2e}
#  Domain Counts:
#               Total : {:6d}
#               NC : {:3d} ({:.2%})
#               CP : {:3d} ({:.2%})
#               IS : {:3d} ({:.2%})
#  Property Definitions:
#               NC = Non-Contiguous Domain
#               CP = Circular Permutant Domain
#               IS = InSertional Domain
#==========================================================================
# Accession\tE-value\tReside Range\tPropery\tArchitecture\tX-group\tT-group\tF-group\tF-id
""".format(time, in_file, out_file, gap_tol, overlap_tol, eval_tol, Tot_cnt, NC_cnt, float(NC_cnt)/float(Tot_cnt), CP_cnt, float(CP_cnt)/float(Tot_cnt), IS_cnt, float(IS_cnt)/float(Tot_cnt),)
    return fileHeader

def ErrorMsg(msg):
    ErrMsg = 'ERROR: ' + msg + '\n' + 'System Exiting...\n'
    sys.stderr.write(ErrMsg)
    sys.stderr.flush()
    return sys.exit()

def WarningMsg(msg):
    WrnMsg = 'WARNING: ' + msg + '\n'
    sys.stderr.write(WrnMsg)
    return sys.stderr.flush()