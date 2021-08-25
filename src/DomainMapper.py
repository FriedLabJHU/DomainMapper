import os
import re
import sys
import argparse
import LatestDomains
from datetime import datetime
from Bio.SearchIO import parse

# TODO
# Bioinformatics
# Make README / Documentations
# sh script / example script
# Finalize fig 2 
# Include refs in papers -- bioinfo
# Include CP portion
# Implementation of DM of how it is available and how it operates
#   - Needs hmmer and assumes certain level of bioinfo knowl of user


def ErrorMsg(msg):
    ErrMsg = 'ERROR: ' + msg + '\n' + 'System Exiting...'
    print(ErrMsg)

descriptionText="""
DomainMapper is an HMM parser designed to annotate proteins' domain structures from sequence alone using the ECOD database.  
DomainMapper annotates non-contiguous, insertional, and circularly permuted domains as well.

To use DomainMapper, the full output from an HMMR3 hmmscan is required as an input.  
Here is how to generate this file:
1) Get a file containing protein sequences in .fasta format (for instance, a proteome file from uniprot), save this in DomainMapper
2) Obtain the ECOD HMM profile database at http://prodata.swmed.edu/ecod/distributions/ecodf.hmm.tar.gz
3) Place this file in the same directory as DomainMapper, and unzip (tar -xzvf ecodf.hmm.tar.gz)

At this point, you will need to run an hmmscan.  The user needs to have separately installed HMMR3.
HMMR3 is pretty easy to install. On Mac, brew install hmmer. On Ubuntu, apt install hmmer. 
Visit http://hmmer.org/documentation.html for details.

4) hmmpress ecodf.hmm
5) hmmscan -o DomainMapper/your_hmmscan_output.hmm.out ecodf.hmm DomainMapper/your_fasta_file.fasta

6) python DomainMapper/src/DomainMapper.py -f your_hmmscan_output.hmm.out -o your_hmmscan_output.mapped.out

"""
argparser = argparse.ArgumentParser(description=descriptionText)

argparser.add_argument('-f', type=str, default='NULL', help='Input path to file from \'hmmscan\'')
argparser.add_argument('-o', type=str, default='NULL', help='Output path for mapped domains')
argparser.add_argument('--ecod_domains', default='NULL', type=str, help='Path to ECOD \'Latest Domains\' text file  (default = file is automatically downloaded [165 MB Free Space Required (deleted after parsing)] [2 MB File Saved])')
argparser.add_argument('--gap', '--domain_gap_tolerance', type=int, default=40, help='Optional gap size between HMM domain sequence and fasta aligment  (default = 40)')
argparser.add_argument('--overlap', '--domain_overlap_tolerance', type=int, default=15, help='Optional overlap between HMM domain sequence and fasta aligment in consecutive or split domains  (default = 15)')
argparser.add_argument('--eval_cutoff', type=float, default=1e-5, help='Optional upper bound tolerance of the E-value  (default = 1e-5)')
argparser.add_argument('-v', help='Verbosity', action="store_true")

args = argparser.parse_args()

if len(sys.argv) < 1:
    ErrorMsg("No Arguments Passed. View help page with \'DomainMapper.py -h\'.")
    quit()
if args.f == 'NULL':
    ErrorMsg("No Input hmmscan file provided. View help page with \'DomainMapper.py -h\'")
    quit()
if args.o == 'NULL':
    ErrorMsg("No Output path provided. View help page with \'DomainMapper.py -h\'")
    quit()

if args.ecod_domains == 'NULL':
    LatestDomains.fetch()
    ecod_domain_dict = LatestDomains.load()
else:
    if os.path.isfile(args.ecod_domains):
        ecod_input_check = re.search('ecod.latest.domains.txt', args.ecod_domains)
        if ecod_input_check.group() == 'ecod.latest.domains.txt':
            print("Parsing user input: \'{}\'".format(args.ecod_domains))
            ecod_domain_dict = LatestDomains.parse_user_ecod_domain_txt(args.ecod_domains)
        else:
            ErrorMsg("Input \'{}\' does not satify \'ecod.latest.domains.txt\' naming requirement. Change file name to continue. View help page with \'DomainMapper.py -h\'".format(args.ecod_domains))
            quit()
    else:
        ErrorMsg("Input \'{}\' cannot be found. View help page with \'DomainMapper.py -h\'".format(args.ecod_domains))
        quit()

# maybe worth refining this
if args.gap < 0 or args.overlap < 0 or args.eval_cutoff < 0:
    ErrorMsg("Non-positive option detected. Please ensure all numerical arguments are positive numbers. View help page with \'DomainMapper.py -h\'")
    quit()

fileShebang = """#==========================================================================
#  DOMAIN MAPPER v0.0.1
#  Johns Hopkins Univeristy - August 5th, 2021
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
#               gap = {}
#               overlap = {}
#               E-value cutoff = {}
#  Property Definitions:
#               NC = Non-contiguous
#               CP = Circular Permutation
#               IV = Intervening domain followed by a number of the domain it is between (starting from 0)
#==========================================================================
# Accession\tE-value\tReside Range\tPropery\tArchitecture\tX-group\tT-group\tF-group
""".format(datetime.now(), args.f, args.o, args.gap, args.overlap, args.eval_cutoff)

# Non-contiguous domain counter
NC_domain_cnt = 0

# Circular Permutant domain counter
CP_domain_cnt = 0

# Interveneing domain counter
IV_domain_cnt = 0

# Total domain counter
Tot_domain_cnt = 0

hmmscan = parse(args.f,'hmmer3-text')
with open(args.o, 'w') as mapped_domains_file:
    mapped_domains_file.write(fileShebang)

    for p_idx, protein in enumerate(hmmscan):
        accession = protein.id
        potential_domain_mappings = list()
        for hit in protein.hits:
            # Single sequence alignment hit
            if len(hit.hsps) == 1:
                query_range = list(range(hit.hsps[0].query_start, hit.hsps[0].query_end))
                hmm_range = list(range(hit.hsps[0].hit_start, hit.hsps[0].hit_end))
                evalue = hit.hsps[0].evalue_cond
                F_group = hit.hsps[0].hit_id

                # Checking for gap in the alignment in the hmm vs the query (particular protein sequence)
                if abs(len(query_range) - len(hmm_range)) > args.gap:
                    alignment = hit.hsps[0].aln
                    alignment_gap = re.search('\.{'+str(args.gap)+',}', str(alignment[1].seq))
                    
                    # Finding the range of the gap in the alignment
                    if alignment_gap:
                        gap_start_aln_index = alignment_gap.start()
                        gap_end_aln_index = alignment_gap.end()
                        # we should make this comment comprehensible
                        # this next line looks confusing,
                        # but all it is saying is that to go from column #  (what is column number ??) in alignment to actual residue # in query sequence, 
                        # take column number PLUS resi where alignment starts on query MINUS the number of gaps in the query up to that point
                        gap_start_query_index = gap_start_aln_index + query_range[0] - len(re.findall('-',str(alignment[0].seq[:gap_start_aln_index])))
                        gap_end_query_index = gap_end_aln_index + query_range[0] - len(re.findall('-',str(alignment[0].seq[:gap_start_aln_index])))
                        gap_range = list(range(gap_start_query_index, gap_end_query_index))
                        query_range = [x for x in query_range if x not in gap_range]
                
                if evalue < args.eval_cutoff: #only bother keeping it if it is a decent match
                    potential_domain_mappings.append([F_group, evalue, query_range, hmm_range, []])
            
            # Multiple sequence alignment hits
            # There is potential for these to be non-contiguous, circular permutants, or interveneing domains
            if len(hit.hsps) > 1:
                potential_noncontig_domains = list()
                for hsp in hit.hsps:
                    query_range = list(range(hsp.query_start, hsp.query_end))
                    hmm_range = list(range(hsp.hit_start, hsp.hit_end))
                    evalue = hsp.evalue_cond
                    F_group = hsp.hit_id

                    if abs(len(query_range) - len(hmm_range)) > args.gap: #check if there's a big gap in the alignment
                        alignment = hsp.aln #grab the alignment from the hmm out file
                        alignment_gap = re.search('\.{'+str(args.gap)+',}', str(alignment[1].seq)) #find a gap of 40 or more
                        
                        if alignment_gap:
                            gap_start_aln_index = alignment_gap.start()
                            gap_end_aln_index = alignment_gap.end()
                            #this next line looks confusing, but all it is saying is that to go from column # in alignment to actual residue # in query sequence, take column number PLUS resi where alignment starts on query MINUS the number of gaps in the query up to that point
                            gap_start_query_index = gap_start_aln_index + query_range[0] - len(re.findall('-',str(alignment[0].seq[:gap_start_aln_index])))
                            gap_end_query_index = gap_end_aln_index + query_range[0] - len(re.findall('-',str(alignment[0].seq[:gap_start_aln_index])))
                            gap_range = list(range(gap_start_query_index, gap_end_query_index))
                            query_range = [x for x in query_range if x not in gap_range] #remove the gap from the range, so that it doesn't hog up those resis (there may be another domain inside)
                    
                    if evalue < 1e-2: #this is a relaxed cut-off, since it could benefit by 'joining' to  more confident fragment, but let's filter out stuff that is so crappy, it shouldn't even be considered further
                        potential_noncontig_domains.append([F_group, evalue, query_range, hmm_range, []])
                
                #now we need to see if some of these domains need to be combined, into say, a split-domain
                for a,domain_A in enumerate(potential_noncontig_domains):
                    for b,domain_B in enumerate(potential_noncontig_domains):
                        if b > a: #only do unique pairwise comparisons
                            #In a split domain, what happens is that the various AA-ranges correspond to distinct portions in the query sequence and in the HMM model
                            #If this is a split domain, combine them 
                            F_grp_A, eval_A, query_rng_A, hmm_rng_A = domain_A[0:4]
                            F_grp_B, eval_B, query_rng_B, hmm_rng_B = domain_B[0:4]

                            if len(set(query_rng_A).intersection(query_rng_B)) < 10 and len(set(hmm_rng_A).intersection(hmm_rng_B)) < 10 and F_grp_A != 'null' and F_grp_B != 'null':
                                if (query_rng_A[0] < query_rng_B[0] and hmm_rng_A[0] < hmm_rng_B[0]) or (query_rng_A[0] > query_rng_B[0] and hmm_rng_A[0] > hmm_rng_B[0]):
                                    potential_noncontig_domains[b][4].append('NC')
                                    NC_domain_cnt += 1
                                    pass #consistent order suggests it's a normal domain
                                elif (query_rng_A[0] < query_rng_B[0] and hmm_rng_A[0] > hmm_rng_B[0]) or (query_rng_A[0] > query_rng_B[0] and hmm_rng_A[0] < hmm_rng_B[0]):
                                    CP_domain_cnt += 1
                                    potential_noncontig_domains[b][4].append('CP')
                                
                                potential_noncontig_domains[a] = ['null',1,[],[]]
                                potential_noncontig_domains[b][1] = min(eval_A, eval_B)
                                potential_noncontig_domains[b][2] = query_rng_A + query_rng_B
                                potential_noncontig_domains[b][3] = hmm_rng_A + hmm_rng_B

                for pot_nc_dom in potential_noncontig_domains:
                    F_group, eval = pot_nc_dom[0:2]
                    if F_group != 'null' and eval < args.eval_cutoff:
                        potential_domain_mappings.append(pot_nc_dom)

        mapped_domains = []
        for a,domain_A in enumerate(potential_domain_mappings):
            for b,domain_B in enumerate(potential_domain_mappings):
                if b > a:
                    F_grp_A, eval_A, query_rng_A, hmm_rng_A = domain_A[0:4]
                    F_grp_B, eval_B, query_rng_B, hmm_rng_B = domain_B[0:4]
                    
                    if len(set(query_rng_A).intersection(query_rng_B)) > args.overlap: #if there's more than 15 (default) residues of overlap, there's a conflict
                        if eval_A < eval_B: #if A is the 'better' hit for that stretch of residues
                            potential_domain_mappings[b] = ['null',1,[],[]] #cancel B. A stays
                        else: #if B is the 'better' hit for that stretch of residues
                            potential_domain_mappings[a] = ['null',1,[],[]] #cancel A. B stays.

        for pot_dom_map in potential_domain_mappings:
            F_group, eval = pot_dom_map[0:2]
            if F_group != 'null' and eval < args.eval_cutoff:
                mapped_domains.append(pot_dom_map)
        
        #label the intervening domains (domains that lie within split domains)
        for a,domain_A in enumerate(mapped_domains):
            for b,domain_B in enumerate(mapped_domains):
                #check if domain b resides within domain a.
                F_grp_A, eval_A, query_rng_A, hmm_rng_A = domain_A[0:4]
                F_grp_B, eval_B, query_rng_B, hmm_rng_B = domain_B[0:4]

                alignment_gap_A = set([x for x in range(query_rng_A[0],query_rng_A[-1]) if x not in query_rng_A]) #for domain a, create the set of all 'missing residues'
                num_res_B_in_aln_gap_A = len(alignment_gap_A.intersection(query_rng_B))
                if num_res_B_in_aln_gap_A > len(query_rng_B) - args.overlap: #if all of B's residues (with 15 (default) allowed as exception) lie within a gap of A
                    mapped_domains[b][4].append(' IV-'+str(a))
                    IV_domain_cnt += 1


        #Now just output this to a file
        for a,domain in enumerate(mapped_domains):
            #reformat the residue range into a nice tidy little string
            F_group, eval, query_range, hmm_range, domain_properties = domain[0:5]

            start = query_range[0] 
            residue_range_as_string = str(start+1)  #the +1 is an artifact caused by BioPython re-indexing 1 to 0.
            for i in range(1,len(query_range)):
                if query_range[i] == query_range[i-1] + 1:
                    pass #they're consecutive
                else:
                    residue_range_as_string += ('-{},{}').format(str(query_range[i-1]+1), str(query_range[i]+1)) #Again, +1's is because of BioPython indexing
            residue_range_as_string += ('-{}').format(str(query_range[-1]+1))
            
            #try to find the domain in the domain dict else output the F group from the hmmscan
            if F_group in ecod_domain_dict.keys():
                f_id = ecod_domain_dict[F_group][0]
                arch = ecod_domain_dict[F_group][1]
                x_group = ecod_domain_dict[F_group][2]
                t_group = ecod_domain_dict[F_group][3]
                output_line = ('{}\t'*9).format(accession, str(eval), residue_range_as_string, ''.join(str(s) for s in domain_properties), arch, x_group, t_group, F_group, f_id) + '\n'
            else:
                output_line = ('{}\t'*9).format(accession, str(eval), residue_range_as_string, ''.join(str(s) for s in domain_properties), "N/A", "N/A", "N/A", F_group, "N/A") + '\n'
            
            Tot_domain_cnt += 1
            mapped_domains_file.write(output_line)
            if args.v:
                print(output_line)
