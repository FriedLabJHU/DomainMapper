import re
import sys
import argparse
from datetime import datetime
from Bio.SearchIO import parse, to_dict
from DomainMapper import LatestDomains, IODomains

#TODO
# Add compatibility for multiple Properties (NC, CP, IS)
# Add CP identification in single HSPs
# Print domains in order of the first index in the list

# returns the length of the set A∩B, can also return an array that includes all intersecting residue indices
def RangeIntersection(range_A, range_B, array = False):
    intersection = set(range_A).intersection(range_B)
    length = len(intersection)
    if not array:
        return length
    else:
        return length, intersection

def hsp_range_finder(hsp, query_range_idx, hmm_range_idx):
    query_start, query_end = query_range_idx
    hmm_start, hmm_end = hmm_range_idx

    alignment = hsp.aln
    query_aln = alignment[0]
    hmm_aln = alignment[1]

    gap_ranges = []

    if (query_end - query_start) - (hmm_end - hmm_start) > args.gap:
        alignment_gap = re.finditer('\.{'+str(args.gap)+',}', str(hmm_aln.seq))
        
        # Finding the range of the gap in the hmm alignment
        # Match the hmm gap range start/end index to the query sequence index
        if alignment_gap:
            for gap in alignment_gap:
                gap_start_idx = gap.start()
                gap_end_idx = gap.end()
                gap_start_query_idx = gap_start_idx + query_start - len(re.findall('-',str(query_aln.seq[:gap_start_idx])))
                gap_end_query_idx = gap_end_idx + query_start - len(re.findall('-',str(query_aln.seq[:gap_start_idx])))
                gap_ranges += list(range(gap_start_query_idx, gap_end_query_idx))
                
    # Remove the gap from the range, so that it doesn't hog up those resi's (there may be another domain inside)
    # If no gaps are found the whole will be returned
    domain_map_range = [x for x in range(query_start, query_end) if x not in gap_ranges]

    return domain_map_range

# mark_for_deletion is a global variable
def eliminate_overlapping_domains(idx, eval, overlap_map):
        if mark_for_deletion[idx]:
            return

        # save all idx and eval for residues with which the reference domain overlaps with
        tmp_idx = [idx]
        tmp_eval = [eval]
        for i,overlap in enumerate(overlap_map):
            if overlap and not mark_for_deletion[i]: # skips over i=idx since that always maps to 0
                F_grp_ol, eval_ol = potential_domain_mappings[i][0:2]
                tmp_idx.append(i)
                tmp_eval.append(eval_ol)

        # check if the reference domain is the best fit for it's sequence coverage, if not check if any of the overlaping domains are... 
        min_eval_idx = tmp_idx[tmp_eval.index(min(tmp_eval))]
        if  min_eval_idx == idx:
            for j in tmp_idx:
                if j != idx:
                    mark_for_deletion[j] = 1
            return
        else:
            eliminate_overlapping_domains(min_eval_idx, min(tmp_eval), overlap_logic_array[min_eval_idx])

descriptionText = """"""
argparser = argparse.ArgumentParser(description=descriptionText)

argparser.add_argument('-f', type=str, default='NULL', help='Input path to file from \'hmmscan\'')
argparser.add_argument('-o', type=str, default='NULL', help='Output path for mapped domains')
argparser.add_argument('--ecod_domains', default='NULL', type=str, help='Path to ECOD \'Latest Domains\' text file  (default = file is automatically downloaded [165 MB Free Space Required (deleted after parsing)] [2 MB File Saved])')
argparser.add_argument('--gap', '--domain_gap_tolerance', type=int, default=35, help='Optional gap size between HMM domain sequence and fasta aligment  (default = 35)')
argparser.add_argument('--overlap', '--domain_overlap_tolerance', type=int, default=15, help='Optional overlap between HMM domain sequence and fasta aligment in consecutive or split domains  (default = 15)')
argparser.add_argument('--eval_cutoff', type=float, default=1e-5, help='Optional upper bound tolerance of the E-value  (default = 1e-5)')
argparser.add_argument('--update', help='Update ECOD \'Latest Domains\'', default=False, action="store_true")

args = argparser.parse_args()

if len(sys.argv) < 1:
    IODomains.ErrorMsg("No Arguments Passed. View help page with \'DomainMapper.py -h\'.")

if args.f == 'NULL':
    IODomains.ErrorMsg("No Input hmmscan file provided. View help page with \'DomainMapper.py -h\'")

if args.o == 'NULL':
    IODomains.ErrorMsg("No Output path provided. View help page with \'DomainMapper.py -h\'")


if args.ecod_domains == 'NULL':
    if args.update:
        LatestDomains.update()
    ecod_domain_dict = LatestDomains.load()
else:
    ecod_domain_dict = LatestDomains.load(args.ecod_domains)

if args.gap < 0:
    IODomains.ErrorMsg("Non-positive option detected for gap size tolerance. Please ensure all numerical arguments are positive numbers. View help page with \'DomainMapper.py -h\'")

if args.overlap < 0:
    IODomains.ErrorMsg("Non-positive option detected for overlap tolerance. Please ensure all numerical arguments are positive numbers. View help page with \'DomainMapper.py -h\'")

if args.eval_cutoff < 0:
    IODomains.ErrorMsg("Non-positive option detected for E-value cuttoff. Please ensure all numerical arguments are positive numbers. View help page with \'DomainMapper.py -h\'")


# Non-contiguous domain counter
NC_domain_cnt = 0

# Circular Permutant domain counter
CP_domain_cnt = 0

# Interveneing domain counter
IS_domain_cnt = 0

# Total domain counter
Tot_domain_cnt = 0

output_lines = list()

hmmscan = parse(args.f,'hmmer3-text')
# Very inefficient step but worth it
num_proteins = len(to_dict(hmmscan).keys())
IODomains.printProgressBar(0, num_proteins, prefix = 'Mapping:', suffix = 'Complete', length = 50)

# counting the protiens in the hmmscan generator deletes it for some reason, parsing twice is not nice but not that inefficient
hmmscan = parse(args.f,'hmmer3-text')
for p_idx, protein in enumerate(hmmscan):
    accession = protein.id
    potential_domain_mappings = list()
    for hit in protein.hits:
        # Single sequence alignment hit
        if len(hit.hsps) == 1:
            query_range = hit.hsps[0].query_range
            hmm_range = hit.hsps[0].hit_range
            evalue = hit.hsps[0].evalue_cond
            F_group = hit.hsps[0].hit_id
            query_property = list()

            # Keep if the E-value is lower than the cutoff
            if evalue < args.eval_cutoff:
                map_range = hsp_range_finder(hit.hsps[0], query_range, hmm_range)
                potential_domain_mappings.append([F_group, evalue, map_range, hmm_range, query_property])
        
        # Multiple sequence alignment hits
        # There is potential for these to be non-contiguous, circular permutants, or insertional domains
        if len(hit.hsps) > 1:
            potential_noncontig_domains = list()
            for hsp in hit.hsps:
                query_range = hsp.query_range
                hmm_range = hsp.hit_range
                evalue = hsp.evalue_cond
                F_group = hsp.hit_id
                query_property = list()

                # Keep if the E-value is lower than the cutoff
                if evalue < args.eval_cutoff:
                    map_range = hsp_range_finder(hit.hsps[0], query_range, hmm_range)
                    potential_noncontig_domains.append([F_group, evalue, map_range, hmm_range, query_property])

            # Check if any potential non-contig. domains must be combined
            # By referencing domain_A from [:-1] (all but the last) and domain_B from [a+1:] (from index one more than "A" to the end)
            # We are guaranteed to only check unique pairs of domans against each other
            # In a non-contiguous domain, what happens is that the various AA-ranges correspond to distinct portions in the query sequence and in the HMM model
            
            if len(potential_noncontig_domains) > 1:
                for a,domain_A in enumerate(potential_noncontig_domains[:-1]):
                    for b,domain_B in enumerate(potential_noncontig_domains[a+1:]):
                        
                        F_grp_A, eval_A, map_rng_A, hmm_rng_A = domain_A[0:4]
                        F_grp_B, eval_B, map_rng_B, hmm_rng_B = domain_B[0:4]
                        
                        # If the ranges of the hmm are indexed first for the downsteam query sequence and then the upstream query sequence, consider it a circular permutant
                        ### Here test
                        if RangeIntersection(map_rng_A,map_rng_B) <= args.overlap and RangeIntersection(list(range(hmm_rng_A[0],hmm_rng_A[1])),list(range(hmm_rng_B[0],hmm_rng_B[1]))) < args.overlap and ((map_rng_A[0] < map_rng_B[0] and hmm_rng_A[0] > hmm_rng_B[0]) or (map_rng_A[0] > map_rng_B[0] and hmm_rng_A[0] < hmm_rng_B[0])) and F_grp_A != 'null' and F_grp_B != 'null' and F_grp_A == F_grp_B:
                            # Mark as Circ. Permut.
                            potential_noncontig_domains[b+a+1][4].append('CP')
                            
                            # If there is a small space between the circ. permut. domains close it off and consider it a single contiguous domain
                            # otherwise, if the space is larger than the overlap tolerance it can also be marked as non-contig. later on in this program
                            if map_rng_A[-1] < map_rng_B[0] and (map_rng_B[0] - map_rng_A[-1]) < args.overlap:
                                for x in range(map_rng_A[-1], map_rng_B[0]):
                                    if x not in map_rng_A:
                                        map_rng_A.append(x)

                            potential_noncontig_domains[a] = ['null',1,[],[]]
                            potential_noncontig_domains[b+a+1][1] = min(eval_A, eval_B)
                            potential_noncontig_domains[b+a+1][2] = map_rng_A + map_rng_B
                            potential_noncontig_domains[b+a+1][3] = (hmm_rng_A[0], hmm_rng_B[1])

                        if RangeIntersection(map_rng_A,map_rng_B) <= args.gap and RangeIntersection(list(range(hmm_rng_A[0],hmm_rng_A[1])),list(range(hmm_rng_B[0],hmm_rng_B[1]))) < args.overlap and F_grp_A != 'null' and F_grp_B != 'null' and F_grp_A == F_grp_B:
                            # If two domain mappings are considerably close, consider them one domain
                            # otherwise, it is a non-contig. domain which will be marked as such later in this program
                            if map_rng_A[-1] < map_rng_B[0] and (map_rng_B[0] - map_rng_A[-1]) < args.gap:
                                for x in range(map_rng_A[-1], map_rng_B[0]):
                                    if x not in map_rng_A:
                                        map_rng_A.append(x)

                            potential_noncontig_domains[a] = ['null',1,[],[]]
                            potential_noncontig_domains[b+a+1][1] = min(eval_A, eval_B) # maybe find a way to combine these values? Same up top
                            potential_noncontig_domains[b+a+1][2] = map_rng_A + map_rng_B
                            potential_noncontig_domains[b+a+1][3] = hmm_rng_A + hmm_rng_B

            # Remove any domains with low E-values
            for pot_nc_dom in potential_noncontig_domains:
                F_group, eval = pot_nc_dom[0:2]

                if F_group != 'null' and eval < args.eval_cutoff:
                    potential_domain_mappings.append(pot_nc_dom)

    # Mark all overlapping domains
    mapped_domains = []
    overlap_logic_array = [[0 for i in range(len(potential_domain_mappings))] for j in range(len(potential_domain_mappings))]
    for a,domain_A in enumerate(potential_domain_mappings[:-1]):
        for b,domain_B in enumerate(potential_domain_mappings[a+1:]):
            F_grp_A, eval_A, map_rng_A, hmm_rng_A = domain_A[0:4]
            F_grp_B, eval_B, map_rng_B, hmm_rng_B = domain_B[0:4]
            
            # If there's more than 15 (default) residues of overlap, there might be a conflict
            if RangeIntersection(map_rng_A,map_rng_B) > args.overlap:
                
                # Here we check for situtations were a domain might overlap by <=15 residues on each side of other domains, allowing for "feathering" of domains.
                # This could be bugged for any domains under twice the overlap length
                if RangeIntersection(map_rng_A,map_rng_B) <= 2*args.overlap and len(map_rng_A) >= 2*args.overlap:
                    mid_rng_idx_B = len(map_rng_B)//2
                    # More than 15 residue overlaps on either side of a domain and it wil be marked overlapping
                    if RangeIntersection(map_rng_A,map_rng_B[:mid_rng_idx_B]) >= args.overlap or RangeIntersection(map_rng_A,map_rng_B[mid_rng_idx_B:]) >= args.overlap:
                        overlap_logic_array[a][b+a+1] = 1
                        overlap_logic_array[b+a+1][a] = 1

                # More than twice the overlap tolerance and it will be marked overlapping
                else:
                    overlap_logic_array[a][b+a+1] = 1
                    overlap_logic_array[b+a+1][a] = 1
            
            # Smaller domains (< gap size) must be treated differently since their overlap could be 100% but smaller than our tolerances
            if len(map_rng_A) < args.gap:
                if float(RangeIntersection(map_rng_A,map_rng_B))/float(len(map_rng_A)) > 0.7:
                    overlap_logic_array[a][b+a+1] = 1
                    overlap_logic_array[b+a+1][a] = 1

    # Recursively check all overlap maps and mark those with poor overlap and high E-vlaues for deletion
    mark_for_deletion = [0 for i in range(len(potential_domain_mappings))]
    for a,dom_overlap_map in enumerate(overlap_logic_array):
        F_grp_A, eval_A = potential_domain_mappings[a][0:2]
        eliminate_overlapping_domains(a,eval_A,dom_overlap_map)

    # Deleting
    for a,delete in enumerate(mark_for_deletion):
        if delete:
            potential_domain_mappings[a] = ['null',1,[],[]]

    # Final domains
    for pot_dom_map in potential_domain_mappings:
        F_group, eval = pot_dom_map[0:2]
        if F_group != 'null' and eval < args.eval_cutoff:
            # Count the total number of NC and CP domains after filtering
            try:
                dom_prop = pot_dom_map[4]
                if 'CP' in dom_prop:
                    CP_domain_cnt += 1
            except IndexError:
                pass

            mapped_domains.append(pot_dom_map)
    
    # Label the intervening domains (domains that lie within non-contiguous domains)
    for a,domain_A in enumerate(mapped_domains):
        for b,domain_B in enumerate(mapped_domains):
            F_grp_A, eval_A, map_rng_A, hmm_rng_A = domain_A[0:4]
            F_grp_B, eval_B, map_rng_B, hmm_rng_B = domain_B[0:4]

            # Check if domain B resides within domain A
            # First, create the set of all 'missing residues' in A
            alignment_gap_A = [x for x in range(map_rng_A[0],map_rng_A[-1]) if x not in map_rng_A] 
            num_res_B_in_aln_gap_A = RangeIntersection(alignment_gap_A,map_rng_B)

            # Then, if all of B's residues (with 15 (default) allowed as exception) lie within a gap of A mark as IS of that domain
            if num_res_B_in_aln_gap_A > len(map_rng_B) - args.overlap: 
                if str(a) == "0":
                    mapped_domains[b][4].append('IS-'+str(a))
                    IS_domain_cnt += 1
                else:
                    mapped_domains[b][4].append('IS-'+str(a))

    #Now just output this to a file
    for a,domain in enumerate(mapped_domains):
        #reformat the residue range into a nice tidy little string
        F_group, eval, query_range, hmm_range, domain_properties = domain[:]

        start = query_range[0] 
        residue_range_as_string = str(start+1)  #the +1 is an artifact caused by BioPython re-indexing 1 to 0.
        for i in range(1,len(query_range)):
            if query_range[i] == query_range[i-1] + 1:
                pass #they're consecutive
            else:
                residue_range_as_string += ('-{},{}').format(str(query_range[i-1]+1), str(query_range[i]+1)) #Again, +1's is because of BioPython indexing
                # need to catch these all before
                if "NC" not in mapped_domains[a][4]:
                    mapped_domains[a][4].append("NC")
                    NC_domain_cnt += 1
        residue_range_as_string += ('-{}').format(str(query_range[-1]+1))
        
        #try to find the domain in the domain dict else output the F group from the hmmscan
        if F_group in ecod_domain_dict.keys():
            f_id = ecod_domain_dict[F_group][0]
            arch = ecod_domain_dict[F_group][1]
            x_group = ecod_domain_dict[F_group][2]
            t_group = ecod_domain_dict[F_group][3]
            output_lines.append(('{}\t'*9).format(accession, str(eval), residue_range_as_string, ' '.join(str(s) for s in domain_properties), arch, x_group, t_group, F_group, f_id)+'\n')
        else:
            output_lines.append(('{}\t'*9).format(accession, str(eval), residue_range_as_string, ' '.join(str(s) for s in domain_properties), "N/A", "N/A", "N/A", F_group, "N/A")+'\n')

        Tot_domain_cnt += 1 # yes I know taking the len of mapped_domains is fast but this shows future readers of this code where we are taking final count

    IODomains.printProgressBar(p_idx + 1, num_proteins, prefix = 'Mapping:', suffix = 'Complete', length = 50)

with open(args.o, 'w') as mapped_domains_file:
    mapped_domains_file.write(
        IODomains.FileHeader(datetime.now(), args.f, args.o, args.gap, args.overlap, args.eval_cutoff, Tot_domain_cnt, NC_domain_cnt, CP_domain_cnt, IS_domain_cnt))

    for line in output_lines:
        mapped_domains_file.write(line)