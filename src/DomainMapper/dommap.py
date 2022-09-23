import re
import sys
import argparse
from datetime import datetime
from Bio.SearchIO import parse
from DomainMapper import dommap_io, dommap_tools
from DomainMapper.dommap_data_structures import *

# Parsing Arguments
argparser = argparse.ArgumentParser(description=dommap_io.descriptionText)

argparser.add_argument("-f", type=str, default="NULL", help="Input path to file from \'hmmscan\'")

argparser.add_argument("-o", type=str, default="NULL", help="Output path for mapped domains")

argparser.add_argument("--dom_def", default="NULL", type=str, help="Path to ECOD \'Latest Domains\' text file  (default = file is automatically downloaded [165 MB Free Space Required (deleted after parsing)] [2 MB File Saved])")

argparser.add_argument("--intra_gap", "--intra_domain_gap_tolerance", type=int, default=30, help="Optional minimum gap size within a high-scoring pair for those residues to be carved out, generating a non-contiguous hit (default = 30)")

argparser.add_argument("--inter_gap", "--inter_domain_gap_tolerance", type=int, default=30, help="Optional minimum gap size between two high-scoring pairs for the residues inbetween to be left out, generating a non-contiguous hit (default = 30)")

argparser.add_argument("--overlap", "--domain_overlap_tolerance", type=int, default=40, help="Optional overlap between high-scoring pairs to mandate an elimination  (default = 40)")

argparser.add_argument("--frac_overlap", "--fractional_domain_overlap_tolerance", type=float, default=0.7, help="Optional fractional overlap between high-scoring pairs to mandate an elimination (0.0 - 1.0) (default = 0.7)")

argparser.add_argument("--eval_cutoff", type=float, default=1e-5, help="Optional upper bound tolerance of the E-value (default = 1e-5)")

argparser.add_argument("--update", help="Update ECOD \'Latest Domains\'", default=False, action="store_true")

args = argparser.parse_args()

# Checking if the minimum number of arguments required have been passed
if len(sys.argv) < 2:
    dommap_io.error_msg("No Arguments Passed. View help page with \'dommap -h\'")

# Checking which domain definitions to use
if args.dom_def == "NULL":

    #Dommap will use built in domain definitions
    if args.update and len(sys.argv) < 3:

        # Update built in domain definitions, recommended if they are older than 2 months
        dommap_tools.update()

        dommap_io.notice_msg("Latest domain definitions updated.")

    if args.update and len(sys.argv) > 2:
        # Update built in domain definitions, recommended if they are older than 2 months
        dommap_tools.update()

    # Read in built in domain definitions
    ecod_domain_dict = dommap_tools.load()

else:

    # User would like to use their own definitions
    ecod_domain_dict = dommap_tools.load(args.ecod_domains)

if args.f == "NULL":
    dommap_io.error_msg("No Input hmmscan file provided. View help page with \'dommap -h\'")

if args.o == "NULL":
    dommap_io.error_msg("No Output path provided. View help page with \'dommap -h\'")

if args.intra_gap < 0:
    dommap_io.error_msg("Non-positive option detected for the intra-gap size tolerance. Please ensure all numerical arguments are positive numbers. View help page with \'dommap -h\'")

if args.inter_gap < 0:
    dommap_io.error_msg("Non-positive option detected for the inter-gap size tolerance. Please ensure all numerical arguments are positive numbers. View help page with \'dommap -h\'")

if args.overlap < 0:
    dommap_io.error_msg("Non-positive option detected for overlap tolerance. Please ensure all numerical arguments are positive numbers. View help page with \'dommap -h\'")

if args.frac_overlap < 0:
    dommap_io.error_msg("Non-positive option detected for fractional overlap tolerance. Please ensure all numerical arguments are positive numbers. View help page with \'dommap -h\'")

if args.eval_cutoff < 0:
    dommap_io.error_msg("Non-positive option detected for E-value cuttoff. Please ensure all numerical arguments are positive numbers. View help page with \'dommap -h\'")


# Non-contiguous domain counter
NC_domain_cnt = 0

# Circular Permutant domain counter
CP_domain_cnt = 0

# Insertional domain counter
IS_domain_cnt = 0

# Total domain counter
Tot_domain_cnt = 0

# Final formatted output
output_lines = list()

# Initial reading of the input hmm file
hmmscan = parse(args.f,"hmmer3-text")

# Counting total number of proteins from input hmm file
num_proteins = sum(1 for _ in hmmscan)

# If proteins were not detected from the input hmm file, then send error message
# Usually, this is because `--domtblout` was used in HMMER3 instead of `-o`
if not num_proteins:

    err_msg = "Input hmmscan file could not be read. View help page with \'dommap -h\'.\n\nOne common reasons for this error is providing a `hmmscan --domtblout` file instead of a `hmmscan -o` file."
    
    dommap_io.error_msg(err_msg)

# Initial progess bar
dommap_io.progress_bar(0, num_proteins, prefix = "Mapping:", suffix = "Complete", length = 50)

# This is here to make the creation of Domain objects easier
params = [args.intra_gap, args.inter_gap, args.overlap, args.frac_overlap]

# Final reading of the input hmm file
hmmscan = parse(args.f,"hmmer3-text")

for p_idx, protein in enumerate(hmmscan):

    accession = protein.id

    potential_domain_mappings = DomainMap()

    mapped_domains = DomainMap()

    for hit in protein.hits:

        # Single high-scoring pair
        if len(hit.hsps) == 1:

            # Save as Domain() object
            domain = Domain(hit.hsps[0], *params)

            # Keep domain if the E-value is less than or equal to cutoff
            if domain.e_val <= args.eval_cutoff:

                potential_domain_mappings.append(domain)
        
        # Multiple high-scoring pairs
        # When an alignment has multiple HSPs, these domains can have complex topologies (e.g. non-contiguouity, circular permutant, or both, or just repetitive domains)
        if len(hit.hsps) > 1:

            # Create a temporary DomainMap() of domain annotations with multiple HSPs
            # This will then be recursively eliminated so that only the lowest E-Values are retained
            multi_hsps_domains = DomainMap()

            for hsp in hit.hsps:
                
                # Save as Domain() object
                domain = Domain(hsp, *params)

                # Keep domain if the E-value is less than or equal to cutoff
                if domain.e_val <= args.eval_cutoff:

                    multi_hsps_domains.append(domain)

            # Eliminate overlapping HSP's
            multi_hsps_domains.update_overlap_matrix()

            multi_hsps_domains.eliminate_overlapping_domains()

            # Check if any potential non-contig. domains must be combined
            # By referencing domain_A from [:-1] (all but the last) and domain_B from [a+1:] (from index one more than "A" to the end)
            # We are guaranteed to only check unique pairs of domans against each other
            # In a non-contiguous domain, what happens is that the various AA-ranges correspond to distinct portions in the query sequence and in the HMM model
            
            if len(multi_hsps_domains) > 1:

                for a,domain_A in enumerate(multi_hsps_domains[:-1]):

                    for b,domain_B in enumerate(multi_hsps_domains[a+1:]):

                        # If the domains have not been eliminated, check if they can be merged
                        if domain_A and domain_B and domain_A.f_group == domain_B.f_group:

                            # Merge domains if their query (map) ranges do not overlap
                            # And if their hmm ranges do not overlap (70% for small domains)
                            if domain_A.map_intersection(domain_B) <= domain_A.overlap \
                                and domain_A.map_intersection(domain_B)/float(domain_A.map_len) < args.frac_overlap and domain_A.map_intersection(domain_B)/float(domain_B.map_len) < args.frac_overlap \
                                    and domain_A.hmm_intersection(domain_B) <= domain_A.overlap \
                                        and domain_A.hmm_intersection(domain_B)/float(domain_A.hmm_len) < args.frac_overlap and domain_A.hmm_intersection(domain_B)/float(domain_B.hmm_len) < args.frac_overlap:                                
                                
                                # Check to see if this is CP
                                if ((domain_A.map_range[0] < domain_B.map_range[0] and domain_A.hmm_range[0] > domain_B.hmm_range[0]) or (domain_A.map_range[0] > domain_B.map_range[0] and domain_A.hmm_range[0] < domain_B.hmm_range[0])) \
                                        and domain_A and domain_B and domain_A.f_group == domain_B.f_group:
                            
                                    domain_B.update_topology(f"CP")
                                
                                domain_B.merge(domain_A)

                                # Remove domain_A
                                multi_hsps_domains[a] = None

                            else:
                                # Otherwise, retain them, as they could be repetative domains
                                pass

            # Remove any domains with low E-values
            for dom in multi_hsps_domains:

                if dom and dom.e_val < args.eval_cutoff:

                    potential_domain_mappings.append(dom)

    # Eliminate overlapping HITs
    potential_domain_mappings.update_overlap_matrix()

    potential_domain_mappings.eliminate_overlapping_domains()

    # Final domains

    for pot_dom_map in potential_domain_mappings:

        if pot_dom_map:

            mapped_domains.append(pot_dom_map)
    
    # Label the insertional domains (domains that lie within non-contiguous domains)
    if len(mapped_domains) > 1: # only proteins with multiple domains can contain insertional domains

        for a,domain_A in enumerate(mapped_domains):

            for b,domain_B in enumerate(mapped_domains):

                # Check if domain B resides within domain A
                # First, create the set of all "missing residues" in A
                alignment_gap_A = [x for x in range(domain_A.map_range[0],domain_A.map_range[-1]) if x not in domain_A.map_range]

                num_res_B_in_aln_gap_A = len(set(alignment_gap_A).intersection(set(domain_B.map_range)))

                # Then, if all of B"s residues (with 15 (default) allowed as exception) lie within a gap of A mark as IS of that domain
                if num_res_B_in_aln_gap_A > (domain_B.map_len - args.overlap) and (domain_B.map_len - args.overlap) > 0: 

                    mapped_domains[b].update_topology(f"IS")

    #Now just output this to a file
    domain_info = dict()

    for a,domain in enumerate(mapped_domains):

        #reformat the residue range into a nice tidy little string
        start = domain.map_range[0]

        residue_range_as_string = str(start+1)  #the +1 is an artifact caused by BioPython re-indexing 1 to 0.

        for i, mr in enumerate(domain.map_range[:-1]):

            if mr + 1  == domain.map_range[i+1]:

                pass #they"re consecutive

            else:
                residue_range_as_string += ("-{},{}").format(str(mr+1), str(domain.map_range[i+1]+1)) #Again, +1"s is because of BioPython indexing

                domain.update_topology(f"NC")

        residue_range_as_string += ("-{}").format(str(domain.map_range[-1]+1))
        
        domain.res_str = residue_range_as_string

        #try to find the domain in the domain dict else output the F group from the hmmscan
        if domain.f_group in ecod_domain_dict.keys():

            domain.f_id = ecod_domain_dict[domain.f_group][0]

            domain.arch = ecod_domain_dict[domain.f_group][1]

            domain.x_group = ecod_domain_dict[domain.f_group][2]

            domain.t_group = ecod_domain_dict[domain.f_group][3]

        Tot_domain_cnt += 1 # yes I know taking the len of mapped_domains is fast but this shows future readers of this code where we are taking final count

    # print domains out in order of the first index that appears for a given annotation
    final_mapped_domains = sorted(mapped_domains, key = lambda dom: int(dom.map_range[0]))
    
    for dom in final_mapped_domains:

        if "CP" in dom.topology:
            CP_domain_cnt += 1
        if "NC" in dom.topology:
            NC_domain_cnt += 1
        if "IS" in dom.topology:
            IS_domain_cnt += 1

        output_lines.append(("{}\t"*9).format(accession, f"{dom.e_val:3.2e}", dom.res_str, " ".join([top for top in dom.topology]), dom.arch, dom.x_group, dom.t_group, dom.f_group, dom.f_id)+"\n")

    dommap_io.progress_bar(p_idx + 1, num_proteins, prefix = "Mapping:", suffix = "Complete", length = 50)

with open(args.o, "w") as mapped_domains_file:

    mapped_domains_file.write(
        dommap_io.file_header(datetime.now(), args.f, args.o, args.intra_gap, args.inter_gap, args.overlap, args.eval_cutoff, num_proteins, Tot_domain_cnt, NC_domain_cnt, CP_domain_cnt, IS_domain_cnt))

    for line in output_lines:
        
        mapped_domains_file.write(line)
