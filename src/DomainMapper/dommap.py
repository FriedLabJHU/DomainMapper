import re
import sys
import argparse
from datetime import datetime
from Bio.SearchIO import parse
from DomainMapper import dommap_io, dommap_tools
from DomainMapper.dommap_data_structures import *

argparser = argparse.ArgumentParser(description=dommap_io.descriptionText)

argparser.add_argument("-f", type=str, default="NULL", help="Input path to file from \'hmmscan\'")

argparser.add_argument("-o", type=str, default="NULL", help="Output path for mapped domains")

argparser.add_argument("--ecod_domains", default="NULL", type=str, help="Path to ECOD \'Latest Domains\' text file  (default = file is automatically downloaded [165 MB Free Space Required (deleted after parsing)] [2 MB File Saved])")

argparser.add_argument("--intra_gap", "--intra_domain_gap_tolerance", type=int, default=30, help="Optional gap size between HMM sequence and query sequence for non-contiguous alignment within a domain (default = 30)")

argparser.add_argument("--inter_gap", "--inter_domain_gap_tolerance", type=int, default=30, help="Optional gap size between two domains sequences for non-contiguous merging (default = 30)")

argparser.add_argument("--overlap", "--domain_overlap_tolerance", type=int, default=40, help="Optional overlap between HMM domain sequence and fasta aligment in consecutive or split domains  (default = 40)")

argparser.add_argument("--eval_cutoff", type=float, default=1e-5, help="Optional upper bound tolerance of the E-value  (default = 1e-5)")

argparser.add_argument("--update", help="Update ECOD \'Latest Domains\'", default=False, action="store_true")

args = argparser.parse_args()

if len(sys.argv) < 2:
    dommap_io.error_msg("No Arguments Passed. View help page with \'dommap -h\'")

if args.ecod_domains == "NULL":
    if args.update and len(sys.argv) < 3:
        dommap_tools.update()
        dommap_io.notice_msg("Latest domain definitions updated.")
    if args.update and len(sys.argv) > 2:
        dommap_tools.update()
    ecod_domain_dict = dommap_tools.load()
else:
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

## DM 3.0.0 should change to a more robust object
# Final formatted output
output_lines = list()

# Counting total number of proteins from input hmm file
hmmscan = parse(args.f,"hmmer3-text")

num_proteins = sum(1 for _ in hmmscan)

if not num_proteins:
    err_msg = "Input hmmscan file could not be read. View help page with \'dommap -h\'.\n\nOne common reasons for this error is providing a `hmmscan --domtblout` file instead of a `hmmscan -o` file."
    dommap_io.error_msg(err_msg)

# Initializing progess bar
dommap_io.progress_bar(0, num_proteins, prefix = "Mapping:", suffix = "Complete", length = 50)

# This is here to make the creation of Domain objects easier
params = [args.intra_gap, args.inter_gap, args.overlap]

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
            potential_noncontig_domains = DomainMap()
            overlapping_potential_noncontig_domains = DomainMap()

            for hsp in hit.hsps:
                domain = Domain(hsp, *params)

                # Keep domain if the E-value is less than or equal to cutoff
                if domain.e_val <= args.eval_cutoff:
                    potential_noncontig_domains.append(domain)

            # Check if any potential non-contig. domains must be combined
            # By referencing domain_A from [:-1] (all but the last) and domain_B from [a+1:] (from index one more than "A" to the end)
            # We are guaranteed to only check unique pairs of domans against each other
            # In a non-contiguous domain, what happens is that the various AA-ranges correspond to distinct portions in the query sequence and in the HMM model
            
            if len(potential_noncontig_domains) > 1:
                for a,domain_A in enumerate(potential_noncontig_domains[:-1]):
                    for b,domain_B in enumerate(potential_noncontig_domains[a+1:]):

                        # If the ranges of the hmm are indexed first for the downsteam query domain and then the upstream query domain, consider it a circular permutant
                        if domain_A.map_intersection(domain_B) <= args.overlap \
                            and not domain_A.hmm_intersection(domain_B) \
                                and ((domain_A.map_range[0] < domain_B.map_range[0] and domain_A.hmm_range[0] > domain_B.hmm_range[0]) or (domain_A.map_range[0] > domain_B.map_range[0] and domain_A.hmm_range[0] < domain_B.hmm_range[0])) \
                                    and domain_A and domain_B and domain_A.f_group == domain_B.f_group:
                            
                            # Mark as circular permutant
                            domain_B.update_topology(f"CP")

                            # Merge the domain data
                            domain_B.merge(domain_A)

                            # Remove domain_A
                            potential_noncontig_domains[a] = None

                        # If the query ranges do not overlap, 
                        if domain_A.map_intersection(domain_B) <= 0 \
                            and domain_A and domain_B and domain_A.f_group == domain_B.f_group:

                            # If two query domains are considerably close, merge them
                            if domain_A.map_range[-1] < domain_B.map_range[0] and (domain_B.map_range[0] - domain_A.map_range[-1]) < args.inter_gap\
                                and not domain_A.hmm_intersection(domain_B):
                                # Merge the domain data
                                domain_B.merge(domain_A)

                                # Remove domain_A
                                potential_noncontig_domains[a] = None
                            
                            # If two query domains cannot merge, they are non-contigous and could potentially be host an insertional domain with overlap
                            elif domain_A.map_range[-1] < domain_B.map_range[0] and (domain_B.map_range[0] - domain_A.map_range[-1]) >= args.inter_gap\
                                and not domain_A.hmm_intersection(domain_B):
                                # Prevent double saving of domains in `overlapping_potential_noncontig_domains`
                                if domain_A not in overlapping_potential_noncontig_domains:

                                    overlapping_potential_noncontig_domains.append(domain_A)
                                    
                                    potential_noncontig_domains[a] = None
                                # If domain_B was added but more domains follow, delete the last added domain in list `potential_noncontig_domains`
                                elif domain_A in overlapping_potential_noncontig_domains:
                                    
                                    potential_noncontig_domains[a] = None 


                                # Prevent double saving of domains in `overlapping_potential_noncontig_domains`
                                if domain_B not in overlapping_potential_noncontig_domains:

                                    overlapping_potential_noncontig_domains.append(domain_B)

                                # Delete the last added domain in list `potential_noncontig_domains`
                                elif b+a+1 == len(potential_noncontig_domains):
                                    
                                    potential_noncontig_domains[b+a+1] = None 

                # preparing to remove overlapping domains that are non-contiguous
                if len(overlapping_potential_noncontig_domains) > 1:
                    
                    overlapping_potential_noncontig_domains.update_overlap_matrix()
                    overlapping_potential_noncontig_domains.eliminate_overlapping_domains()
                    
                    nonoverlapping_potential_noncontig_domains = DomainMap()
                    for domain in overlapping_potential_noncontig_domains:
                        if domain:
                            nonoverlapping_potential_noncontig_domains.append(domain)
                    
                    if len(nonoverlapping_potential_noncontig_domains) > 1:
                        for a, domain_A in enumerate(nonoverlapping_potential_noncontig_domains[:-1]):
                            domain_B = nonoverlapping_potential_noncontig_domains[a+1]

                            domain_B.merge(domain_A)

                            nonoverlapping_potential_noncontig_domains[a] = None

                        potential_noncontig_domains = DomainMap()
                        potential_noncontig_domains.append(nonoverlapping_potential_noncontig_domains[-1])

                    else:
                        potential_noncontig_domains = nonoverlapping_potential_noncontig_domains

            # Remove any domains with low E-values
            for pot_nc_dom in potential_noncontig_domains:
                if pot_nc_dom and pot_nc_dom.e_val < args.eval_cutoff:
                    potential_domain_mappings.append(pot_nc_dom)

    potential_domain_mappings.update_overlap_matrix()
    potential_domain_mappings.eliminate_overlapping_domains()

    # print(potential_domain_mappings._overlap_matrix.tolist())
    # for i,dom in enumerate(potential_domain_mappings):
    #     if dom:
    #         print(i, potential_domain_mappings._overlap_matrix[i].tolist())

    ## Change too io

    # Final domains
    for pot_dom_map in potential_domain_mappings:
        if pot_dom_map:
            mapped_domains.append(pot_dom_map)
    
    # Label the intervening domains (domains that lie within non-contiguous domains)
    if len(mapped_domains) > 1: # only proteins with multiple domains can contain insertional domains
        for a,domain_A in enumerate(mapped_domains):
            for b,domain_B in enumerate(mapped_domains):

                # Check if domain B resides within domain A
                # First, create the set of all "missing residues" in A
                alignment_gap_A = [x for x in range(domain_A.map_range[0],domain_A.map_range[-1]) if x not in domain_A.map_range] 
                num_res_B_in_aln_gap_A = len(set(alignment_gap_A).intersection(set(domain_B.map_range)))

                # Then, if all of B"s residues (with 15 (default) allowed as exception) lie within a gap of A mark as IS of that domain
                if num_res_B_in_aln_gap_A > domain_B.map_len - args.overlap and domain_B.map_len - args.overlap > 0: 
                    if "IS" not in mapped_domains[b].topology: # only count insertional domains once even if they are found to be insertional to multiple host domains
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
                # need to catch these all before
                if "NC" not in domain.topology:
                    """Fix counter"""
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
