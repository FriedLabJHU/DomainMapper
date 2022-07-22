# dommmap_data_structures.py 
# This file contains all the data structures and supporting functions for these structures used in DomainMapper

from scipy.stats import combine_pvalues

from numpy import zeros

from re import findall, finditer

from Bio.SearchIO._model.hsp import HSP

# Initializing for type annotations
class Domain:
    pass


# Initializing for type annotations
class DomainMap(list):
    pass


class Domain:
    """
    Data structure for individual domains annotated (mapped) from a HMMER3 output
    """

    def __init__(self, hsp: HSP, intra_gap: int, inter_gap: int, overlap: int):

        self.map_range = self.__map_range_finder(hsp, intra_gap)

        self.map_len = len(self.map_range)

        self.hmm_range = hsp.hit_range

        self.hmm_len = self.hmm_range[1] - self.hmm_range[0] + 1

        self.e_val = hsp.evalue_cond

        self.topology = set()

        self.f_group = hsp.hit_id

        self.intra_gap = intra_gap

        self.inter_gap = inter_gap

        self.overlap = overlap

    def map_intersection(self, dom: Domain, start = 0, end = None):
        """
        Returns the number of residues overlapping between two domains
        """

        return len(set(self.map_range).intersection(set(dom.map_range[start:end])))

    def hmm_intersection(self, dom: Domain):
        """
        Returns the number of residues overlapping between two hmm
        """

        start_A, end_A = self.hmm_range[:2]
        start_B, end_B = dom.hmm_range [:2]
        return len(set(range(start_A, end_A)).intersection(set(range(start_B, end_B))))

    def update_topology(self, topo: str):
        """
        Updates the topology types this domain contains
        """

        self.topology.add(topo)
    
    def update_map_range(self, range: list):
        """
        Updates the map range this domain contains
        """

        self.map_range += range
        tmp_set = list(set(self.map_range))
        self.map_range = sorted(tmp_set)

    def merge(self, dom: Domain):
        """
        Merges two domains into one of the same F-group
        """
        if self.f_group == dom.f_group:

            self.__map_range_filler(dom)

            self.__avg_e_val(dom)

            self.hmm_range = self.hmm_range + dom.hmm_range
            
            for topo in dom.topology:
                self.update_topology(topo)

        else:
            raise AttributeError("Error merging domain.")

    def __map_range_finder(self, hsp: HSP, intra_gap: int):
        """
        This function finds the all protein residue indices which were aligned to an ECOD HMM that define a domain.
        If aligment gaps are less than `intra_gap` they will be filled, else gaps greater than or equal to `intra_gap`
        will remain due to the poteinal of non-contiguous domain topology.

        Parmeters
        ------------
        hsp : Bio.SearchIO._model.hsp.HSP
        High-scoring Pair from an HMM alignment to a query (protein) sequence.

        intra_gap : int
        Maximum tolerated alignment gap.

        Returns
        ------------
        map_ranges : list
        A list of residue indicies which aligned to an ECOD HMM
        """

        gap_ranges = list()

        query_start, query_end = hsp.query_range
        hmm_start, hmm_end = hsp.hit_range

        query_aln = hsp.aln[0]
        hmm_aln = hsp.aln[1]

        if (query_end - query_start) - (hmm_end - hmm_start) >= intra_gap:
            alignment_gap = finditer(r'\.{'+str(intra_gap)+',}', str(hmm_aln.seq))
            if alignment_gap:
                for gap in alignment_gap:
                    gap_start = gap.start()
                    gap_end = gap.end()

                    gap_start_query_idx = gap_start + query_start - len(findall('-',str(query_aln.seq[:gap_start])))
                    gap_end_query_idx = gap_end + query_start - len(findall('-',str(query_aln.seq[:gap_start])))

                    gap_ranges += list(range(gap_start_query_idx, gap_end_query_idx))
                    
        map_ranges = sorted([x for x in range(query_start, query_end) if x not in gap_ranges])
        
        return map_ranges

    def __avg_e_val(self, dom_B: Domain):
        """
        This function calculates the average of two E-Values (treated as psuedo P-Values) with Tippett's method.

        Parameters
        ------------
        eval_A : float
        E-value of Domain A

        eval_B : float
        E-value of Domain B

        Returns
        ------------
        avg_eval : float
        Average E-value using Fisher's method
        """

        eval_A = self.e_val

        eval_B = dom_B.e_val

        try:
            # if none of the E-values are exactly zero
            self.e_val = combine_pvalues([eval_A, eval_B], "fisher")[1]
        except:
            # else if an E-value is zero
            self.e_val = min(eval_A, eval_B)

    def __map_range_filler(self, dom_B: Domain):
        """
        This function fill gap ranges between domains which are less than `inter_gap`.

        Parmeters
        ------------
        domA : Domain
        Domain to-be filled

        domB : Domain
        Domain ranges filled up-to

        Returns
        ------------
        None
        """

        map_rng_A = self.map_range

        map_rng_B = dom_B.map_range

        tmp_save = list()
        # if the range is non-overlapping and contains a gap less than `inter_gap` fill it in
        if map_rng_A[-1] < map_rng_B[0] and (map_rng_B[0] - map_rng_A[-1]) < self.inter_gap:
            for x in range(map_rng_A[-1], map_rng_B[0]):
                if x not in map_rng_A:
                    tmp_save.append(x)
            
            self.update_map_range(tmp_save)

        # else maintain any overlaps/gaps and simply merge
        else:
            for x in dom_B.map_range:
                if x not in map_rng_A:
                    tmp_save.append(x)
            
            self.update_map_range(tmp_save)


class DomainMap(list):
    """
    Data structure with inheritance from <list> with additional functionality for handling <Domains>
    """

    def overlap_matrix(self):
        """
        Returns overlap matrix
        """

        try:
            return self._overlap_matrix
        except AttributeError:
            return self.__update_overlap_matrix()
    
    def __init_overlap_matrix(self):
        """
        Initializes new overlap matrix for DomainMap list
        """

        self._overlap_matrix = zeros(shape=(len(self), len(self)))

    def __map_range_overlapper(self):
        """
        This method fills a logical matrix of all overlapping domains in a DomainMap.
        A square matrix contains rows which are indexed to the "centeral" domain being considered,
        with columns indicating which domains overlap with it.
        For a given row `i`, r_i will always contains a logical `False` at column i, c_i;
        this is the centeral domain. All over columns in r_i could be either `True` or `False`.

        Paramters
        ------------
        self : DomainMap
        Modified list object

        Returns
        ------------
        None 
        """

        for a, dom_A in enumerate(self[:-1]):
            for b, dom_B in enumerate(self[a+1:]):
                
                # Domains that overlap more than the tolerated `overlap` could be allowed as long as the overlap is less than `overlap`
                # on both the N- and C-terminal 
                if dom_A.map_intersection(dom_B) > dom_A.overlap:
                    
                    # Check for situtations were a domain might overlap by greater than or equal to `overlap` number of residues at domain flanks
                    if dom_A.map_intersection(dom_B) <= 2*dom_A.overlap and len(dom_A.map_range) >= 2*dom_A.overlap:
                        mid_rng_idx_B = len(dom_A.map_range)//2
                        
                        if dom_A.map_intersection(dom_B, end = mid_rng_idx_B) >= dom_A.overlap or dom_A.map_intersection(dom_B, start = mid_rng_idx_B) >= dom_A.overlap:
                            self._overlap_matrix[a][b+a+1] = 1
                            self._overlap_matrix[b+a+1][a] = 1

                    # More than twice the `overlap`` and it will be marked overlapping
                    else:
                        self._overlap_matrix[a][b+a+1] = 1
                        self._overlap_matrix[b+a+1][a] = 1
                
                # Small domains (less than `overlap`) must be treated differently since their overlap could be a larger fraction of their length
                if dom_A.map_intersection(dom_B)/float(dom_A.map_len) > 0.7 or dom_A.map_intersection(dom_B)/float(dom_B.map_len) > 0.7:
                    self._overlap_matrix[a][b+a+1] = 1
                    self._overlap_matrix[b+a+1][a] = 1

    def __update_overlap_matrix(self):
        """
        Return initialized overlap matrix
        """

        self.__init_overlap_matrix()
        self.__map_range_overlapper()
        return self._overlap_matrix

    def eliminate_overlapping_domains(self):
        """
        This function organizes the overlapping domain elimination scheme
        """

        def _recursive_elimination(overlap_map, i):
            """
            This support function recursively eliminates domains that overlap

            Parameters
            ------------
            overlap_map : list
            Logical array denoting the overlapping domain of the reference domain

            i : int
            Index of the reference domain

            Returns
            ------------
            None
            """

            # If this domain has already been eliminated, skip
            if not self[i]:
                return

            # Temporarily save all index and E-value for domains that overlap
            tmp_idx = [i]
            tmp_eval = [self[i].e_val]
            for j, overlap in enumerate(overlap_map):
                if overlap and self[j] and j != i:
                    tmp_idx.append(j)
                    tmp_eval.append(self[j].e_val)

            # If the reference domain has the lowest E-value eliminate all overlapping domains 
            min_eval_idx = tmp_idx[tmp_eval.index(min(tmp_eval))]
            if  min_eval_idx == i:
                for j in tmp_idx:
                    if j != i:
                        self[j] = None
                return
            # else recursively check the overlapping domains of the lowest E-value domain
            else:
                _recursive_elimination(self._overlap_matrix[min_eval_idx], min_eval_idx)
        
        # Check all rows of the overlap matrix
        try:
            for i, overlap_row in enumerate(self._overlap_matrix):
                _recursive_elimination(overlap_row, i)
        except AttributeError:
            self.overlap_matrix()
            for i, overlap_row in enumerate(self._overlap_matrix):
                _recursive_elimination(overlap_row, i)