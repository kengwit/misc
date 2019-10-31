class REF_FEN_SRC:
    
    def __init__(self):
        # typedef map <CONN_CODE, conn_ref_list_map_t *> conn_ref_list_map_map_t;
        # map of CONN_CODE to a MAP "conn_ref_list_map_t"
        self._crlmm = {} 

    '''
     Get a refinement node.  All connectivities defined so far
     are searched for a match.  If the match is found, the refinement
     node index indx is mapped into an index corresponding to the
     re-oriented connectivity, and the refinement node is returned;
     otherwise the connectivity
     is defined to be recognized by any subsequent searches, a new refinement
     node is generated and returned.
    '''
    def get_ref_fen (conn, indx, ref_fen_maker):
        raise NotImplementedError
         