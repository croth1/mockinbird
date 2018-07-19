def extract_ieclip_cy(line):
    
    chrom, pos, ref_nuc, cov, cov_str, qual_str = line.split()

    cdef int plus_pos
    cdef int minus_pos
    
    cdef int n_plus = 0
    cdef int k_plus = 0
    cdef int n_minus = 0
    cdef int k_minus = 0
    
    plus_pos = int(pos) - 1
    minus_pos = int(pos) + 1
    if cov == '0':
        return [(chrom, plus_pos, '+', n_plus, k_plus), (chrom, minus_pos, '-', n_minus, k_minus)]
    
    cdef int i = 0
    cdef Py_UNICODE s

    cdef int FWD_STRAND = 1
    cdef int REV_STRAND = 2
    cdef int END = 3
    cdef int START = 4
    cdef int state = 0
    
    cdef int cov_len = len(cov_str)

    cdef Py_UNICODE* cov_arr = <Py_UNICODE*> cov_str
    while i < cov_len:
        s = cov_arr[i]
        if s == u'^':
            i += 1
            state = START   
        elif s == u'$':
            if state == REV_STRAND:
                k_minus += 1
            state = END
        elif s in u'ACGT.':
            if state == START:
                k_plus += 1
            n_plus += 1
            state = FWD_STRAND
        elif s in u'acgt,':
            n_minus += 1
            state = REV_STRAND
        i += 1

    return [(chrom, plus_pos, '+', n_plus, k_plus), (chrom, minus_pos, '-', n_minus, k_minus)]
