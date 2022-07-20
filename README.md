A minimal hitting set solver.
Given a set M and a collection of its subsets, named N,
we need to find all the minimal hitting sets.
A hitting set (hs) is a subset of M such that its
intersection with all the sets of N is non-empty.
A minimal hitting set (mhs) is a hs such that it does
not contain smaller hs

The problem in the application is formulated as follows:

A: boolean input matrix with |N| rows and |M| columns where  a_{ij} = 1  if N_i, the i_th set belonging to the collection N, has the j_th element of M (0 otherwise)
MHS: output sequence of all the MHS found given the matrix A, generated in lexicographical order.
Extra information: time required for computing the 'slow' version of the algorithm and the 'fast' one, in which some rows and cols are dropped in pre-processing
