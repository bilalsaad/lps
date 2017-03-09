import numpy as np
import networkx as nx
import gmpy2 as gp

def validate(p, q):
    assert p != q, "primes must be different"
    assert p % 4 == 1, "p must be 1 mod 4"
    assert q % 4 == 1, "q must be 1 mod 4"

# Returns integer solutions <a,b,c,d> for a^2 + b^2 + c^2 + d^2 = p
# and a > 0 odd and rest are even.
def create_solutions(p):
    def viable_solution(a, b, c, d):
        return a**2 + b**2 + c**2 + d**2 == p and a > 0 and \
                a % 2 != 0 and  b % 2 == 0 and c % 2 == 0 and d % 2 == 0
    res = []
    end = p
    start = -(p - 1)
    for a in xrange(start, end):
        for b in xrange(start, end):
            for c in xrange(start, end):
                for d in xrange(start, end):
                    if viable_solution(a, b, c, d):
                        res.append((a,b,c,d))
    return res


# Given a solution created by create_solutions (a,b,c,d), i (i^2 = -1 mod q), q
# the a matrix is created as described in the article.
def create_matrix(sol, i, q):
   return np.matrix([[(sol[0] + sol[1] * i) % q, (sol[2] + i * sol[3]) % q],
       [(-1 * sol[2] + i * sol[3]) % q, (sol[0] - i * sol[1]) % q]])

def findi(q):
    for i in xrange(0, q):
        if (i**2) % q == (-1) % q:
            return i

def create_matrices(sols, q):
    i = findi(q)
    return [ create_matrix(sol, i, q) for sol in sols ]


def create_solutions_and_matrices(p, q):
    sols = create_solutions(p)
    matrices = create_matrices(sols, q)
    return sols, matrices

# returns GL(2, q). The group of 2x2 matrices with non zero determinants, over
# Fq.
def GL_2(q):
    res = []
    from itertools import product
    # Go over all possible matrices in Fq and keep those that are non singular.
    for a, b, c, d in product(xrange(0, q),
                          xrange(0, q), 
                          xrange(0, q), 
                          xrange(0, q)):
        # determinant of ( (a b) (c d) ).
        if (a*d - b*c) % q != 0:
            res.append(np.matrix([ [a, b], [c, d] ]))
    return res 

# Returns true if \exists c in Fq s.t m1 = cm2.
def PGL_equiv(m1, m2, q):
    for c in xrange(1, q):
        if ((m1 * c) % q == m2 % q).all():
            return True
    return False

def special_gl(q):
    def mat_filter(mat):
        # These are matrices : [ [1, _] [_, _] ]
        # These are matrices : [ [0, 1] [_, _] ]
        # These are matrices : [ [0, 0] [1, _] ]
        # These are matrices : [ [0, 0] [0, 1] ]
        for i in [0, 1, 2, 3]:
            if mat.item(i) > 1:
                return False
            if mat.item(i) == 1:
                return True

        return false
    # now we do some smart filtering on these.
    return filter(mat_filter, GL_2(q))

def first_non_zero_entry(m):
    for i in xrange(0, 4):
        if m.item(i) != 0:
            return m.item(i)
    assert False, "called first_non_zero_entry on zero matrix"
# This will return the LPS(2, Fq) group and an operation for it.
def PGL_2_opt(q):
    def PGL_2_op(m1, m2):
        product = np.dot(m1, m2) % q 
        normalizer = long(gp.invert(first_non_zero_entry(product), q))
        return (product * normalizer) % q
        # normalize product:
    return special_gl(q), PGL_2_op

def matrix_to_tuple(m):
    return m.item(0), m.item(1), m.item(2), m.item(3)

def create_cayley_graph(G, S, op):
    expander = nx.Graph()
    expander.add_nodes_from(matrix_to_tuple(m) for m in G)
    for mat in G:
        for s in S:
            expander.add_edge(matrix_to_tuple(mat),
                    matrix_to_tuple(op(mat, s)))
    return expander 

def create_lps(p, q):
    validate(p, q)
    solutions, matrices = create_solutions_and_matrices(p, q)
    vertices, op = PGL_2_opt(q)
    return create_cayley_graph(vertices, matrices, op)
