# Least L1-norm solution for ill-posed problem
# Author: Zhang Yunjun, 10 Jan 2019
#
# The least-norm problem is interesting only when m < n, i.e. when the equation Ax = b is underdetermined.
# It uses prior information and guess that x is more likely to be small (as measured by ||*||) than large.
# Tge least-norm problem chooses the estimate of x the one that is the smallest among all solutions that
# are consistent with measurements Ax = b.
#
# The least L1-norm problem tends to produce sparse solutions of Ax = b, often with m nonzero components.
#
# Reference:
# Boyd, S., and L. Vandenberghe (2004), Convex optimization, Cambridge university press. Chap 6.2, page 304.
# StackExchange: https://math.stackexchange.com/questions/1639716/how-can-l-1-norm-minimization-with-linear-equality-constraints-basis-pu
#


from cvxopt import glpk, matrix, sparse, spmatrix


def lst_l1(A, y, integer=False, xmax=1000):
    """

    Returns the solution of least L1-norm problem

        minimize    || x ||_1.
        subject to  A'x == y

    x can be float or integer.
    Return None if no optimal/feasible solution found.

    This problem can be converted into a linear programming problem
    by setting v = |u|, x = [u' v']',

        minimize    [0]' [u]
                    [1]  [v]

        subject to  [ I  -I]       [ 0  ]
                    [-I  -I] [u] = [ 0  ]
                    [ 0  -I] [v]   [ 0  ]
                    [ I   0]       [xmax]

                    [A]' [u] = [y]
                    [0]  [v]
    """
    m, n = A.size

    c = matrix(0.0, (2*n,1))
    c[n:] = 1.0

    # inequality constraint
    I = spmatrix(1.0, range(n), range(n))
    O = matrix(0.0, (n,n))
    G = sparse(matrix([[I, -I, O, I], [-I, -I, -I, O]]))
    h = matrix(0.0, (4*n,1))
    h[3*n:] = xmax

    # equality constraint
    Al = sparse(matrix([[A], [matrix(0.0, (m,n))]]))
    bl = y

    # solve the linear programming problem
    if integer:
        (status, x) = glpk.ilp(c, G, h, Al, bl)[0:2]
    else:
        (status, x) = glpk.ilp(c, G, h, Al, bl)[0:2]

    return x
