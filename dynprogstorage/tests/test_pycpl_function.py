from dynprogstorage.Wrapper_dynprogstorage import Pycplfunction, Pycplfunctionvec
import math
from numpy import random


def test_pycplfunction():
    A = Pycplfunction([-1., 1.], [-10., 0.], 0.)
    A.getBreakPoints()
    A.evalf(1)
    A.evalf(0)
    A.AddSimple(1., 2., 3., 4.)
    A.getBreakPoints()

    B = Pycplfunction([-1, 1], [-math.inf, 0], 0.)
    B.getBreakPoints()


def test_pycplfunctionvec():
    C = Pycplfunctionvec([-1., 1.], [-10., 0.], [0, 0])
    C.size()
    C.Maxf_1Breaks_withO([-1., 3.], [-12., 2.], [0, 0])
    E = Pycplfunctionvec([-2., 2.], [-12., 2.], [0, 0])
    # E = E.Max_(D, C)


def test_one_breakpoint():
    n = 100
    S1 = random.uniform(low=1, high=10, size=n).tolist()
    E = Pycplfunctionvec(S1, [-math.inf] * n, [0] * n)
    PP = 1
    CC = 5 * PP
    xEtoile = E.OptimMargInt([-PP] * n, [PP] * n, [0] * n, [CC] * n)
    # print(xEtoile)


def test_two_breakpoints():
    n = 100
    S1 = random.uniform(low=1, high=10, size=n)
    S2 = S1 + 2
    E = Pycplfunctionvec(S1.tolist(), [-math.inf] * n, [0] * n,  # S1, B1, f0
                                        S2.tolist(), [0] * n)  # S2, B2
    PP = 1
    CC = 5 * PP
    xEtoile = E.OptimMargInt([-PP] * n, [PP] * n, [0] * n, [CC] * n)
