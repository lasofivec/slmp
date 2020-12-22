import slmp.derivatives as der
from random import random

def test_compute_weights2d():
    h = [0.1, 1., 2., random()]
    for th in h:
        res = comp
        assert(res[0, 0] == 0.)
        assert(res[1, 1] == - 0.0141843971631206 / th)
        assert(res[3, 2] == 0.)
        assert(res[4, 2] == 0.00797872340425532 * 0.5 / th)
    return
