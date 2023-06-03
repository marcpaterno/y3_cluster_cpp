import numpy as np

def bias_at_nu(nu, odelta=200):
    """ Bias Tinker et a. 2010 Eqn 6
    """
    A, a, B, b, C, c = get_tinker_pars(odelta)
    bias = _bias_at_nu(nu, A, a, B, b, C, c, deltac=1.686)
    return bias

def get_tinker_pars(odelta=200):
    y = np.log10(odelta)
    # for delta=200
    tinker_best_fit = {
        'A': 1.0 + 0.24*y*np.exp(- (4/y)**4),
        'a': 0.44*y - 0.88,
        'B': 0.183,
        'b': 1.5,
        'C': 0.019+0.107*y+0.19*np.exp(-(4/y)**4),
        'c': 2.4
    }
    return [tinker_best_fit[col] for col in ['A','a','B','b','C','c']]

def _bias_at_nu(nu, A, a, B, b, C, c, deltac=1.686):
    """ Bias Tinker et a. 2010 Eqn 6
    """
    res = 1.0 - A * nu**a/ (nu**a + deltac**a)
    res+= B * nu**b
    res+= C* nu**c
    return res