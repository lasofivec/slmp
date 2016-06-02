import numpy as np
from numpy import pi
from globals_variables import center_adv_x
from globals_variables import center_adv_y

def initialize_distribution(name):
    """ Initializes distribution function, with type defined in name.
    Args:
      name (str) : type of the distribution function."""
    
    if (name == "CONSTANT_DISTRIBUTION") :
        func_init = lambda x,y : constant_distribution(x,y)
    elif (name == "GAUSSIAN_DISTRIBUTION") :
        func_init = lambda x,y : gaussian_pulse(x, y, \
                                                xc=center_adv_x, \
                                                yc=center_adv_y, \
                                                sig = 0.04, \
                                                max = 1.0)
    elif (name == "SINUSOIDAL_WAVES") :
        func_init = lambda x,y : sinusoidal_waves(x, y, maxval=1., along="x")
    elif (name == "DIOCOTRON_DISTRIBUTION") :
        func_init = lambda x,y : sinusoidal_waves(x, y, eta=0.1, l=7., \
                                                  rmin=0.7, rmax=0.8,  \
                                                  maxval=1.1)
    elif (name == "COS_SIN_DISTRIBUTION") :
        func_init = lambda x, y : cossin_distribution(x,y)
    else:
        raise SystemExit("ERROR initialize_distribution():" + \
                         " Undefined type of distribution function")

    return func_init


def constant_distribution(x,y,val=0.5):
    """ Constant distribution function.
    Args:
      x   (real) : x-coordinate of point where function will be evaluated.
      y   (real) : y-coordinate of point where function will be evaluated.
      val (real) : value of distribution. Default =0.5. """

    res = np.zeros_like(x)
    if (np.size(res) > 1):
        res[:] = val
    else:
        res = val
    return res


def gaussian_pulse(x, y, xc, yc, sig, max):
    """ Distribution function as a gaussian pulse.
    Args:
      x   (real) : x-coordinate of point where function will be evaluated.
      y   (real) : y-coordinate of point where function will be evaluated.
      xc  (real) : x-coordinate of center of gaussian.
      yc  (real) : y-coordinate of center of gaussian.
      sig (real) : width of gaussian.
      max (real) : amplitude of gaussian. """

    res = max * np.exp(-0.5 * ((np.mod(x+1.0, 2.)-1.-xc)**2/sig**2 + (y-yc)**2/sig**2))
    return res


def sinusoidal_waves(x, y, maxval=1., along="x"):
    """ Distribution function in form of sinusoidal waves either along x or y.
    Args:
      x   (real) : x-coordinate of point where function will be evaluated.
      y   (real) : y-coordinate of point where function will be evaluated.
      maxval (real) : maxvalue of distribution. Default =1.0
      along (string) : either "x" or "y" direction. Default "x". """
    if (along == "x"):
        res = maxval/2.*np.sin(2.*pi*x) + maxval/2.
    elif (along == "y"):
        res = maxval/2.*np.sin(2.*pi*y) + maxval/2.
    elif (along == "xy"):
        res = maxval/2.*np.sin(2.*pi*y)*np.cos(2.*pi*x) + maxval/2.
    else:
        raise SystemExit("ERROR sinusoidal_waves(): undefined parameter along")
    return res

def diocotron_distribution(x, y, eta, l, rmin, rmax, maxval):
    r = lambda x,y : np.sqrt(x**2+y**2)
    th = lambda x,y : np.arctan2(y,x) % (2.*pi)
    if (np.size(x) > 1) :
        res = np.zeros_like(x)
        res = 1.+eta*np.cos(l*th(x,y))
        res[np.where((r(x,y) > rmax) | (r(x,y) < rmin)) ] = 0.
    else :
        res = 1.+eta*np.cos(l*th(x,y))
        if (r(x,y) > rmax) or (r(x,y) < rmin) :
            res = 0.
    return res

def cossin_distribution(x,y):
    res = np.cos(2.0*pi*x)*np.sin(2.0*pi*y)
    return res