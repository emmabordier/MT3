import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits
import os as os
from astropy import wcs
from scipy.optimize import curve_fit
from astropy.utils import NumpyRNGContext
from astropy.stats import bootstrap


def mybootstrap(Nb,f):  #Nb= nombre de bootstrap à faire
	data=fits.open('../data/Patch_KCMB_'+str(f)+'.fits')[0].data
	with NumpyRNGContext(1):
		bootresult = bootstrap(data, Nb)
	return bootresult


