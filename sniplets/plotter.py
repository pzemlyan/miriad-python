#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 18:06:03 2017

@author: pete
"""

from astropy import units as u;
from astropy.io.fits.hdu.image import PrimaryHDU
from astropy.wcs import WCS
from scipy.constants import c,h,k,pi
import numpy as np
from lmfit import minimize, Parameters, Parameter, report_fit, report_errors
from math import log, atan, sin
from matplotlib.widgets import Button
from enum import Enum
from astropy.io import fits
import pickle

from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import gcf, draw_if_interactive, delaxes

from astropy import wcs
w = wcs.WCS(recs[0].header)
sp=121
subplot(sp,projection=w,slices=['x','y',0,0])
