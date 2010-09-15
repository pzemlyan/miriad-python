'''mirtask.uvdat - wrappers for the MIRIAD UV-data streaming API'''

# Copyright 2009, 2010 Peter Williams
#
# This file is part of miriad-python.
#
# Miriad-python is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Miriad-python is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with miriad-python.  If not, see <http://www.gnu.org/licenses/>.

import lowlevel as ll
import numpy as N
from mirtask import MiriadError, UVDataSet
from miriad import VisData

__all__ = []

# Python 2.4 compatibility: no try/finally in generators. Trying
# this raises a SyntaxError, which can't be caught within a module,
# but can be on import

try:
    from _uvdat_compat_default import _inputSets, _readFileLowlevel_gen
except SyntaxError:
    import sys
    v = sys.version_info[0] * 1000 + sys.version_info[1]
    if v >= 2005:
        # Genuine syntax error!
        raise
    del v, sys
    from _uvdat_compat_24 import _inputSets, _readFileLowlevel_gen


def init (flags, keyword='vis'):
    """Initialize standard UV data reading subsystem. If you are
writing a standalone task, you should use keys.doUvdat() rather
than this function to give the use control over whether calibration
corrections are applied or not.

Parameters:

flags - A sequence of characters giving options to the UV reading
subsystem. Possible contents are:

  r - Get a reference linetype specification via the 'ref' keyword
  s - Get Stokes parameters / polarizations via the 'stokes' keyword
  d - Perform input selection via the 'select' keyword
  l - Get a data linetype specification via the 'line' keyword
  p - Apply planet rotation and scaling
  w - Return U and V values in wavelengths
  1 - Default number of channels is 1
  x - Data must be cross-correlation data
  a - Data must be auto-correlation data
  b - Input must be a single file
  c - Apply gain/phase and delay corrections
  e - Apply polarization leakage corrections
  f - Apply bandpass corrections
  3 - Always return a 5-element 'preamble' with UVW coordinates

keyword - The keyword from which to get one or more UV dataset names.
Defaults to 'vis', the usual value.
"""
    ll.uvdatinp (keyword, flags)

class UVDatDataSet (UVDataSet):
    def __init__ (self, tno):
        self.tno = tno
        self.name = getCurrentName ()
        self.refobj = VisData (self.name)

    def _close (self):
        ll.uvdatcls ()


def inputSets ():
    # Implemented in one of _uvdat_compat_default or _uvdat_compat_24
    # thanks to compatibility issues with Python 2.4
    return _inputSets (UVDatDataSet)


def singleInputSet ():
    """Get a single DataSet object representing the visdata input set.
You should only use this function if you pass the 'b' option to
init ().
"""

    # In certain cases (e.g. an empty dataset, which yields an
    # "Invalid preamble time/baseline" MiriadError if there's no coord
    # UV variable), UVDATOPN can fail with a bug('f') call but still
    # leave the dataset opened. This means that too many successive
    # such failures will fill up MIRIAD's static buffer of UV dataset
    # information and make it impossible to open more datasets. This
    # could be solved by catching the exception and calling
    # UVDATCLS().
    #
    # On the other hand, if, say, the desired dataset doesn't exist,
    # then we get a MiriadError without having the dataset opened, and
    # if we call UVDATCLS, we segfault.
    #
    # I haven't yet figured out a way to handle both of these cases
    # cleanly. With some hacking of UVDATOPN or the uvdat module
    # in general, something could be worked out. In the meantime, we
    # prefer to avoid the segfault.
    #
    # Similar code is relevant in _uvdat_compat_*:inputSets.

    (status, tin) = ll.uvdatopn ()

    if not status:
        raise RuntimeError ('No input datasets?!')

    # Count on the user or the __del__ to close this.
    return UVDatDataSet (tin)

def readData (maxchan = 4096):
    """Generate a sequence of (preamble, data, flags) tuples representing
the visibility data in the current file. Must be called with a UVDatDataSet having
been opened, such as from calling singleInputSet() or inputSets()."""
    
    preamble = N.zeros (5, dtype=N.double)
    data = N.zeros (maxchan, dtype=N.complex64)
    flags = N.zeros (maxchan, dtype=N.int32)

    while True:
        nread = ll.uvdatrd (preamble, data, flags, maxchan)

        if nread == 0: break

        yield preamble, data[:nread], flags[:nread]

def readAll (maxchan = 4096):
    """Yield the data from all of the input datasets sequentially."""
    
    for ds in inputSets ():
        for t in readData (maxchan=maxchan):
            yield (ds, ) + t


def readFileLowlevel (fn, saveFlags, nopass=False, nocal=False, nopol=False,
                      select=None, line=None, stokes=None, ref=None,
                      maxchan=4096):
    import keys

    # Set up args

    args = ['miriad-python', 'vis=' + fn]
    flags = 'wb3'

    if select is not None:
        flags += 'd'
        args.append ('select=' + select)
    if line is not None:
        flags += 'l'
        args.append ('line=' + line)
    if stokes is not None:
        flags += 's'
        args.append ('stokes=' + stokes)
    if ref is not None:
        flags += 'r'
        args.append ('ref=' + ref)
    if not nopass: flags += 'f'
    if not nocal: flags += 'c'
    if not nopol: flags += 'e'

    # Do the actual reading -- copy what readData does for
    # greater speed.
        
    keys.init (args)
    init (flags)
    inp = singleInputSet ()

    preamble = N.zeros (5, dtype=N.double)
    data = N.zeros (maxchan, dtype=N.complex64)
    flags = N.zeros (maxchan, dtype=N.int32)

    # Python 2.4 compat ...
    return _readFileLowlevel_gen (inp, saveFlags, ll.uvdatrd, preamble,
                                  data, flags, maxchan, inp.rewriteFlags)


# Variable probes

def _getOneInt (kw):
    a = N.zeros (1, dtype=N.int32)
    ll.uvdatgti (kw, a)
    return a[0]

def getNPol ():
    """Return the number of simultaneous polarizations being returned by readData.
Zero indicates that this number could not be determined.
"""    
    return _getOneInt ('npol')

def getPols ():
    """Return the polarizations being returned by readData, an array of the size
returned by getNPol (). Zeros indicate an error. Polarization values are as in FITS
and are defined in mirtask.util.POL_??. """

    a = N.zeros (getNPol (), dtype=N.int32)
    ll.uvdatgti ('pols', a)
    return a

def getPol ():
    """Return the last Stokes parameter returned by readData. May vary from one
visibility to another."""

    return _getOneInt ('pol')

def getNChan ():
    """Return the number of channels being processed."""
    return _getOneInt ('nchan')

def getNFiles ():
    """Return the number of files being processed."""
    return _getOneInt ('nfiles')

def getVisNum ():
    """Return the current visibility number."""
    return _getOneInt ('visno')

def getVariance ():
    """Return the variance of the current visibility."""
    return ll.uvdatgtr ('variance')

def getJyPerK ():
    """Return the Jansky-per-Kelvin value of the current visibility."""
    return ll.uvdatgtr ('jyperk')

def getCurrentName ():
    """Return the name of the file currently being processed."""
    return ll.uvdatgta ('name')

def getLinetype ():
    """Return the linetype of the current visibility."""
    return ll.uvdatgta ('ltype')
