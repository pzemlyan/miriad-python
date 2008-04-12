"""A module that provides an object that reads a UV data set's
gains table conveniently.
"""

import lowlevel as ll
import numpy as N
from lowlevel import MiriadError

__all__ = []

# Very simple wrapper classes. Shouldn't necessarily be used,
# given that there are standard APIs like uvdat*

class DataSet (object):
    """A generic Miriad data-set. Subclasses must implement a _close()
    method."""
    
    def __del__ (self):
        # tno can be None if we got an exception inside hopen,
        # or if we are deleteAll'ed

        if ll is None or not hasattr (self, 'tno'): return

        self._close ()
        
    def __repr__ (self):
        if hasattr (self, 'name'):
            return 'DataSet (%s)' % (repr (self.name))
        return 'DataSet (<unknown filename>)'

    def __str__ (self):
        if hasattr (self, 'name'):
            return '<DataSet \"%s\" handle %d>' % (self.name, self.tno)
        return '<DataSet [unknown filename] handle %d>' % (self.tno, )

    def isOpen (self):
        return hasattr (self, 'tno')
    
    def close (self):
        """Close the dataset."""

        if not hasattr (self, 'tno'): raise RuntimeError ('Trying to re-close a dataset')

        if self._histOpen: self.closeHistory ()
        
        self._close ()
        delattr (self, 'tno')
    
    def flush (self):
        """Write any changed items in the data set out to disk."""
        
        ll.hflush (self.tno)

    def deleteAll (self):
        """Completely delete this data set. After calling this function,
        this object cannot be used."""
        
        ll.hrm (self.tno)
        delattr (self, 'tno') # make any further use of this item fail
        
    def deleteItem (self, name):
        """Delete an item from this data-set."""

        ll.hdelete (self.tno, name)
    
    MODE_UNKNOWN, MODE_RD, MODE_RDWR = range (0, 3)

    def getMode (self):
        """Return the access mode of this data-set: readonly or
        read-write. See the MODE_X fields of this class for possible
        return values."""
        
        mode = ll.hmode (self.tno)

        if mode == '': return self.MODE_UNKNOWN
        elif mode == 'r': return self.MODE_RD
        elif mode == 'rw': return self.MODE_RDWR
        else: raise ValueError ('Unexpected value for "mode" argument: ' + mode)
        
        raise MiriadError ('Unknown hio mode type: ' + mode)
    
    # Data items

    def hasItem (self, name):
        """Return whether this data-set contains an item with the given name."""
        
        return ll.hexists (self.tno, name)

    def getItem (self, keyword, mode):
        """Return a DataItem object representing the desired item
        within this dataset. See the documentation of the DataItem
        constructor for the meaning of the 'keyword' and 'mode'
        parameters.
        """

        if keyword == '.': raise ValueError ("Use itemNames() instead.")
        
        return DataItem (self, keyword, mode)

    def itemNames (self):
        """Generate a list of the names of the data items contained in
        this data set."""
        
        ilist = DataItem (self, '.', 'r')
        s = ilist.getSize ()
        
        while ilist.getPosition () < s:
            yield ilist.seqReadString ()

        del ilist

    # History

    _histOpen = False
    
    def openHistory (self, mode='a'):
        """Open the history item of this data set. 'mode' may be 
        'r' if the history is being read, 'w' for truncation and writing,
        and 'a' for appending. The default is 'a'.
        """
        
        if mode == 'r': modestr = 'read'
        elif mode == 'w': modestr = 'write'
        elif mode == 'a': modestr = 'append'
        else: raise ValueError ('Unexpected value for "mode" argument: ' + mode)

        ll.hisopen (self.tno, modestr)
        self._histOpen = True

    def writeHistory (self, text):
        """Write text into this data set's history file."""
        
        ll.hiswrite (self.tno, text)

    def logInvocation (self, taskname, args=None):
        """Write text into this data set's history file logging the invocation
        of this task: when it was run and what parameters it was given. Can
        optionally be given an argument list if that contained in sys.argv
        does not represent this task."""

        ll.hisinput (self.tno, taskname, args)
    
    def closeHistory (self):
        """Close this data set's history item."""

        ll.hisclose (self.tno)
        self._histOpen = False

    # Header variables

    def getHeaderFloat (self, keyword, default):
        """Retrieve the value of a float-valued header variable."""

        return ll.rdhdr (self.tno, keyword, float (default))

    def getHeaderInt (self, keyword, default):
        """Retrieve the value of an int-valued header variable."""

        return ll.rdhdi (self.tno, keyword, int (default))

    def getHeaderBool (self, keyword, default):
        """Retrieve the value of a bool-valued header variable."""

        return bool (ll.rdhdl (self.tno, keyword, int (default)))

    def getHeaderDouble (self, keyword, default):
        """Retrieve the value of a double-valued header variable."""

        return ll.rdhdd (self.tno, keyword, float (default))

    def getHeaderComplex (self, keyword, default):
        """Retrieve the value of a complex-valued header variable."""

        dc = complex (default)
        out = ll.rdhdc (self.tno, keyword, (dc.real, dc.imag))
        return complex (out[0], out[1])
    
    def getHeaderString (self, keyword, default):
        """Retrieve the value of a string-valued header variable.
        Maximum value length is 512."""

        return ll.rdhda (self.tno, keyword, str (default))

    def copyHeader (self, dest, keyword):
        """Copy a header variable from this data-set to another."""

        ll.hdcopy (self.tno, dest.tno, keyword)

    # skip hdprsnt: same thing as hexists
    
    def getHeaderInfo (self, keyword):
        """Return the characteristics of the header variable. Returns:
        (desc, type, n), where 'desc' describes the item or gives its value
        if it can be expressed compactly; 'type' is one of 'nonexistant',
        'integer*2', 'integer*8', 'integer', 'real', 'double', 'complex',
        'character', 'text', or 'binary'; and 'n' is the number of elements
        in the item. If 'n' is 1, then 'desc' encodes the item's value.
        """

        (desc, type, n) = ll.hdprobe (self.tno, keyword)

        if n == 0: raise MiriadError ('Error probing header ' + keyword)

        return (desc, type, n)

class DataItem (object):
    """An item contained within a Miriad dataset."""
    
    def __init__ (self, dataset, keyword, mode):
        self.dataset = dataset
        self.name = keyword

        if mode == 'r': modestr = 'read'
        elif mode == 'w': modestr = 'write'
        elif mode == 'a': modestr = 'append'
        elif mode == 's': modestr = 'scratch'
        else: raise ValueError ('Unexpected value for "mode" argument: ' + mode)

        self.itno = ll.haccess (dataset.tno, keyword, modestr)

    def __del__ (self):
        # itno can be None if we got an exception inside hopen.

        if ll is None or not hasattr (self, 'itno'): return

        ll.hdaccess (self.itno)

    def getSize (self):
        """Return the size of this data item."""

        return ll.hsize (self.itno)

    def seek (self, offset):
        """Seek to the specified position within this data item."""

        ll.hseek (self.itno, int (offset))

    def getPosition (self):
        """Retrieve the current position within this data item."""

        return ll.htell (self.itno)

    def seqReadString (self):
        """Read until newline from the current position within this
        data item. Maximum string length of 512."""

        return ll.hreada (self.itno)

    def seqWriteString (self, line, length=None):
        """Write a textual string into the data item, terminating
        the string with a newline. If desired, only a subset of the
        string can be written out; the default is to write the
        entire string."""

        if length is None: length = len (line)
        ll.hwritea (self.itno, str (line), length)

    # Reading buffers
    
    def readBytes (self, buf, offset, length=None):
        """Read an array of bytes from the given location in the data
        item. The default read length is the size of the array."""

        if length is None: length = darray.size
        ll.hreadb (self.itno, buf, offset, length)

    def readInts (self, buf, offset, length=None):
        """Read an array of integers from the given location in the data
        item. The default read length is the size of the array."""

        if length is None: length = darray.size
        ll.hreadi (self.itno, buf, offset, length)

    def readShorts (self, buf, offset, length=None):
        """Read an array of 16-bit integers from the given location in the data
        item. The default read length is the size of the array."""

        if length is None: length = darray.size
        ll.hreadj (self.itno, buf, offset, length)

    def readLongs (self, buf, offset, length=None):
        """Read an array of 64-bit integers from the given location in the data
        item. The default read length is the size of the array."""

        if length is None: length = darray.size
        ll.hreadl (self.itno, buf, offset, length)

    def readFloats (self, buf, offset, length=None):
        """Read an array of floats from the given location in the data
        item. The default read length is the size of the array."""

        if length is None: length = darray.size
        ll.hreadr (self.itno, buf, offset, length)

    def readDoubles (self, buf, offset, length=None):
        """Read an array of doubles from the given location in the data
        item. The default read length is the size of the array."""

        if length is None: length = darray.size
        ll.hreadd (self.itno, buf, offset, length)

    def readComplex (self, buf, offset, length=None):
        """Read an array of complexes from the given location in the data
        item. The default read length is the size of the array."""

        if length is None: length = darray.size
        ll.hreadc (self.itno, buf, offset, length)

    # Writing
    
    def writeBytes (self, buf, offset, length=None):
        """Write an array of bytes to the given location in the data
        item. The default write length is the size of the array."""

        if length is None: length = darray.size
        ll.hwriteb (self.itno, buf, offset, length)

    def writeInts (self, buf, offset, length=None):
        """Write an array of integers to the given location in the data
        item. The default write length is the size of the array."""

        if length is None: length = darray.size
        ll.hwritei (self.itno, buf, offset, length)

    def writeShorts (self, buf, offset, length=None):
        """Write an array of 16-bit integers to the given location in the data
        item. The default write length is the size of the array."""

        if length is None: length = darray.size
        ll.hwritej (self.itno, buf, offset, length)

    def writeLongs (self, buf, offset, length=None):
        """Write an array of 64-bit integers to the given location in the data
        item. The default write length is the size of the array."""

        if length is None: length = darray.size
        ll.hwritel (self.itno, buf, offset, length)

    def writeFloats (self, buf, offset, length=None):
        """Write an array of floats to the given location in the data
        item. The default write length is the size of the array."""

        if length is None: length = darray.size
        ll.hwriter (self.itno, buf, offset, length)

    def writeDoubles (self, buf, offset, length=None):
        """Write an array of doubles to the given location in the data
        item. The default write length is the size of the array."""

        if length is None: length = darray.size
        ll.hwrited (self.itno, buf, offset, length)

    def writeComplex (self, buf, offset, length=None):
        """Write an array of complexes to the given location in the data
        item. The default write length is the size of the array."""

        if length is None: length = darray.size
        ll.hwritec (self.itno, buf, offset, length)

__all__ += ['DataSet', 'DataItem']

class UserDataSet (DataSet):
    def __init__ (self, fname, create=False):
        if create: mode = 'new'
        else: mode = 'old'

        self.tno = ll.hopen (fname, mode)
        self.name = fname
        
    def _close (self):
        ll.hclose (self.tno)

__all__ += ['UserDataSet']

class UVDataSet (DataSet):
    def __init__ (self, fname, mode):
        if mode == 'r': modestr = 'old'
        elif mode == 'w': modestr = 'new'
        elif mode == 'a': modestr = 'append'
        else: raise ValueError ('Unexpected mode string ' + mode)

        self.tno = ll.uvopen (fname, modestr)
        self.name = fname

    def _close (self):
        ll.uvclose (self.tno)

    # These override the basic DataSet operations

    def flush (self):
        """Write out any unbuffered changes to the UV data set."""
        
        ll.uvflush (self.tno)

    # UV-specific operations

    def next (self):
        """Skip to the next UV data record. On write, this causes an
        end-of-record mark to be written."""

        ll.uvnext (self.tno)

    def rewind (self):
        """Rewind to the beginning of the file, allowing the UV data to
        be reread from the start."""

        ll.uvrewind (self.tno)

    def write (self, preamble, data, flags, length=None):
        """Write a visibility record consisting of the given preamble,
        data, flags, and length. Length defaults to the length of the
        flags array."""

        if length is None: length = flags.size

        ll.uvwrite (self.tno, preamble, data, flags, length)

    def rewriteFlags (self, flags):
        """Rewrite the channel flagging data for the current
        visibility record. 'flags' should be a 1D integer ndarray of the
        same length and dtype returned by a uvread call."""

        ll.uvflgwr (self.tno, flags)
        
    # uvset exploders

    def _uvset (self, object, type, n, p1, p2, p3):
        ll.uvset (self.tno, object, type, n, p1, p2, p3)

    def setPreambleType (self, *vars):
        """Specify up to five variables to put in the preamble block.
        Should be given a list of variable names; 'uv' and 'uvw' are
        a special expansion of 'coord' that expand out to their
        respective UV coordinates. Default list is 'uvw', 'time',
        'baseline'."""
        
        self._uvset ('preamble', '/'.join (vars), 0, 0., 0., 0.)

    def setSelectAmplitude (self, selamp):
        """Specify whether selection based on amplitude should be
        performed."""

        if selamp: val = 1
        else: val = 0
        
        self._uvset ("selection", "amplitude", val, 0., 0., 0.,)
        
    def setSelectWindow (self, selwin):
        """Specify whether selection based on window should be
        performed."""

        if selwin: val = 1
        else: val = 0
        
        self._uvset ("selection", "window", val, 0., 0., 0.,)

    def setPlanetParams (self, major, minor, angle):
        """Set the reference parameters for planet scaling and
        rotation."""

        self._uvset ("planet", "", 0, major, minor, angle)
    
    def setWavelengthMode (self, wlmode):
        """Specify that UV coordinates should be returned in units
        of wavelength. Otherwise, they are returned in nanoseconds."""

        if wlmode:
            self._uvset ("coord", "wavelength", 0, 0., 0., 0.)
        else:
            self._uvset ("coord", "nanosec", 0, 0., 0., 0.)

    def setCorrelationType (self, type):
        """Set the correlation type that will be used in this
        vis file."""

        self._uvset ("corr", type, 0, 0., 0., 0.)
    
    # oh god there are a bunch more of these: data linetype, refernce
    # linetype, gflag, flags, corr
    
    # Variable handling

    def copyMarkedVars (self, output):
        """Copy variables in this data set to the output data set. Only
        copies those variables which have changed and are marked as
        'copy'."""

        ll.uvcopyvr (self.tno, output.tno)

    def updated (self):
        """Return true if any user-specified 'important variables' have
        been updated in the last chunk of data read."""

        return bool (ll.uvupdate (self.tno))

    def initVarsAsInput (self, linetype):
        """Initialize the UV reading functions to copy variables from
        this file as an input file. Linetype should be one of 'channel',
        'wide', or 'velocity'. Maps to Miriad's varinit() call."""

        ll.varinit (self.tno, linetype)

    def initVarsAsOutput (self, input, linetype):
        """Initialize this dataset as the output file for the UV
        reading functions. Linetype should be one of 'channel', 'wide',
        or 'velocity'. Maps to Miriad's varonit() call."""

        ll.varonit (input.tno, self.tno, linetype)

    def copyLineVars (self, output):
        """Copy UV variables to the output dataset that describe the
        current line in the input set."""

        ll.varcopy (self.tno, output.tno)

    def makeVarTracker (self):
        """Create a UVVarTracker object, which can be used to track
        the values of UV variables and when they change."""
        
        return UVVarTracker (self)

    def probeVar (self, varname):
        """Get information about a given variable. Returns (type, length,
        updated) or None if the variable is undefined.

        type - The variable type character: a (text), r ("real"/float),
        i (int), d (double), c (complex)

        length - The number of elements in this variable; zero if unknown.

        updated - True if the variable was updated on the last UV data read.
        """

        (type, length, updated) = ll.uvprobvr (self.tno, varname)

        if type == '' or type == ' ': return None
        return (type, length, updated)

    def getVarString (self, varname):
        """Retrieve the current value of a string-valued UV
        variable. Maximum length of 512 characters."""

        return ll.uvgetvra (self.tno, varname)
    
    def getVarInt (self, varname, n=1):
        """Retrieve the current value or values of an int-valued UV
        variable."""

        ret = ll.uvgetvri (self.tno, varname, n)

        if n == 1: return ret[0]
        return N.asarray (ret, dtype=N.int)
        
    def getVarFloat (self, varname, n=1):
        """Retrieve the current value or values of a float-valued UV
        variable."""

        ret = ll.uvgetvrr (self.tno, varname, n)

        if n == 1: return ret[0]
        return N.asarray (ret, dtype=N.float32)

    def getVarDouble (self, varname, n=1):
        """Retrieve the current value or values of a double-valued UV
        variable."""

        ret = ll.uvgetvrd (self.tno, varname, n)
    
        if n == 1: return ret[0]
        return N.asarray (ret, dtype=N.double64)

    def getVarComplex (self, varname, n=1):
        """Retrieve the current value or values of a complex-valued UV
        variable."""

        ret = ll.uvgetvrc (self.tno, varname, n)
    
        if n == 1: return ret[0]
        return N.asarray (ret, dtype=N.complex64)

    def getVarFirstString (self, varname, dflt):
        """Retrieve the first value of a string-valued UV
        variable with a default if the variable is not present.
        Maximum length of 512 characters."""

        return ll.uvrdvra (self.tno, varname, dflt)
    
    def getVarFirstInt (self, varname, dflt):
        """Retrieve the first value of an int-valued UV
        variable with a default if the variable is not present."""

        return ll.uvrdvri (self.tno, varname, dflt)
    
    def getVarFirstFloat (self, varname, dflt):
        """Retrieve the first value of a float-valued UV
        variable with a default if the variable is not present."""

        return ll.uvrdvrr (self.tno, varname, dflt)
    
    def getVarFirstDouble (self, varname, dflt):
        """Retrieve the first value of a double-valued UV
        variable with a default if the variable is not present."""

        return ll.uvrdvrd (self.tno, varname, dflt)
    
    def getVarFirstComplex (self, varname, dflt):
        """Retrieve the first value of a complex-valued UV
        variable with a default if the variable is not present."""

        dflt = complex (dflt)
        retval = ll.uvrdvrd (self.tno, varname, (dflt.real, dflt.imag))
        return complex (retval[0], retval[1])
    
    def trackVar (self, varname, watch, copy):
        """Set how the given variable is tracked. If 'watch' is true, updated()
        will return true when this variable changes after a chunk of UV data
        is read. If 'copy' is true, this variable will be copied when
        copyMarkedVars() is called.
        """

        switches = ''
        if watch: switches += 'u'
        if copy: switches += 'c'

        ll.uvtrack (self.tno, varname, switches)

    def scanUntilChange (self, varname):
        """Scan through the UV data until the given variable changes. Reads
        to the end of the record in which the variable changes. Returns False
        if end-of-file was reached, True otherwise."""

        return ll.uvscan (self.tno, varname) == 0

    def writeVarInt (self, name, val):
        """Write an integer UV variable. val can either be a single value or
        an ndarray for array variables."""
        
        if not isinstance (val, N.ndarray):
            v2 = N.ndarray (1, dtype=N.int32)
            v2[0] = int (val)
            val = v2

        ll.uvputvri (self.tno, name, val)
    
    def writeVarFloat (self, name, val):
        """Write an float UV variable. val can either be a single value or
        an ndarray for array variables."""
        
        if not isinstance (val, N.ndarray):
            v2 = N.ndarray (1, dtype=N.float32)
            v2[0] = float (val)
            val = v2

        ll.uvputvrr (self.tno, name, val)
    
    def writeVarDouble (self, name, val):
        """Write a double UV variable. val can either be a single value or
        an ndarray for array variables."""
        
        if not isinstance (val, N.ndarray):
            v2 = N.ndarray (1, dtype=N.float64)
            v2[0] = float (val)
            val = v2

        ll.uvputvrd (self.tno, name, val)
    
class UVVarTracker (object):
    def __init__ (self, owner):
        self.dataset = owner
        self.vhnd = ll.uvvarini (owner.tno)

    def track (self, *vars):
        """Indicate that the specifieds variable should be tracked by this
        tracker."""

        for var in vars:
            ll.uvvarset (self.vhnd, var)

    def copyTo (self, output):
        """Copy the variables tracked by this tracker into the output
        data set."""

        ll.uvvarcpy (self.vhnd, output.tno)

    def updated (self):
        """Return true if one of the variables tracked by this tracker
        was updated in the last UV data read."""

        return bool (ll.uvvarupd (self.vhnd))

__all__ += ['UVDataSet', 'UVVarTracker']
