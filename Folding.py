from __future__ import with_statement
from numpy import array, append

#####################################################
##  Base Pair Class
##
##  Basic container class for base pair info
#####################################################
class BasePair:
    def __init__(self, base, pair=0, dG=0):
	self.base = base
	self.pair = int(pair)
	self.dG = float(dG)

	if type(self.base) != str and len(self.base) != 1:
	    raise TypeError, "First argument base must be a single string character"
	else:
	    self.base = self.base.upper()

    def __str__(self):
	return self.base

    def __repr__(self):
	if self.pair == 0:
	    return self.base
	else:
	    return self.base + " => " + str(self.pair)

#####################################################
##  NAStruct Class
##
##  An array-like class for loading and manipulating
##  UNAfold data
##
#####################################################
class NAStruct:
    """
    An array-like class for loading and manipulating
    UNAfold data.

    Usage:
    
    Load a RNA structure and energy mapping
    >>> rna = NAStruct(hpin.ct, hpin.energy) # Include optional .energy file for basepairing energy

    OR

    Load theh structure, then assign bp energies
    >>> rna = NAStruct(hpin.ct)
    >>> rna.assignEnergy(hpin.energy)

    Get bp info (returns BasePair object)
    >>> rna[3]
    'C' => 'G'
    >>> rna[3].dG
    -2.7

    Get energy of basepairing as a list
    >>> rna.ddG
    [ -1.3, -2.7, -1.4, -3.1 ]

    Total hairpin energy
    >>> rna.dG
    -15.3

    """
    def __init__(self, ct_data, energy_file=None):
	self.size = None
	self.sequence = []
	self.ddG = None
	for line in ct_data:
	    data = line.split('\t')
	    if len(data) == 3:
		if self.size == None:
		    self.size = int(data[0])
		    self.dG = float(data[1].partition(' = ')[2])
		    self.name = data[2]
		else:
		    # Don't try to parse a file beyond the first structure
		    # Use parse_ct_file() for that
		    break
	    else:
		self.sequence.append( BasePair(data[1], data[4]) )

	if energy_file:
	    energy_file.seek(0)
	    self.assignEnergy(energy_file)

    def __len__(self):
	return self.size

    def assignEnergy(self, energy_file):
	self.ddG = []
    	for line in energy_file:
	    print line
    	    if not line.isspace():
    		energy = line.split('=')
    		if len(energy) == 1:
    		    (description, bases, energy) = map(str.strip, line.split(':'))

    		    if description == "Exterior":
    			continue
    		    elif description == "Stack" or description == "Hairpin":
    			bp = bases.split()[0]
    			(index, base) = bp.split('-')
    			index = int(index)

    			if str(self.__getitem__(index)) != base:
    			    raise ValueError, \
    				"Base #%d (%s) in ct file does not match base (%s) in energy file"\
    				    % (index, self.__getitem__(index), base)
    			else:
			    energy = float(energy)
    			    self.__getitem__(index).dG = energy
			    if self[index].pair > index or self[index].pair == 0:
				self.ddG += [energy]
    		else:
    		    energy = float(energy[1].strip())
    		    if energy != self.dG:
    			raise RuntimeError, "Free energy of structure in energy output file (%.2f) does not match ct file (%.2f)" % (self.dG, energy)

    		    break

    def __getitem__(self, x):
	if x<1 or x>self.size:
	    raise IndexError, "Sequence goes from 1 to " + str(self.size)

	return self.sequence[x-1]

    def details(self):
	output=''
	for (index, base) in enumerate(self.sequence):
	    output += '\t'.join( (str(index+1).zfill(3), str(base), str(base.dG)) ) + '\n'
	return output

    def __repr__(self):
	
	gap=5
	header = ''
	seq = ''
	legend = ''
	legend_chars = '*^#@$'
	stem_char = ''
	hairpin = 0
	openstem = 0
	closestem = 0

	for (index, base) in enumerate(self.sequence):
	    # Add the base to the sequence
	    seq += str(base)

	    if index % gap == 0:
		spacer = str(index+1)
		header += spacer.ljust(5)

	    if base.pair != 0:
		# when entering a new region of basepairing
		if index==0 or index > 0 and self.sequence[index-1].pair == 0:
		#    hairpin += base.pair>index+1 and 1 or -1
		    if base.pair>index+1:
			openstem += 1
			stem_char = legend_chars[openstem-1]
		    else:
			stem_char = legend_chars[openstem-closestem-1+hairpin]
			closestem += 1
			if closestem == openstem:
			    hairpin += 1
		#elif index==0:
		#    openstem += 1
		#    stem_char = legend_chars[openstem-1]

		legend += stem_char
	    else:
		legend += ' '

	return '\n'.join([header, seq, legend])

class NAList(list):
    """
    A list of NAStruct objects. Useful for loading multiple RNA structures from
    .ct and .energy files containing multiple folds.
    """
    def __init__(self, arg):
	if isinstance(arg, NAStruct):
	    arg = [ arg]
	if isinstance(arg, list):
	    super(NAList, self).__init__(arg)
	elif type(arg) == file and not arg.closed:
	    self._parse_ct_file(arg)
	else:
	    raise TypeError, "Argument to constructor must be an open file handle or a list of NAStruct %s" % type(arg)

    def _parse_ct_file(self, fh):
	file_data = []
	for line in fh:
	    data = line.split('\t')
	    if len(data) == 3 and file_data:
		self.append( NAStruct(file_data) )
		file_data = []

	    file_data += [line]

	if file_data:
	    self.append( NAStruct(file_data) )
	
    def energize(self, energy_file):
	for NA in self:
	    NA.assignEnergy(energy_file)

    def ddG(self):
	return [ hpin.ddG for hpin in self ]

    def dG(self):
	G = []
	for NA in self:
	    G += [ NA.dG ]
	return G

    def unfold(self):
	dG = array(self.dG())
	G = dG - append(dG[1:],[0])
	return G[:-1].tolist()

    def __getslice__(self,i,j):
	return NAList( super(NAList,self).__getslice__(i,j) )

    def __getitem__(self, i):
	return super(NAList,self).__getitem__(i)

    def __setitem__(self, i, y):
	if not type(y) == NAStruct:
	    raise TypeError, "Only NAStruct objects allowed in NAList"
	else:
	    return super(NAList, self).__setitem__(i,y)

    def view(self):
	output = ''
	for NA in self:
	    output += "%.2f\n" % NA.dG + NA.details() + '\n'
	return output

    def __repr__(self):
	output = ''
	for (index,item) in enumerate(self):
	    try:
		output += "\n" + str(index+1) + ": %s dG = %.1f\n" % (item.name.strip(), item.dG)
	    except AttributeError:
		output += "\n" + str(index+1)
		
	    output += '-'*len(item) + '\n'
	    output += repr(item) + '\n'
	return output

#####################################################
#####################################################
def load_rna(basename, path="./"):
    """Loads ct file and energy into an NAStructList, assuming
    the files are called basename.ct and basename.energy
    """
    if path[-1] != '/':
	path = path + '/'

    with open(path + basename+'.ct','r') as ctfile:
	rna = NAList(ctfile)

    with open(path + basename+'.energy','r') as efile:
	rna.energize(efile)

    return rna

def parse_ct_file(fh):
    """Parses ct file into list of NAStruct objects"""

    NA_list = NAList([])
    file_data = []
    for line in fh:
	data = line.split('\t')
	if len(data) == 3 and file_data:
	    NA_list.append( NAStruct(file_data) )
	    file_data = []

	file_data += [line]

    if file_data:
	NA_list.append( NAStruct(file_data) )

    return NA_list

def unfold_ddG(hpins):
    if type(hpins) != NAList:
	raise TypeError, "Argument to unfold_ddG must be an NAList"

    dG = array(hpins.dG())
    G = append(dG[1:],[0]) - dG
    return G[:-1].tolist()
