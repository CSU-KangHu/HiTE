# Supplemental Data. Wenke et al. (2011) Plant Cell 10.1105/tpc.111.088682.

# Supplemental Dataset 1. Python script of the SINE-Finder.

# sine_finder

# meta

__program_name__ = "sine_finder.py"
__version_no__ = "1.0.1"
__version_date__ = "2010-09-28"


#################
# CONFIGURATION #
#################

# Default values for commandline options

options = {
	'verbose': 0,
	'RUNTYPE': 'seqwise',                   # possible: 'chunkwise'
	'TSD_MISMATCH_TOLERANCE': 2,
	'TSD_MIN_WORDSIZE': 5,
	'TSD_MISMATCH_PENALTY': 1,
	'TSD_SCORE_CUTOFF': 10,
	'TSD_ORIENTATION': 'F',                 # possible: 'FR':
	                                        # 'R' equal to 'F' but as 
	                                        # a result of reverse search  
	'CHUNKSIZE': 100000,
	'OVERLAP': 8000,
	'OUTTYPE': 'fasta',
	}
EXITCHAR = "qQ"

# File extensions shown in the "file browser"

extensions = ('fas','fasta','mfa',)


###########
# IMPORTS #
###########

import string, re, os, sys, time
from functools import reduce


###############################
# CONSTANTS, GLOBAL VARIABLES #
###############################

# Detailed description of command line arguments

__description__ = """
  sine_finder [options] <fastafile_name>

  use -h for a detailed description of command line arguments

  DESCRIPTION:

      sine_finder.py is a script for the analysis of sequence 
      files for putative members of the SINE repeat family.

      The program is written in Python (http://www.python.org/).
      Besides an installation of Python (at least version 2.5)
      there are no other requirements. The script contains a
      basic MotifFinder class, which performs a pattern search
      by regular expressions. A derived SINEFinder class
      extends the functionality to examine TSDs (target site
      duplications). Search patterns can easily be modified.

      The sine_finder may be used in 3 ways:
         - command line mode by calling with arguments
         - interactive mode by calling without arguments
         - as a module by including it in other scripts

      Sequences to be analysed should be in FASTA format: 
         - a record starts with an ID line having '>' as 
           first character,
         - descriptive data may follow the ID delimited 
           by spaces,
         - sequence might be wrapped and displayed over 
           several lines,
         - multiple sequence entries are allowed,
         - empty lines between records are allowed.
      For datasets with very long sequences the option 'chunk-
      wise' may be used. Then at once sequences are only 
      read and processed in segments of a specified length 
      (chunks). To avoid loss of matches at split sites chunks
      can overlap (size of overlap is definable, duplicate
      matches are removed).

  STRATEGY

      The script performs a search in two steps: 

         1. A pattern search. The pattern is defined as a 
            regular expression, which is a textual description
            for a set of strings. It consists of groups or
            subpatterns so that the result contains the
            sequence matches for these as well. 

            As example: the a-box pattern '[GA][CGA]TGG'
            may result in ACTGG or GATGG whereas a spacer
            pattern like '.{25,50}' matches any sequence as
            long as it has a length between 25 and 50 bp.

         2. A TSD cannot be found with a regular expression, 
            so in the pattern search regions have to be defined
            where the TSDs are expected. In a second step these 
            matches were examined for TSDs.

            Each result of the pattern search produces two
            TSD expectation regions which are compared by a fuzzy
            matching. Doing this a window of the size of MIN_WORDSIZE 
            shifts over sequence 1 and a search for exact 
            counterparts of these 'words' in sequence 2 is
            performed. Matches and their coordinates are 
            returned. These serve as seed for the following
            procedure where the given matches are extended to the 
            right and to the left until MISMATCH_TOLERANCE is reached.
            Terminal mismatches are clipped.

      (For additional information on regular expressions 
      and their syntax: http://docs.python.org/library/re.html)
"""

__usage__ = """
  sine_finder [options] <fastafile_name>

  OPTIONS:

     -h     (help:) display this message.
     -d     (description:) display a short description of the 
            program.
     -v     (version:) print version number.

     EVALUATION PARAMETERS:

     -T (seqwise|chunkwise)
            (run type:) the way how infiles are processed (default:
            %(RUNTYPE)s).

            'seqwise': each sequence is loaded and searched for patterns.
               For large sequences 'chunkwise' is the better choice.
            'chunkwise': only fragments (chunks) of the sequence are
               loaded and processed. The size is defined by option
               -C. To ensure that matches do not get lost by sequence 
               splitting an overlap should be specified (option -O).

     -t <integer>
            TSD mismatch tolerance (default:%(TSD_MISMATCH_TOLERANCE)s).

     -w <integer>
            word size TSD seed starts search with (default:%(TSD_MIN_WORDSIZE)s).

     -p <integer>
            penalty for a nucleotide mismatch in TSD search (default:
            %(TSD_MISMATCH_PENALTY)s).

     -s <integer>
            TSD score cutoff (default:%(TSD_SCORE_CUTOFF)s).

     -o (F|R|FR) 
            direction of TSD search, allowed orientation (default:%(TSD_ORIENTATION)s).

     (only chunkwise processing:)

     -C <integer>
            (chunksize:) size of each fragment loaded and processed
            individually. Use only when -T is set to chunkwise (default: 
            %(CHUNKSIZE)s).

     -O <integer>
            (overlap:) overlap of fragments treated by chunkwise processing.
            Use only when -T is set to chunkwise (default:%(OVERLAP)s).

     OTHER OPTIONS:

     -f (fasta|csv|both)
           (file type:) file type result is written to (default:%(OUTTYPE)s).

     -V    (verbose:) display program call.

"""

# String formats for version and error messages

__version__ = "\nVERSION: %s (%s)\n\n"
__error_msg__ = "\nERROR %02d: %s\n"

# Messages for defined errors

__errorcodes__ = {
	0: "unrecognized error",
	1: "'%s' is not a file",
	2: "wrong or missing option for '%s'",
	3: "'%s' seems to be no fasta file",
	12: "during blast",
	13: "two filenames required",
	51: "argument for '%s' must be %s",
	53: "overlap should be smaller than chunk size",
	100: "not yet supported"
	}


#############
# FUNCTIONS #
#############

def exit_on_error(number, *args):
	"""Output on error."""

	sys.stderr.write("sine_finder_version=%s\n" % __version_no__)
	if args:
		error_txt = __errorcodes__[number] % tuple(args)
	else:
		error_txt = __errorcodes__[number]
	sys.stderr.write(__error_msg__ % (number, error_txt))

	sys.exit()
	
def prompter(oname, **kwargs):
	"""
	For printing of a prompt under interactive mode.
	
	The function supports:
	   - check for required arguments
	   - validation of input
	   - selection of a set of options
	   - supplementation of options
	"""

	# Prompting for a valid input
	
	while 1:

		# Construction of a note for special information
		
		l = []
		if kwargs.setdefault('input_required', False) is True:
			l.append("input_required")
		if kwargs.setdefault('note', ''):
			l.append(kwargs['note'])
		if l:
			app = " (%s)" % ', '.join(l)
		else:
			app = ""

		# Prompt with options: display the set of options
		
		if kwargs.setdefault('options', []):
			print()
			for i in range(len(kwargs['options'])):
				if not kwargs['options'][i]:
					print()
				print("\t(% 3d) %s" % (i, kwargs['options'][i][1]))

		# The data input.
		# For to leave the program you should interrupt by CTRL+c!
		print()
		try:
			choice = input("%s%s: " % (oname, app))
			if choice and choice in EXITCHAR:
				#print "'%s'" % choice, "'%s'" % EXITCHAR 
				sys.exit("Bye!")
		except KeyboardInterrupt:
			choice = input("\nquit (q) OR continue (any input)? ")
			if choice == 'q':
				sys.exit("Bye!")
			continue
		if not choice and 'default' in kwargs:
			choice = kwargs['default']

		# Evaluation of chosen option
		
		if kwargs.setdefault('options', []):
			if kwargs.setdefault('allow_raw', False) is True and choice.startswith('='):
				choice = choice[1:]
			else:
				try:
					choice = kwargs['options'][int(choice)][0]
				except ValueError:
					continue
				except IndexError:
					continue

		if not kwargs.setdefault('validate', lambda v: True)(choice):
			continue

		# The while-loop will be left when the choice is accepted or no input is required.

		if choice or kwargs['input_required'] is False:
			return choice

# Lambda-functions

is_integer = lambda v: type(v) == int or re.match("^\d+$", v)
is_integer_argv = lambda o: re.match("^\d+$", sys.argv[sys.argv.index(o)+1])
is_float = lambda v: type(v) in (int, float) or re.match("^\d+\.?\d*$", v)
is_float_argv = lambda o: re.match("^\d+\.?\d*$", sys.argv[sys.argv.index(o)+1])
get_extension = lambda f: f.split('.')[-1]
filter_extension = lambda f: get_extension(f) in extensions
seems_to_be_fasta = lambda f: open(f).read(10).strip().startswith('>')

###########
# CLASSES #
###########

class FastaIterator:
	
	config = {
		'RUNTYPE': 'seqwise',
		'CHUNKSIZE': 100000,
		'OVERLAP': 8000,
		}

	def __init__(self, filename, **kwargs):
		
		self.filename = filename
		self.fi = open(self.filename)
		for k in self.__class__.config:
			if k in kwargs:
				setattr(self, k, kwargs[k])
			else:
				setattr(self, k, self.__class__.config[k])
		self.id = ''
		self.seq = ''
		self.cache = ''
		self.startpos = 0
	
	def _remaining_data(self):
		
		remain_len = self.CHUNKSIZE - len(self.seq)
		if remain_len > 0:
			self.seq += re.sub("\s", "", self.fi.readline(remain_len))
			if self.CHUNKSIZE - len(self.seq) > 0:
				return True
		return False

	def _read_chunk(self):
		
		# id
		# - get first character of next line
		if self.cache:
			fchar = self.cache
			self.cache = ''
		else:
			while 1:
				fchar = self.fi.readline(1)
				if not fchar:
					return False # EOF, no new sequence
				if not fchar.strip():
					continue
				break
		# - check if a new entry starts
		if fchar == '>':
			# new entry -> empty vars
			id_line = self.fi.readline()
			self.id = id_line.strip().split()[0]
			self.seq = ''
			self.startpos = 0
		# sequence
		else:
			# continue with old entry ->
			# - use self.id
			# - put overlap in seq
			# - calculate pos with overlap
			self.startpos = self.startpos + self.CHUNKSIZE - self.OVERLAP
			self.seq = self.seq[-self.OVERLAP:] + fchar
			# add rest of the line
			if not self._remaining_data():
				# chunksize reached (overlap + 1 = chunksize!)
				return True
		# now read next lines until line starts with ">" or chunksize reached
		while 1:
			fchar = self.fi.readline(1)
			if not fchar:
				return True # EOF, evaluate sequence
			if not fchar.strip():
				continue
			if fchar == '>':
				# new entry -> put to cache
				# sequence entry can be evaluated
				self.cache = fchar
				return True
			else:
				self.seq += fchar
				if not self._remaining_data():
					# chunksize reached
					return True
	
	def _read_sequence(self):
		
		# id
		while 1:
			if self.cache:
				line = self.cache
				self.cache = ''
			else:
				line = self.fi.readline()
			if not line:
				return False # EOF, no new sequence
			elif not line.strip():
				continue
			m = re.match("^>(\S+)", line)
			if m:
				self.id = m.group(1)
				break
			else:
				raise Exception("file not in valid format")
		self.seq = ''
		# sequence
		while 1:
			line = self.fi.readline()
			if not line:
				return True # EOF, evaluate sequence
			if not line.strip():
				continue
			if line.startswith('>'):
				# new entry -> put to cache
				# sequence entry can be evaluated
				self.cache = line
				return True
			else:
				self.seq += re.sub("\s", "", line)
	
	def __next__(self):
		
		if self.RUNTYPE == 'seqwise':
			if self._read_sequence():
				return self.id, self.seq, self.startpos
		elif self.RUNTYPE == 'chunkwise':
			if self._read_chunk():
				return self.id, self.seq, self.startpos

	def close(self):
		
		self.fi.close()


class MotifFinder:
	"""
	Base class for motif finders.
	
	The sought-after motif can be divided into several 
	submotifs. The patterns (submotifs) can be defined
	as regular expressions and must be named each to 
	address them as group. Every segment between the
	submotifs must be defined as spacer: '.{minlen,maxlen}'.
	
	Patterns should be passed like this:
	   p = ((name1, pattern1), (name2, pattern2), ...)
	If you pass only one pattern a comma is required:
	   p = ((name1, pattern1),)
	
	Useage:
	    mf = MotifFinder()
	    mf.set_pattern(p)
	    mf.set_sequence(seq)
	    mf.search('FR')
	    => results are stored in mf.mm
	"""

	config = {
		'verbose': 0,
		}

	pattern = (("complete sequence", "^.*$"),)

	outformat = {
		'fasta': ">%(name)s %(direct)s %(start)s-%(end)s\n%(seq)s\n",
		'csv': "%(name)s\t%(direct)s\t%(start)s\t%(end)s\t%(seq)s\n",
		}
	prepend = {
		'fasta': "",
		'csv': "name\tdirect\tstart\tend\tsequence\n",
		}

	# Constructor of the class
	
	def __init__(self, **kwargs):

		# Configuration

		self.cfg = self.__class__.config.copy()
		for k in kwargs:
			if k in self.cfg:
				self.cfg[k] = kwargs[k]

		# Pattern

		self.set_pattern(self.__class__.pattern)

		# Variables

		self.seq = ''
		self.rseq = ''
		self.seqlen = 0
		self.mm = []
		self.status = 0
		self.header_written = False
	
	# PRIVATE METHODS

	def _construct_pattern(self, p):
		"""Concatenate the given re-patterns to one."""

		return reduce(lambda a, b: a + "(?P<%s>%s)" % (b[0], b[1]), p, "")

	def _revcomp(self, seq):
		"""Reverse-complement a sequence."""

		complstr = (
        	'GATCRYMKSWHBVDNXgatcrymkswhbvdnx-.',
        	'CTAGYRKMSWDVBHNXctagyrkmswdvbhnx-.'
        	)
		table = str.maketrans(*complstr)
		seq = str.translate(seq, table)
		return reduce(lambda a, b: b + a, seq, '')

	def _revcomp_m(self, coor1, coor2):
		"""Calculate coordinates on the reverted strand."""

		return (self.seqlen - coor2, self.seqlen - coor1)

	# PUBLIC METHODS

	def __len__(self):
		"""Amount of matches."""
	
		return len(self.mm)
	
	def set_pattern(self, p=[]):
		"""Definition of a search pattern and result groups."""

		if p:
			self.p = self._construct_pattern(p)
			self.groups = [v[0] for v in p]
		else:
			self.p = None
			self.groups = None


	def set_sequence(self, name='', seq='', offset=0):
		"""Setting of several attributes."""

		self.name = name
		self.seq = seq
		self.rseq = ''
		self.seqlen = len(seq)
		self.mm = []
		self.status = 0
		self.offset = offset

	#def findall(self, mm=[], offset=0, direct='F'):
	def findall(self, offset=0, direct='F'):
		"""
		Search all matches of the defined pattern and return their coordinates.
		This method 
		"""

		if direct == 'R':
			if not self.rseq:
				self.rseq = self._revcomp(self.seq)
		elif direct != 'F':
			raise ValueError("wrong argument for direction: '%s'" % direct)

		if direct == 'F':
			m = re.search(self.p, self.seq[offset:], re.I)
		elif direct == 'R':
			m = re.search(self.p, self.rseq[offset:], re.I)
		if m:
			start_coor = offset + m.start(1)
			end_coor = start_coor
			d = {
				'direct': direct, 
				'start': start_coor,
				}
			for g in self.groups:
				d["%s_start" % g] = offset + m.start(g)
				d["%s_end" % g] = offset + m.end(g)
				d["%s_seq" % g] = m.group(g)
				end_coor = max((end_coor, d["%s_end" % g]))
			d['end'] = end_coor
			if direct == 'F':
				d['seq'] = self.seq[start_coor:end_coor]
			elif direct == 'R':
				d['seq'] = self.rseq[start_coor:end_coor]
			offset = end_coor + 1
			self.mm.append(d)
			self.findall(offset, direct)
		self.status |= 1
	
	def search(self, direct='F'):
		"""A wrapper for the findall method, to perform search on both strands."""

		self.mm = []
		if 'F' in direct:
			if self.cfg['verbose'] == 1:
				print("   searching forward")
			self.findall(0, 'F')
		if 'R' in direct:
			if self.cfg['verbose'] == 1:
				print("   searching reverse")
			self.findall(0, 'R')

	def _prepare_matches_for_output(self, groups, formattype, as_motif):

		if not self.header_written:
			mm = [self.__class__.prepend[formattype]]
			self.header_written = True
		for m in self.mm:
			args = {} 
			if groups:
				if m['direct'] == 'F':
					seq = self.seq.lower()
				else:
					seq = self._revcomp(self.seq).lower()
				for g in groups:
					start, end = m['%s_start' % g], m['%s_end' % g]
					seq = seq[:start] + m['%s_seq' % g].upper() + seq[end:]
				args['seq'] = seq[m['start']:m['end']]
			else:
				args['seq'] = m['seq']
			args['direct'] = m['direct']
			if m['direct'] == 'F':
				args['start'] = m['start']+self.offset
				args['end'] = m['end']+self.offset
			else:
				args['start'] = (self.offset+self.seqlen)-m['end']
				args['end'] = (self.offset+self.seqlen)-m['start']
			args['name'] = "%s_%s-%s" % (self.name, args['start'], args['end'])
			mm.append(self.__class__.outformat[formattype] % args)
		return ''.join(mm)

	def matches_as_fasta(self, groups=[], as_motif=False):
		"""Output matches, coordinates and other data in FASTA format."""
		
		return self._prepare_matches_for_output(groups, 'fasta', as_motif)
	
	def matches_as_csv(self, groups=[], as_motif=False):
		"""Output matches, coordinates and other data in FASTA format."""
		
		return self._prepare_matches_for_output(groups, 'csv', as_motif)
	
	# DEBUGGING FUNCTIONS

	def show_matches(self, groups=[], wrap=100):
		"""
		Show matches aligned to the source sequence.
		You have to specify at least one group, otherwise
		nothing will be shown!
		
		TO DO:
		- one line for forward matches and one for reverse!
		- display as index not as character
		"""

		mm = [' '] * self.seqlen
		for m in self.mm:
			for g in groups:
				if m['direct'] == 'F':
					start, end = m['%s_start' % g], m['%s_end' % g]
				else:
					start, end = self._revcomp_m(m['%s_start' % g], m['%s_end' % g])
				mm[start:end] = (end - start) * [g[0]]
		out = []
		mm = ''.join(mm)
		for i in range(0, self.seqlen, wrap):
			out.append(("% 7d  % -"+str(wrap)+"s  % 7d\n") % (i+1+self.offset, self.seq[i:i+wrap], i+wrap+self.offset))
			out.append(("% 7d  % -"+str(wrap)+"s  % 7d\n") % (i+1+self.offset, mm[i:i+wrap], i+wrap+self.offset))
			out.append("\n")
		print("".join(out))

	def print_mm(self):
		"""Print data (strand, coordinates, sequence) of all matches."""

		for m in self.mm:
			print('direct', m['direct'])
			print('start', m['start']+self.offset)
			print('end', m['end']+self.offset)
			print('seq', m['seq'])
			for g in self.groups:
				for t in ('start', 'end', 'seq'):
					print("%s_%s" % (g, t), m["%s_%s" % (g, t)])
			print("-" * 20)


class SINEFinder(MotifFinder):

	config = {
		'verbose': 0,
		'TSD_MISMATCH_TOLERANCE': 2,
		'TSD_MIN_WORDSIZE': 5,
		'TSD_MISMATCH_PENALTY': 1,
		'TSD_SCORE_CUTOFF': 10,
		'TSD_ORIENTATION': 'F', # possible: 'FR'
		}

	pattern = (
		('TSD_region_1', ".{,40}"), 
		('a_box', "[GA][CGA]TGG"),
		('spacer_1', ".{25,50}"), 
		('b_box', "GTTC[AG]A"), 
		('spacer_2', ".{20,500}?"), 
		('polyA', "A{6,}|T{6,}"), 
		('TSD_region_2', ".{,40}"), 
		)

	outformat = {
		'fasta': ">%(name)s %(direct)s %(start)s:%(end)s TSD-len=%(length)s;TSD-score=%(score)s;TSD-mism=%(mm_count)s\n%(seq)s\n",
		'csv': "%(name)s\t%(direct)s\t%(start)s\t%(end)s\t%(length)s\t%(score)s\t%(mm_count)s\t%(seq)s\n",
		}
	prepend = {
		'fasta': "",
		'csv': "name\tdirect\tstart\tend\tTSD-len\tTSD-score\tTSD-mism\tsequence\n",
		}

	def __len__(self):
		"""Amount of matches with TSD-status True. Overrides method of base class."""
		
		if self.status & 2:
			return len([d for d in self.mm if d['TSD_status'] == True])
		else:
			return 0

	def evaluate_TSD(self):
		"""Extension of the base class. Evaluating TSD_regions for TSDs. Must be performed after findall-method."""

		for i in range(len(self.mm)):
			tsd = self.get_best_match(self.mm[i]['TSD_region_1_seq'], self.mm[i]['TSD_region_2_seq'])
			if tsd and tsd[3] >= self.cfg['TSD_SCORE_CUTOFF']:
				# TSD found
				start_coor = self.mm[i]['TSD_region_1_start'] + tsd[0][0]
				end_coor = self.mm[i]['TSD_region_2_start'] + tsd[1][1] + 1
				self.mm[i]['start'] = start_coor
				self.mm[i]['end'] = end_coor
				self.mm[i]['seq'] = self.seq[start_coor:end_coor]
				self.mm[i]['TSD_1_start'] = start_coor
				self.mm[i]['TSD_1_end'] = self.mm[i]['TSD_region_1_start'] + tsd[0][1] + 1
				self.mm[i]['TSD_2_start'] = self.mm[i]['TSD_region_2_start'] + tsd[1][0]
				self.mm[i]['TSD_2_end'] = end_coor
				self.mm[i]['TSD_score'] = tsd[3]
				self.mm[i]['TSD_mismatches'] = tsd[2]
				self.mm[i]['TSD_ORIENTATION'] = tsd[4]
				if self.mm[i]['direct'] == 'F':
					self.mm[i]['TSD_1_seq'] = self.seq[self.mm[i]['TSD_1_start']:self.mm[i]['TSD_1_end']]
					self.mm[i]['TSD_2_seq'] = self.seq[self.mm[i]['TSD_2_start']:self.mm[i]['TSD_2_end']]
				elif self.mm[i]['direct'] == 'R':
					if not self.rseq:
						self.rseq = self._revcomp(self.seq)
					self.mm[i]['TSD_1_seq'] = self.rseq[self.mm[i]['TSD_1_start']:self.mm[i]['TSD_1_end']]
					self.mm[i]['TSD_2_seq'] = self.rseq[self.mm[i]['TSD_2_start']:self.mm[i]['TSD_2_end']]
				self.mm[i]['TSD_status'] = True
			else:
				self.mm[i]['TSD_status'] = False
		self.status |= 2

	def _prepare_matches_for_output(self, groups, formattype, as_motif):

		mm = []
		if self.__class__.prepend[formattype] and not self.header_written:
			mm.append(self.__class__.prepend[formattype])
			self.header_written = True
		args = {'name': self.name} 
		for m in self.mm:
			args['direct'] = m['direct']
			if m['direct'] == 'F':
				args['start'] = m['start']+self.offset
				args['end'] = m['end']+self.offset
			else:
				args['start'] = (self.offset+self.seqlen)-m['end']
				args['end'] = (self.offset+self.seqlen)-m['start']
			try:
				if m['TSD_status'] == True:
					if groups:
						if m['direct'] == 'F':
							seq = self.seq.lower()
						else:
							seq = self._revcomp(self.seq).lower()
						for g in groups:
							start, end = m['%s_start' % g], m['%s_end' % g]
							seq = seq[:start] + m['%s_seq' % g].upper() + seq[end:]
						args['seq'] = seq[m['start']:m['end']]
					else:
						args['seq'] = m['seq']
					args['length'] = len(m['TSD_1_seq'])
					args['score'] = m['TSD_score']
					args['mm_count'] = m['TSD_mismatches']
					mm.append(self.__class__.outformat[formattype] % args)
					if self.cfg['verbose'] == 1:
						print("   writing %(name)s %(direct)s %(start)s:%(end)s: TSD found!" % args)
				else:
					if self.cfg['verbose'] == 1:
						print("   skipping %(name)s %(direct)s %(start)s:%(end)s: no TSD!" % args)
					continue
			except KeyError:
				if self.cfg['verbose'] == 1:
					print("   skipping %(name)s %(direct)s %(start)s:%(end)s: ERROR!" % args)
		return ''.join(mm)

	def pattern_matches(self, groups=[]):
		"""Output matches, coordinates and other data in FASTA format."""
	
		mm = []
		args = {'name': self.name} 
		for m in self.mm:
			if groups:
				if m['direct'] == 'F':
					seq = self.seq.lower()
				else:
					seq = self._revcomp(self.seq).lower()
				for g in groups:
					start, end = m['%s_start' % g], m['%s_end' % g]
					seq = seq[:start] + m['%s_seq' % g].upper() + seq[end:]
				#if as_motif:
				#	for k in m:
				#		km = re.match("(spacer_\d)_seq", k)
				#		if km:
				#			start, end = m['%s_start' % km.group(1)], m['%s_end' % km.group(1)]
				#			seq = seq[:start] + len(m['%s_seq' % km.group(1)]) * 'N' + seq[end:]
				args['seq'] = seq[m['start']:m['end']]
			else:
				args['seq'] = m['seq']
			args['direct'] = m['direct']
			#args['length'] = len(m['TSD_1_seq'])
			#args['score'] = m['TSD_score']
			#args['mm_count'] = m['TSD_mismatches']
			if m['direct'] == 'F':
				args['start'] = m['start']+self.offset
				args['end'] = m['end']+self.offset
			else:
				args['start'] = (self.offset+self.seqlen)-m['end']
				args['end'] = (self.offset+self.seqlen)-m['start']
			args['name'] += "_%s-%s" % (args['start'], args['end'])
			mm.append(self.__class__.outformat[formattype] % args)
		return ''.join(mm)
	
	# METHODS FOR TSD RECOGNITION

	def get_best_match(self, seq1, seq2):
		"""Search in two sequences for similar subsequences. Mismatches up to
		the value of MISMATCH_TOLERANCE are accepted."""
	
		rseqlen = len(seq2)
		mm = ()
		if 'F' in self.cfg['TSD_ORIENTATION']:
			offset = 0
			while 1:
				m = self.get_seed(seq1, seq2, offset)
				if not m:
					break
				em = self.extend_match(seq1, seq2, m)
				if not mm or em[3] > mm[3]:
					mm = em
				offset = em[0][1]
		if 'R' in self.cfg['TSD_ORIENTATION']:
			seq3 = self._revcomp(seq2)
			offset = 0
			while 1:
				m = self.get_seed(seq1, seq3, offset)
				if not m:
					break
				em = self.extend_match(seq1, seq3, m, 'R')
				if not mm or em[3] > mm[3]:
					s = rseqlen - em[1][1];
					em[1][1] = rseqlen - em[1][0];
					em[1][0] = s;
					mm = em
				offset = em[0][1]
		return mm

	def get_seed(self, seq1, seq2, offset=0):
		"""A window of the size of MIN_WORDSIZE shifts over sequence 1 and
		a search for all exact matches of this subsequence in sequence 2 is
		performed. All matches and their coordinates are returned. These
		serve as seed for extension."""
	
		seqlen = len(seq1)
		for ws in range(offset, seqlen-self.cfg['TSD_MIN_WORDSIZE'], 1):
			p = seq1[ws:min(ws+self.cfg['TSD_MIN_WORDSIZE'], seqlen)]
			i = seq2.find(p)
			if i > -1:
				return (p, ws, i)
		return None
	
	def extend_match(self, seq1, seq2, m, typ='F'):
		"""Try to extend given matches of size MIN_WORDSIZE to the right 
		and to the left until MISMATCH_TOLERANCE is reached. Terminal
		mismatches will be clipped."""
	
		mismatches = 0
		end_mm = 0
		start_mm = 0
		is_end = 0

		# The seed

		m1 = [m[1], m[1]+self.cfg['TSD_MIN_WORDSIZE']-1]
		m2 = [m[2], m[2]+self.cfg['TSD_MIN_WORDSIZE']-1]
		
		# Now the extension

		while 1:
			#if mismatches > MISMATCH_TOLERANCE:
			#	raise StandardError("SystematicError: mismatches exceed MISMATCH_TOLERANCE!")
			if not is_end & 1:
				m1[1] += 1
				m2[1] += 1
				if m1[1] < len(seq1) and m2[1] < len(seq2):

					# Not at right border
					if seq1[m1[1]].upper() == seq2[m2[1]].upper():

						# Position is equal
						if end_mm:

							# The first fitting base after 1 or more mismatches
							mismatches += end_mm
							end_mm = 0

					else:

						# Position is not equal
						# add a mismatch for this direction
						end_mm +=1

						if mismatches + end_mm > self.cfg['TSD_MISMATCH_TOLERANCE']:

							# This mismatch end exceeds tolerance:
							# end extension and forget about the last pos
							is_end |= 1
							m1[1] -= end_mm
							m2[1] -= end_mm					
							end_mm = 0

				else:

					# Reached the right border.
					is_end |= 1

					# Clean up: 
					# if it ends with mismatches, discard them.
					if end_mm:
						m1[1] -= end_mm
						m2[1] -= end_mm
						end_mm = 0
					else:
						m1[1] -= 1
						m2[1] -= 1

			if not is_end & 2:
				m1[0] -= 1
				m2[0] -= 1
				if m1[0] >= 0 and m2[0] >= 0:

					# Not at left border
					if seq1[m1[0]].upper() == seq2[m2[0]].upper():

						# Position is equal
						if start_mm:

							# The first fitting base after 1 or more mismatches
							mismatches += start_mm
							start_mm = 0

					else:

						# Position is not equal:
						# add a mismatch for this direction
						start_mm +=1

						if mismatches + start_mm > self.cfg['TSD_MISMATCH_TOLERANCE']:

							# This mismatch end exceeds tolerance:
							# end this extension and forget about the last pos.
							is_end |= 2
							m1[0] += start_mm
							m2[0] += start_mm
							start_mm = 0

				else:

					# Reached the right border.
					is_end |= 2

					# Clean up: 
					# if it ends with mismatches, discard them.
					if start_mm:
						m1[0] += start_mm
						m2[0] += start_mm
						start_mm = 0
					else:
						m1[0] += 1
						m2[0] += 1

			if is_end == 3 or mismatches == self.cfg['TSD_MISMATCH_TOLERANCE']:
				break 

		mismatches -= end_mm + start_mm
		m1[0] += start_mm
		m1[1] -= end_mm
		m2[0] += start_mm
		m2[1] -= end_mm
		score = m1[1] - m1[0] + 1 - mismatches*self.cfg['TSD_MISMATCH_PENALTY']
		return (m1, m2, mismatches, score, typ)


###########
# HELPERS #
###########

def run(seqfile, **kwargs):
	"""The run-function is responsible for the analysis."""

	if kwargs.setdefault('verbose', 0):
		print("Run started: %s" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	# check filetype (fasta or csv)

	if not seems_to_be_fasta(seqfile):
		exit_on_error(3, seqfile)

	# Creating the instances, opening the file handle(s) for the result file(s)
	# and defining the variables.
	
	fi = FastaIterator(seqfile, **kwargs)
	sf = SINEFinder(**kwargs)
	
	matchfile = '.'.join(seqfile.split('.')[:-1]) + "-matches"
	fw = {}
	if kwargs['OUTTYPE'] in ('fasta', 'both'):
		fw['fasta'] = open("%s.fasta" % matchfile, 'a')
	if kwargs['OUTTYPE'] in ('csv', 'both'):
		fw['csv'] = open("%s.csv" % matchfile, 'a')

	i, id, seq, offset = 0, '', '', 0

	# start processing

	while 1:
		e = next(fi)
		if e:
			(id, seq, offset) = e
			if kwargs['verbose'] == 1:
				print("\n   processing (%s) '%s', segment %s:%s" % (i, id, offset, offset+len(seq)), 'SEQLEN=', len(seq))
			sf.set_sequence(id, seq, offset)
			sf.search(kwargs['TSD_ORIENTATION'])
			sf.evaluate_TSD()
			if kwargs['OUTTYPE'] in ('fasta', 'both'):
				fw['fasta'].write(sf.matches_as_fasta(['TSD_1', 'a_box', 'b_box', 'polyA', 'TSD_2']))
			if kwargs['OUTTYPE'] in ('csv', 'both'):
				fw['csv'].write(sf.matches_as_csv(['TSD_1', 'a_box', 'b_box', 'polyA', 'TSD_2']))
			i += 1
		else:
			break

	# clean up

	for k in fw:
		fw[k].close()

	if kwargs['verbose'] == 1:
		cfgi = list(options.items())
		cfgi.sort()
		cfgt = "\n".join(["   % -30s = %s" % i for i in cfgi])
		print("\nSINEFinder-run with config: \n%s" % cfgt)
		print("Run terminated: %s" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


##############
#### MAIN ####
##############

def main():
	"""
	The main-function:
	 - utilizing commandline arguments or prompting for arguments in interactive mode
	 - and calling the run-function for analysis.
	"""

	try:

		# Help needed?

		if '-d' in sys.argv:
			sys.stdout.write(__description__)
			sys.exit()

		if '-h' in sys.argv or '-help' in sys.argv:
			sys.stdout.write(__usage__ % options)
			sys.exit()

		if '-v' in sys.argv:
			sys.stdout.write(__version__ % (__version_no__, __version_date__))
			sys.exit()

		# Under which mode should sine_finder run?

		if len(sys.argv) == 1: 

			# No arguments = interactive mode

			title = " %s (version: %s) " % (__program_name__, __version_no__)
			section = "\n*** %s ***"
			subsection = "\n%s"
			argument = "% -30s = %s"

			print()
			print("%s%s%s" % ('#'*5, '#'*len(title), '#'*5))
			print("%s%s%s" % ('#'*5, title, '#'*5))
			print("%s%s%s" % ('#'*5, '#'*len(title), '#'*5))
			print()
			print("Quit program at any prompt with %s" % (" or ".join(list(EXITCHAR))))

			while 1:

				# - Get infiles by browsing through directories.
				#   Only files with supported extensions for fasta files are shown.
	
				while 1:
	
					ff = [f for f in os.listdir('.') if os.path.isfile(os.path.join('.', f)) and not f.startswith('.') and filter_extension(f)]
					dd = [f for f in os.listdir('.') if os.path.isdir(os.path.join('.', f))]
					ddi = ["CHANGE DIR TO..." + d for d in dd]
					files = [('..', 'UP ONE DIR')] + list(zip(dd, ddi)) + list(zip(ff, ff))
					
					print(section % "SEQUENCE FILE")
					print(subsection % ("Fasta files in directory \n\t'%s'" % os.path.abspath('.')))
					fname = prompter(
						'Select the file to be analysed',
						options=files, 
						input_required=True
						)
	
					if os.path.isdir(fname):
						os.chdir(fname)
					else:
						infile = os.path.basename(fname)
						filesize = os.path.getsize(fname)
						break
	
				# - Selection of the way files should be read and processed

				runtypes = [
					('chunkwise', 'chunkwise (process subsequences of defined size)'), 
					('seqwise', 'seqwise (search sequence by sequence)')
					]

				print(section % 'PROCESSING TYPE')
				options['RUNTYPE'] = prompter(
					'Select a type', 
					options=runtypes, 
					input_required=True, 
					default=0, 
					note="default: %s" % runtypes[0][0]
					)

				# For chunkwise processing of files prompting for some parameters.
	
				if options['RUNTYPE'] == 'chunkwise':
				
					options['CHUNKSIZE'] = int(prompter(
						'Define chunk size', 
						default=options['CHUNKSIZE'], 
						validate=is_integer, 
						note="default:%s" % options['CHUNKSIZE'])
						)
					options['OVERLAP'] = int(prompter(
						'Define size of overlap', 
						default=options['OVERLAP'], 
						validate=lambda v: is_integer(v) and int(v) <= options['CHUNKSIZE'], 
						note="default:%s" % options['OVERLAP'])
						)
	
				# Parameters related to pattern search
	
				print(section % 'TSD SEARCH OPTIONS')
				
				options['TSD_MIN_WORDSIZE'] = int(prompter(
					'Minimal wordsize of TSD seed', 
					default=options['TSD_MIN_WORDSIZE'], 
					validate=is_integer, 
					note="default:%s" % options['TSD_MIN_WORDSIZE'])
					)
				options['TSD_MISMATCH_TOLERANCE'] = int(prompter(
					'TSD mismatch tolerance', 
					default=options['TSD_MISMATCH_TOLERANCE'], 
					validate=is_integer, 
					note="default:%s" % options['TSD_MISMATCH_TOLERANCE'])
					)
				options['TSD_MISMATCH_PENALTY'] = int(prompter(
					'TSD mismatch penalty', 
					default=options['TSD_MISMATCH_PENALTY'], 
					validate=is_integer, 
					note="default:%s" % options['TSD_MISMATCH_PENALTY'])
					)
				options['TSD_SCORE_CUTOFF'] = float(prompter(
					'TSD score cutoff', 
					default=options['TSD_SCORE_CUTOFF'], 
					validate=is_float, 
					note="default:%s" % options['TSD_SCORE_CUTOFF'])
					)
				orientations = (('F', 'forward'),('R', 'reverse'),('FR', 'both directions'))
				options['TSD_ORIENTATION'] = prompter(
					'Direction of TSD search', 
					default={'F':'0','R':'1','FR':'2'}[options['TSD_ORIENTATION']], 
					options=orientations, 
					note="default:%s" % {'F':'0','R':'1','FR':'2'}[options['TSD_ORIENTATION']]
					)

				# Other parameters

				verbosity = ((0, 'no'), (1, 'yes'))
				options['verbose'] = int(prompter(
					'Verbose', 
					default=options['verbose'], 
					options=verbosity, 
					note="default:%s" % options['verbose'])
					)
				filetypes = (('fasta', 'fasta file'),('csv', 'tab-delimited file'),('both', 'both types'))
				options['OUTTYPE'] = prompter(
					'Type of result file', 
					default=options['OUTTYPE'], 
					options=filetypes, 
					note="default:%s" % {'fasta':'0','csv':'1','both':'2'}[options['OUTTYPE']]
					)

				# Confirmation
	
				print(section % 'JOB CONFIRMATION')

				for arg in (
					("File", infile),
					("File size", "%s B" % filesize),
					("Processing type", options['RUNTYPE']),
					("Chunk size", ["-", "%s bases" % options['CHUNKSIZE']][options['RUNTYPE'] == 'chunkwise']),
					("Overlap", ["-", "%s bases" % options['OVERLAP']][options['RUNTYPE'] == 'chunkwise']),
					("Minimal wordsize of TSD seed", options['TSD_MIN_WORDSIZE']),
					("TSD mismatch tolerance", options['TSD_MISMATCH_TOLERANCE']),
					("TSD mismatch penalty", options['TSD_MISMATCH_PENALTY']),
					("TSD score cutoff", options['TSD_SCORE_CUTOFF']),
					("Direction of TSD search", options['TSD_ORIENTATION']),
					("Writing result to", '.'.join(infile.split('.')[:-1]) + "-matches + extension(s)"),
					("Type of result file", options['OUTTYPE']),
					("Verbose", ['no', 'yes'][options['verbose']]),
					):
					print(argument % arg) 
				print()
				choice = input("Start job with these arguments: y(es), n(o) or q(uit)? ")
				if choice in 'yY':
					break
				elif choice in 'nN':
					continue
				elif choice in 'qQ':
					sys.exit("\nGood bye!")
				continue

		else:

			# Arguments specified = commandline mode

			# - Get infiles

			if len(sys.argv) < 2 or sys.argv[-1].startswith("-"):
				exit_on_error(13)
	
			if os.path.isfile(sys.argv[-1]):
				infile = sys.argv[-1]
			else:
				exit_on_error(1, sys.argv[-1])

			# - The way infiles should be processed.

			if '-T' in sys.argv:
				try:
					if sys.argv[sys.argv.index('-T')+1] in ('seqwise', 'chunkwise'):
						options['RUNTYPE'] = sys.argv[sys.argv.index('-T')+1]
					else:
						exit_on_error(51, '-T', "in ('seqwise', 'chunkwise')")
				except (KeyError, ValueError):
					exit_on_error(2, '-T')
	
				if options['RUNTYPE'] == 'chunkwise': 
	
					if '-C' in sys.argv:
						try:
							if is_integer_argv('-C'):
								options['CHUNKSIZE'] = int(sys.argv[sys.argv.index('-T')+1])
							else:
								exit_on_error(51, '-C', 'an integer value')
						except (KeyError, ValueError):
							exit_on_error(2, '-C')
		
					if '-O' in sys.argv:
						try:
							if is_integer_argv('-O'):						
								options['OVERLAP'] = int(sys.argv[sys.argv.index('-O')+1])
								if options['OVERLAP'] <= options['CHUNKSIZE']:
									exit_on_error(53)
							else:
								exit_on_error(51, '-O', 'an integer value')
						except (KeyError, ValueError):
							exit_on_error(2, '-O')
		
			# - arguments for searching TSDs
		
			if '-t' in sys.argv:
				try:
					if is_integer_argv('-t'):
						options['TSD_MISMATCH_TOLERANCE'] = int(sys.argv[sys.argv.index('-t')+1])
					else:
						exit_on_error(51, '-t', 'an integer value')
				except (KeyError, ValueError):
					exit_on_error(2, '-t')
		
			if '-w' in sys.argv:
				try:
					if is_integer_argv('-w'):
						options['TSD_MIN_WORDSIZE'] = int(sys.argv[sys.argv.index('-w')+1])
					else:
						exit_on_error(51, '-w', 'an integer value')
				except (KeyError, ValueError):
					exit_on_error(2, '-w')
		
			if '-p' in sys.argv:
				try:
					if is_integer_argv('-p'):
						options['TSD_MISMATCH_PENALTY'] = int(sys.argv[sys.argv.index('-p')+1])
					else:
						exit_on_error(51, '-p', 'an integer value')
				except (KeyError, ValueError):
					exit_on_error(2, '-p')
		
			if '-s' in sys.argv:
				try:
					if is_float_argv('-s'):
						options['TSD_SCORE_CUTOFF'] = float(sys.argv[sys.argv.index('-s')+1])
					else:
						exit_on_error(51, '-s', 'a float value')
				except (KeyError, ValueError):
					exit_on_error(2, '-s')
		
			if '-o' in sys.argv:
				try:
					if sys.argv[sys.argv.index('-o')+1] in 'FR':
						options['TSD_ORIENTATION'] = int(sys.argv[sys.argv.index('-o')+1])
					else:
						exit_on_error(51, '-o', "in ('F', 'R', 'FR')")
				except (KeyError, ValueError):
					exit_on_error(2, '-o')
	
			# - Other options
	
			if '-f' in sys.argv:
				try:
					if sys.argv[sys.argv.index('-f')+1] in ('fasta', 'csv', 'both'):
						options['OUTTYPE'] = sys.argv[sys.argv.index('-f')+1]
					else:
						exit_on_error(51, '-f', "in ('fasta', 'csv', 'both')")
				except (KeyError, ValueError):
					exit_on_error(2, '-f')
	
			if '-V' in sys.argv:
				options['verbose'] = 1
				print(options)
	
		run(infile, **options)

	except Exception as py_msg:

		# catch unrecognized errors and write pythons error message to
		# STDERR
		sys.stderr.write("%s\n" % py_msg)
		exit_on_error(0)


### main ###

if __name__ == "__main__":
	main()

