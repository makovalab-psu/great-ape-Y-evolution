import sys

from struct   import *
from UserDict import DictMixin
from sys      import stderr


class FastaFile(DictMixin):

	def __init__(self,f,unmask=False):
		self.debug = []

		# read magic and determine byte order

		if (type(f) == str):
			f = file(f,"rb")

		self.file   = f
		self.unmask = unmask

		# read the sequences and create a name-to-sequence hash; the same hash
		# also maps number-to-sequence; and if there is only one sequence, we
		# will map None to that sequence

		self.seqCount = 0
		self.index = {}
		self.indexByAlias = {}
		for (name,seq) in self.read_sequences(f):
			if (unmask): seq = seq.upper()
			self.index[name] = seq
			self.indexByAlias[self.seqCount] = seq
			self.seqCount += 1
		self.preloadedCount = self.seqCount
		if (self.seqCount == 1): self.indexByAlias[None] = self.indexByAlias[0]

	def __getitem__(self,name):
		if (name in self.index): return self.index[name]
		else:                    return self.indexByAlias[name]

	def keys(self):
		return self.index.keys()

	def read_sequences(self,f):
		seqName = None
		lineNum = 0
		for line in f:
			lineNum += 1
			line = line.strip()

			if (line.startswith(">")):
				if (seqName != None):
					seq = "".join(seq)
					yield (seqName,seq)
				fields = line[1:].split()
				if (fields == []):
					assert (False), \
					       "sequence has no name (at line %d)" % lineNum
				seqName = fields[0]
				seq = []
			elif (seqName == None):
				assert (False), "first sequence has no header"
			else:
				seq += [line]

		if (seqName != None):
			seq = "".join(seq)
			yield (seqName,seq)

