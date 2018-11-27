from Bio.Seq import Seq
from Bio.Alphabet import generic_alphabet

class Primer(Seq):
    """A class to represent a primer sequence (Extends Bio.Seq.Seq).

        :name - The sequence name
        :data - The sequence to create primer from
        :length - A length to aim for (default 20
        :tm - A melting temperature to aim for. Primer length can be iterated until this is met.
        :direction - Is the primer "F"orward (5'->3'), "R"everse (5'->3') or other (def = None).
        :method - What 'style' of primer to design (inward or outward facing etc)
        :tmtype - What algorithm to use to calculate the melting temperature
                (To do: find a way to pass the optional Tm algorithm arguments for salt concs etc.
        :alphabet - The alphabet for the sequence (kept for Seq compatibility)

    """

    def __init__(self, name, data, length=int(20), tm=float(65), direction=None, method='simple',
                 tmtype='default', alphabet=generic_alphabet):
        """Create and store a primer with a name, sequence, and optionally a direction/alphabet.
           'method' is used as a switch to decide how to design the primer."""

        self._data = data # No idea why Seq decided to use this variable name
        self.name = name
        self.length = length
        self.tm = tm
        self.direction = direction
        self._method = method
        self._tmtype = tmtype
        self.alphabet = alphabet

        self.primerseq = self.create_sequence(self._data, self.direction, self.length, self._method)

#    def __repr__(self, name, _data):
#        return "'Primer("+self.name+":"+self._data+")'"

    def __str__(self):
        return str(self.__class__.__name__ + "(" + self.name + '_' + self.direction + ": " + str(self.primerseq) + ")" )

    def create_sequence(self, _data, direction, length, _method):
        """Return primer sequences for the input data
        """
        if _method == 'simple':     # Make bracketed sequence around provided sequence (this is the standard functionality)
            if direction == 'F':
                return _data[0:length]
            elif direction == 'R':
                return _data[-length:]

#                return super(Primer, self).reverse_complement(_data[-length:])



    def melting_temperature(self, primerseq, tmtype):
        """Calculate a melting temperature for the oligo according to different methods"""

        from Bio.SeqUtils import MeltingTemp as mt
        melt_method = {'wallace': mt.Tm_Wallace,
                       'default': mt.Tm_NN,
                       'gc': mt.Tm_GC,
                       'corrected': mt.Tm_NN}.get(tmtype, lambda: "No such method.")
        melt_method(primerseq)

