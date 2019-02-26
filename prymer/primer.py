from Bio.Seq import Seq
from Bio.Alphabet import generic_alphabet
import sys

class Primer(Seq):
    """A class to represent a primer sequence (Extends Bio.Seq.Seq).

        :name - The sequence name
        :data - The sequence to create primer from
        :length - A length to aim for (default 20
        :tm - A melting temperature to aim for. Primer length can be iterated until this is met.
        :direction - Is the primer "F"orward (5'->3'), "R"everse (5'->3') or other (def = None).
        :alphabet - The alphabet for the sequence (kept for Seq compatibility)

    """

    def __init__(self, name, data, length=int(20), tm=float(65), direction=None,
                 method="simple", alphabet=generic_alphabet):
        """Create and store a primer with a name, sequence, and optionally a direction/alphabet.
           'method' is used as a switch to decide how to design the primer."""

        self._data = data # No idea why Seq decided to use this variable name
        self.name = name
        self.length = length
        self.tm = tm
        self.direction = direction
        self.method = method
        self.alphabet = alphabet

        assert len(data) > len(str(length)), "Primer sequences must be shorter than the sequence they target."
        try:
            direction
        except NameError:
            sys.stderr.write("Primers need a direction specified ('F'/'R').")
        try:
            name
        except NameError:
            sys.stderr.write("Primers need a base name string specified.")
        try:
            data
        except NameError:
            sys.stderr.write("Primers need a sequence to derive from pecified.")

        self.primerseq = self.create_sequence()

#    def __repr__(self, name, _data):
#        return "'Primer("+self.name+":"+self._data+")'"

    def __repr__(self):
        return str(self.__class__.__name__ + "(" + self.name + '_' + self.direction + ": " + str(self.primerseq) + ")" )

    def __str__(self):
        return str(self.primerseq)

    def create_sequence(self):
        """Return primer sequences for the input data
        """
        if self.method == 'simple':     # Make bracketed sequence around provided sequence (this is the standard functionality)
            if self.direction == 'F':
                return self._data[0:self.length]
            elif self.direction == 'R':
                return self._data[-self.length:].reverse_complement()

    #
    # def melting_temperature(self, primerseq, tmtype):
    #     """Calculate a melting temperature for the oligo according to different methods"""
    #
    #     from Bio.SeqUtils import MeltingTemp as mt
    #     melt_method = {'wallace': mt.Tm_Wallace,
    #                    'default': mt.Tm_NN,
    #                    'gc': mt.Tm_GC,
    #                    'corrected': mt.Tm_NN}.get(tmtype, lambda: "No such method.")
    #     melt_method(primerseq)
    #
