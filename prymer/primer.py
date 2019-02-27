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
        :method - This will eventually describe the manner in which to design the primers
        :alphabet - The alphabet for the sequence (kept for Seq compatibility)

    """

    def __init__(self, name, data, length=int(20), tm=float(65), direction=None,
                 alphabet=generic_alphabet):
        """Create and store a primer with a name, sequence, and optionally a direction/alphabet."""

        self._data = data  # No idea why Seq decided to use this variable name
        self.name = name
        self.length = length
        self.direction = direction
        self.alphabet = alphabet

        assert len(data) > len(str(length)), "Primer sequences must be shorter than the sequence they target."
        try:
            self.direction
        except NameError:
            sys.stderr.write("Primers need a direction specified ('F'/'R').")
        try:
            self.name
        except NameError:
            sys.stderr.write("Primers need a base name string specified.")
        try:
            self._data
        except NameError:
            sys.stderr.write("Primers need a sequence to derive from specified, and none was found.")

        self.primerseq = self.create_sequence()
        self.tm = self.melting_temperature()

    def __repr__(self):
        return str(self.__class__.__name__ + "(" + self.name + '_' + self.direction + ": " + str(self.primerseq) + ")")

    def __str__(self):
        return str(self.primerseq)

    def create_sequence(self):
        """Return primer sequences for the input data
        """
        if self.direction == 'F':
            return self._data[0:self.length]
        elif self.direction == 'R':
            return self._data[-self.length:].reverse_complement()

    def melting_temperature(self):
        """
        Calculate assorted melting temperatures for the oligo.
        Returns a dict of method: temp
        """

        from Bio.SeqUtils import MeltingTemp as mt
        temps = {'Wallace': str(round(mt.Tm_Wallace(self.primerseq), 2)),
                 'GC': str(round(mt.Tm_GC(self.primerseq), 2)),
                 'Corrected': str(round(mt.Tm_NN(self.primerseq), 2))
                 }
        temps['Average'] = str(round(
            ((float(temps['Wallace']) + float(temps['GC']) + float(temps['Corrected']))/3), 2))
        return temps


