from unittest import TestCase
import sys, os
import prymer
from prymer.primer import Primer
from Bio import SeqIO

class TestBracketed(TestCase):
    def test_start_and_end_primers(self):
        rec = SeqIO.read(os.path.dirname(os.path.realpath(__file__))+"/dnaK.fasta", 'fasta')
        f = Primer(rec.id, rec.seq, length=20, direction='F')
        r = Primer(rec.id, rec.seq, length=20, direction='R')
        self.assertTrue(isinstance(f, Primer))
        self.assertTrue(isinstance(r, Primer))
        self.assertEqual(str(f), 'ATGGGTAAAATAATTGGTAT')
        self.assertEqual(str(r), 'TTATTTTTTGTCTTTGACTT')

class TestSanger(TestCase):
    def test_sanger_primers(self):
        length = 20
        offset = 100
        tile = 800
        rec = SeqIO.read(os.path.dirname(os.path.realpath(__file__))+"/dnaK.fasta", 'fasta')
        Fprimers = []
        Rprimers = []
        for x, subseq in enumerate([rec.seq[i:i + length] for i in range(offset, len(rec), tile)]):
            Fprimers.append(Primer(rec.id, subseq, length=length, direction='F'))
            Rprimers.append(Primer(rec.id, subseq, length=length, direction='R'))
        self.assertTrue(isinstance(Fprimers[0], Primer))
        self.assertTrue(isinstance(Rprimers[0], Primer))
        self.assertEqual(len(Fprimers), 3)
        self.assertEqual(len(Rprimers), 3)
        for a, b in zip(Fprimers, ['GCACCACGCCTTCTATCATT', 'ACTCGTGCGAAACTGGAAAG', 'TGCGCTGACTGCACTGGAAA']):
            self.assertEqual(str(a), b)
        for a, b in zip(Rprimers, ['AATGATAGAAGGCGTGGTGC', 'CTTTCCAGTTTCGCACGAGT', 'TTTCCAGTGCAGTCAGCGCA']):
            self.assertEqual(str(a), b)


