# A converter function to allow multiple input sequence types
# Big thanks to Peter van Heusden for refactoring the genbank code to yeild generators!
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class InputFormatException(Exception):
    """The exception to return when an input format cannot be identified or coerced (informative name)"""
    def __init__(self, message):
        super().__init__(message)


def yield_from_genbank(infile):
    for rec in SeqIO.parse(infile, 'genbank'):
        cds_features = (f for f in rec.features if f.type == 'CDS')
        for feat in cds_features:
            try:
                header = feat.qualifiers['gene'][0] + '_' + rec.id
            except KeyError:
                header = feat.qualifiers['locus_tag'][0] + '_' + rec.id
            header.replace(' ', '_')
            yield SeqRecord(id=header, seq=feat.location.extract(rec).seq)


def convert_seqs(infile):
    with open(infile, 'r') as ifh:
        firstline = ifh.readline()
    if os.path.splitext(infile)[1] in (".fasta", ".fa", ".fas", ".fna"):
        try:
            assert firstline[0] == '>'
        except AssertionError:
            raise InputFormatException("File extension implies fasta but the first line doesn't look like a header.")
        return SeqIO.parse(infile, 'fasta')
    elif os.path.splitext(infile)[1] in (".gbk", ".gb", ".genbank", ".gbff"):
        try:
            assert firstline.startswith("LOCUS")
        except AssertionError:
            raise InputFormatException("File extension implies genbank, but the first line doesn't look like a header.")

        return yield_from_genbank(infile)
