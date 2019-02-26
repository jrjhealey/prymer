#!/usr/bin/env python3

"""
Design of primer sequences easily from a provided sequence file.
"""
# TODO:
#  - Create a Tm calculator
#  - 'Gradient descent' style primer pair optimisation of Tms etc
#  - Self and hetero-dimer analysis (one day perhaps).

import sys
import argparse
import itertools

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from primer import Primer
from converter import convert_seqs


def get_args():
    """Parse command line arguments"""
    desc = """Return primers for given input sequences."""
    epi = """Easy peasy primer design. Recent updates now include
             the ability to return melting temperatures, as well as design
             overlapping/primer-walking schemes.
          """

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi, prog='prymer.py')
        parser.add_argument('-v', '--verbose', action='store_true',
                            help='Print additional information about the designed sequences (to STDERR')
        parser.add_argument('-l,', '--length', action='store', default=20, type=int,
                            help='Desired length of primer sequences.')
        parser.add_argument('-s', '--separator', action='store', default='fasta', type=str,
                            help='What file type to output the primers as (fasta, or separated by delimiter of choice.)'
                                 'Should be specified in quotes, e.g. \',\' or \';\'.')
        parser.add_argument('-m', '--method', action='store', default='bracketed', type=str,
                            choices=['bracketed', 'sanger'],
                            help='The design scheme for the primers (i.e. how to design them).'
                                 'Current options: bracketed - Designs \'normal\' 5\' and 3\' amplification primers. '
                                 'sanger - Designs \'primer walking\' primers, approximately 800bp apart')
        parser.add_argument('-t', '--tile', action='store', type=int, default=800,
                            help='Use this option to overwrite the default between-primer walking distance '
                                 'to increase the sequencing overlap (e.g. for poor sequencing results. '
                                 'Sanger sequencing typically returns about ~1kb of usable sequence, but '
                                 'the 3\' end will often begin to become error-laden.')
        parser.add_argument('-o', '--offset', action='store', type=int, default=100,
                            help='A small offset from the start of the sequence, to manipulate primer position. '
                                 'This helps increase 5\' junction primer coverage.')
        parser.add_argument('sequences', action='store',
                            help='The sequence file to design primers for. '
                                 'It is assumed all sequences are provided 5\' -> 3\' ')
        parser.add_argument('outfile', action='store', default=None,
                            help='Output file of primers in fasta format.')

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
        sys.exit(1)

    return parser.parse_args()


def main():
    """Easy primer design main function"""
    args = get_args()

    seqs = list(convert_seqs(args.sequences))

    Fprimers = []
    Rprimers = []
    if args.method == 'bracketed':
        for rec in seqs:
            # Create primers
            f = Primer(rec.id, rec.seq, length=args.length, direction='F')
            r = Primer(rec.id, rec.seq, length=args.length, direction='R')
            # Represent primers as SeqRecords to benefit from BioPythons output functions
            Fprimers.append(SeqRecord(f, id=f.name + "_" + f.direction, description=''))
            Rprimers.append(SeqRecord(r, id=r.name + "_" + r.direction, description=''))
    elif args.method == 'sanger':
        for rec in seqs:
            for i, subseq in enumerate([rec.seq[i:i+args.length] for i in range(args.offset, len(rec), args.tile)]):
                f = Primer(rec.id, subseq, length=args.length, direction='F')
                r = Primer(rec.id, subseq, length=args.length, direction='R')
                Fprimers.append(SeqRecord(f, id=f.name + "_" + str(i) + "_" + f.direction, description=''))
                Rprimers.append(SeqRecord(r, id=r.name + "_" + str(i) + "_" + r.direction, description=''))


    if args.verbose:
#            for a, b, c in zip(seqs, Fprimers, Rprimers):
        for a, b, c in itertools.zip_longest(seqs, Fprimers, Rprimers):
            if a is not None:
                sys.stderr.write("Target Sequence " + rec.id + " = " + "\n" + str(a.seq) + "\n")
            sys.stderr.write("Forward Primer (5\'-> 3\') = " + str(b.seq) + "\n")
            sys.stderr.write("Reverse Primer (5\'-> 3\') = " + str(c.seq) + "\n")

    if args.outfile:
        if not args.separator or args.separator == 'fasta':
            SeqIO.write((i for i in itertools.chain(*zip(Fprimers, Rprimers))), args.outfile, 'fasta')
        else:
            with open(args.outfile, 'w') as tfh:
                for i in itertools.chain(*zip(Fprimers, Rprimers)):
                    tfh.write(args.separator.join([str(i.id + i.description), str(i.seq) + "\n"]))

if __name__ == '__main__':
    main()
