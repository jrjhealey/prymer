#!/usr/bin/env python

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
from primer import Primer
from converter import convert_seqs

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    msg = """
Could not import the BioPython module which means it is probably
not installed, or at least not available in the PYTHONPATH for this 
particular binary.

If you have conda (recommended) try running:

 $ conda install -c anaconda biopython

or, alternatively with python/pip:

 $ python -m pip install biopython
"""
    sys.stderr.write(msg)
    sys.exit(1)


def get_args():
    """Parse command line arguments"""
    desc = """Return primers for given input sequences."""
    epi = """Easy peasy primer design."""

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi, prog='prymer.py')
        parser.add_argument('-v', '--verbose', action='store_true',
                            help='Print additional information about the designed sequences (to STDERR')
        parser.add_argument('-l,', '--length', action='store', default=20, type=int,
                            help='Desired length of primer sequences.')
        parser.add_argument('-s', '--separator', action='store', default='fasta',
                            help='What file type to output the primers as (fasta, or separated by delimiter of choice.)'
                                 'Should be specified in quotes, e.g. \',\' or \';\'.')
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
    # If a tempfile was created by convert_seqs, remove it

    Fprimers = []
    Rprimers = []
    for rec in seqs:
        f = Primer(rec.id, rec.seq, length=args.length, direction='F')
        r = Primer(rec.id, rec.seq, length=args.length, direction='R')
        Fprimers.append(SeqRecord(f, id=f.name + "_" + f.direction, description=''))
        Rprimers.append(SeqRecord(r, id=r.name + "_" + r.direction, description=''))

    if args.verbose:
        for a, b, c in zip(seqs, Fprimers, Rprimers):
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
