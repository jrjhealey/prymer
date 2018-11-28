#!/usr/bin/env python

"""
Design of primer sequences easily from a provided sequence file.
"""


import sys
import argparse
from primer import Primer

try:
    from Bio import SeqIO
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
        parser = argparse.ArgumentParser(description=desc, epilog=epi, prog='PrYmer.py')
        parser.add_argument('-l,', '--length', action='store', default=20, type=int,
                             help='Desired length of primer sequences.')
        parser.add_argument('sequences', action='store',
                            help='The sequence file to design primers for. '
                                 'It is assumed all sequences are provided 5\' -> 3\' ')


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
    Fprimers = []
    Rprimers = []
    seqs = list(SeqIO.parse(args.sequences, 'fasta'))
    for rec in seqs:
        Fprimers.append(Primer(rec.id,rec.seq, length=args.length, direction='F'))
        Rprimers.append(Primer(rec.id, rec.seq, length=args.length, direction='R'))

    for a, b, c in zip(seqs, Fprimers, Rprimers):
        print("Target Sequence = " + str(a.seq))
        print("Forward Primer = " + str(b))
        print("Reverse Primer = " + str(c))



if __name__ == '__main__':
    main()
