# A converter function to allow multiple input sequence types

class InputFormatException(Exception):
    """The exception to return when an input format cannot be identified or coerced (informative name)"""
    def __init__(self, message):
        super().__init__(message)


def convert_seqs(infile):
    """Consume an input file of sequences, do basic tests on format and convert to multifasta"""
    from Bio import SeqIO
    from tempfile import NamedTemporaryFile
    import os

    temp = NamedTemporaryFile(mode='w')
    with open(infile, 'r') as ifh:
       firstline = ifh.readline()

    # Begin testing formats:
    if os.path.splitext(infile)[1] in (".fasta", ".fa", ".fas", ".fna"):
        # Give the benefit of the doubt by first testing the extension, but check that its well formed
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
    # Genbank will need features extracted if CDSs

        with open(temp.name, 'w') as tfh:
            for rec in SeqIO.parse(infile, 'genbank'):
                all_features = [feat for feat in rec.features if feat.type == "CDS"]
                for f in all_features:
                    try:
                        header = f.qualifiers['gene'][0]
                        header = header.replace(' ', '_')
                    except KeyError:
                        header = f.qualifiers['locus_tag']
                    tfh.write('>{}\n{}'.format(header, f.location.extract(rec).seq))

        tfh.seek(0)
        return SeqIO.parse(temp, 'fasta')

    return temp

