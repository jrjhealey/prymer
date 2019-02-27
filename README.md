# prymer

Generate primers for input DNA sequences in an automated manner.

# Installation

    $ pip install prymer

Primer's only (running) dependency is `biopython`.

# Usage

The only required arguments are positionals for the input sequence (in fasta or genbank format),
and a file to write the output sequences to.

    $ prymer.py sequences outfile

## Full usage options

    usage: prymer.py [-h] [-v] [-m {bracketed,sanger}] [-t TILE] [-o OFFSET]
                     [-l, LENGTH] [-f FORMAT]
                     [-s {25nm,100nm,250nm,1um,5um,10um,4nmU,20nmU,PU,25nmS}]
                     [-p {STD,PAGE,HPLC,IEHPLC,RNASE,DUALHPLC,PAGEHPLC}]
                     sequences outfile

    Return primers for given input sequences.

    positional arguments:
      sequences             The sequence file to design primers for. It is assumed
                            all sequences are provided 5' -> 3'
      outfile               Output file of primers in the chosen format.

    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         Print additional information about the designed
                            sequences (to STDERR).

    Primer Customisation Options:
      -m {bracketed,sanger}, --method {bracketed,sanger}
                            The design scheme for the primers (i.e. how to design
                            them).Current options: bracketed - Designs 'normal' 5'
                            and 3' amplification primers. sanger - Designs 'primer
                            walking' primers, approximately 800bp apart
      -t TILE, --tile TILE  Use this option to overwrite the default between-
                            primer walking distance to increase the sequencing
                            overlap (e.g. for poor sequencing results. Sanger
                            sequencing typically returns about ~1kb of usable
                            sequence, but the 3' end will often begin to become
                            error-laden.
      -o OFFSET, --offset OFFSET
                            A small offset from the start of the sequence, to
                            manipulate primer position. This helps increase 5'
                            junction primer coverage.
      -l, LENGTH, --length LENGTH
                            Desired length of primer sequences.

    Output Options:
      -f FORMAT, --format FORMAT
                            What file type to output the primers as (fasta,
                            separated by delimiter of choice, or IDT bulk upload
                            format).Delimiter should be specified in quotes, e.g.
                            ',' or ';'.
      -s {25nm,100nm,250nm,1um,5um,10um,4nmU,20nmU,PU,25nmS}, --synthesis {25nm,100nm,250nm,1um,5um,10um,4nmU,20nmU,PU,25nmS}
                            Specify synthesis scale options to add to IDT output
                            (e.g. "25nm")
      -p {STD,PAGE,HPLC,IEHPLC,RNASE,DUALHPLC,PAGEHPLC}, --purification {STD,PAGE,HPLC,IEHPLC,RNASE,DUALHPLC,PAGEHPLC}
                            Purification type for IDT output.

    Easy peasy primer design. Recent updates now include the ability to return
    melting temperatures, as well as design overlapping/primer-walking schemes.

----
Currently, there are 2 main 'modes' of primer design: Bracketed and Sanger.

Bracketed primers are 'normal' amplification primers, matching the start and end
of a given sequence.

e.g, for the sequence below, bracket primers appear as:

    F-TGGGTAAAATAATTGGTATC->
      ATGGGTAAAATAATTGGTATCGACCTGGGTACTACCAACTCTTGTGTAGCGATTATGGATGGCACCACTCCTCGCGTGCTGGAGAACGCCGAA
                                                                             <-AAGCCGCAAGAGGTCGTGCG-R

The reverse primer is of course reverse complemented relative to the forward sequence.

Sanger primers provide 'tiled' primers for primer walking across a sequence (mainly intended
for sequencing applications):

e.g., for the sequence below, 10-mer primers, tiled every 30 bp would appear as (showing forward only):


    ATGGGTAAAA->                  ACTACCAACT->                  GGCACCACTC->
    ATGGGTAAAATAATTGGTATCGACCTGGGTACTACCAACTCTTGTGTAGCGATTATGGATGGCACCACTCCTCGCGTGCTGGAGAACGCCGAA

# Options

 - `-t` | `--tile`  - the separation distance between Sanger sequencing primers
 - `-o` | `--offset` - A 5' offset to 'nudge' all sequencing primers along the sequence (helps to span junctions etc.)
 - `-l` | `--length` - A target length of primer to return. In future, iterative primer pair refinement will have some control over this, but for now primers are a static length.
 - `-f` | `--format` - The output formats for the primers. By default the sequences will be written as fasta, but tabular formats with any custom separator can be specified. If `"IDT"` is specified, the primers will be output in a format for use with IDT, but you will also need to provide `-p` and `-s` (see below).
 - `-s` | `--synthesis` - A string corresponding to one of IDT's synthesis scales. This string will be appended to the output if `"IDT"` was chosen*.
 - `-p` | `--purification` - A string corresponding to one of IDT's purification methods*.

*Note, some purification methods and synthesis scales are mutually exclusive, but this script *DOES NOT CHECK*. Consult IDT's site.
