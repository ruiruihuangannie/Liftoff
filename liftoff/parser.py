import argparse

def parse_args(arglist):
    parser = argparse.ArgumentParser(description='Lift features from one genome assembly to another')
    parser.add_argument('target', help='target fasta genome to lift genes to')
    parser.add_argument('reference', help='reference fasta genome to lift genes from')

    refrgrp = parser.add_argument_group('Required input (annotation)')
    mxgrp = refrgrp.add_mutually_exclusive_group(required=True)
    mxgrp.add_argument(
        '-g', metavar='GFF', help='annotation file to lift over in GFF or GTF format'
    )
    mxgrp.add_argument(
        '-db', metavar='DB', help='name of feature database; if not specified, the -g '
                                  'argument must be provided and a database will be built automatically'
    )

    outgrp = parser.add_argument_group('Output')
    outgrp.add_argument(
        '-o', default='stdout', metavar='FILE',
        help='write output to FILE in same format as input; by default, output is written to terminal (stdout)'
    )
    outgrp.add_argument(
        '-u', default='unmapped_features.txt', metavar='FILE',
        help='write unmapped features to FILE; default is "unmapped_features.txt"',
    )
    outgrp.add_argument(
        '-exclude_partial', action='store_true',
        help='write partial mappings below -s and -a threshold to unmapped_features.txt; if true '
             'partial/low sequence identity mappings will be included in the gff file with '
             'partial_mapping=True, low_identity=True in comments'
    )
    outgrp.add_argument(
        '-dir', default='intermediate_files', metavar='DIR',
        help='name of directory to save intermediate fasta and SAM files; default is "intermediate_files"',
    )

    aligngrp = parser.add_argument_group('Alignments')
    aligngrp.add_argument('-mm2_options', metavar='=STR', type=str, default='-a --end-bonus '
                                                                            '5 --eqx -N 50 '
                                                                            '-p 0.5',
                          help='space delimited minimap2 parameters. By default ="-a --end-bonus 5 --eqx -N 50 -p 0.5"')
    aligngrp.add_argument(
        '-a', default=0.5, metavar='A', type=float,
        help='designate a feature mapped only if it aligns with coverage ≥A; by default A=0.5',
    )
    aligngrp.add_argument(
        '-s', default=0.5, metavar='S', type=float,
        help='designate a feature mapped only if its child features (usually exons/CDS) align '
             'with sequence identity ≥S; by default S=0.5'
    )
    aligngrp.add_argument(
        '-d', metavar='D', default=2.0, type=float,
        help='distance scaling factor; alignment nodes separated by more than a factor of D in '
             'the target genome will not be connected in the graph; by default D=2.0'
    )
    aligngrp.add_argument(
        '-flank', default=0, metavar='F', type=float, help="amount of flanking sequence to align as a "
                                                           "fraction [0.0-1.0] of gene length. This can improve gene "
                                                           "alignment where gene structure  differs between "
                                                           "target and "
                                                           "reference; by default F=0.0")

    pcgrp = parser.add_argument_group('Protein-Coding Genes')
    pcgrp.add_argument(
        '--no_prot_prior', default=False, action='store_true', 
        help='disable heuristics that prioritizes protein-coding genes during lift-over. [default=False]'
    )
    pcgrp.add_argument(
        '--prot_S', default=0.95, metavar='S', type=float,
        help='protein S-Score. When protein prioritization is enabled, the -s score for lifting over proteins. [default=0.95]'
    )

    regiongrp = parser.add_argument_group('Special genomic regions')
    regiongrp.add_argument(
        '--chrY_separate', default=False, action='store_true', 
        help='separate annotation for chromosome Y. Requires --annot_2. Applies ' 
             'protein prioritization heuristics separately to chrY when '
             'protein-prioritize is enabled. [default=False]'
    )
    regiongrp.add_argument(
        '--rDNA_separate', default=False, action='store_true',
        help='separate annotation for rDNA arrays. Requires --annot_2. [default=False]',
    )
    regiongrp.add_argument(
        '--annot_2', metavar='FILE',
        help='secondary annotation file. Required when either (or both) of previous two options is enabled.'
    )

    parser.add_argument('-V', '--version', help='show program version', action='version', version='v1.6.3')
    parser.add_argument(
        '-p', default=1, type=int, metavar='P', help='use p parallel processes to accelerate alignment; by default p=1'
    )
    parser.add_argument('-m', help='Minimap2 path', metavar='PATH')
    parser.add_argument('-f', metavar='TYPES', help='list of feature types to lift over')
    parser.add_argument(
        '-infer_genes', action='store_true',
        help='use if annotation file only includes transcripts, exon/CDS features'
    )
    parser.add_argument(
        '-infer_transcripts', action='store_true', required=False,
        help='use if annotation file only includes exon/CDS features and does not include transcripts/mRNA'
    )
    parser.add_argument(
        '-chroms', metavar='TXT', help='comma seperated file with corresponding chromosomes in '
                                       'the reference,target sequences',
    )
    parser.add_argument(
        '-unplaced', metavar='TXT',
        help='text file with name(s) of unplaced sequences to map genes from after genes from '
             'chromosomes in chroms.txt are mapped; default is "unplaced_seq_names.txt"',
    )
    parser.add_argument('-copies', action='store_true', help='look for extra gene copies in the target genome')
    parser.add_argument(
        '-sc', default=1.0, metavar='SC', type=float,
        help='with -copies, minimum sequence identity in exons/CDS for which a gene is considered '
             'a copy; must be greater than -s; default is 1.0',
    )
    parser.add_argument('-overlap', default=0.1, metavar='O', help="maximum fraction [0.0-1.0] of overlap allowed by 2 "
                                                                   "features; by default O=0.1", type=float)
    parser.add_argument('-mismatch', default=2, metavar='M', help="mismatch penalty in exons when finding best "
                                                                  "mapping; by default M=2", type=int)
    parser.add_argument('-gap_open', default=2, metavar='GO', help="gap open penalty in exons when finding best "
                                                                   "mapping; by default GO=2", type=int)
    parser.add_argument('-gap_extend', default=1, metavar='GE', help="gap extend penalty in exons when finding best "
                                                                     "mapping; by default GE=1", type=int)
    parser.add_argument('-subcommand', required=False,  help=argparse.SUPPRESS)
    parser.add_argument('-polish', required=False, action='store_true', default = False,
                        help="re-align exons to restore proper coding sequences in cases of start/stop codon loss or in-frame stop codons")
    parser.add_argument('-cds', required=False, action="store_true", default=True, help="annotate status of each CDS "
                                                                                        "(partial, missing start, "
                                                                                        "missing stop, inframe stop "
                                                                                        "codon)")
    parser._positionals.title = 'Required input (sequences)'
    parser._optionals.title = 'Miscellaneous settings'
    parser._action_groups = [parser._positionals, refrgrp, outgrp, aligngrp, pcgrp, regiongrp, parser._optionals]
    args = parser.parse_args(arglist)
    if '-a' not in args.mm2_options:
        args.mm2_options += ' -a'
    if '--eqx' not in args.mm2_options:
        args.mm2_options += ' --eqx'
    if '-N' not in args.mm2_options:
        args.mm2_options += " -N 50"
    if '-p' not in args.mm2_options:
        args.mm2_options += " -p 0.5"
    if '--end-bonus' not in args.mm2_options:
        args.mm2_options += "--end-bonus 5"
    if float(args.s) > float(args.sc):
        parser.error("-sc must be greater than or equal to -s")
    if args.chroms is None and args.unplaced is not None:
        parser.error("-unplaced must be used with -chroms")
    if args.chrY_separate and not args.annot_2:
        parser.error("--chrY-separate must be used with --annot-2")
    if args.rDNA_separate and not args.annot_2:
        parser.error("--rDNA_separate must be used with --annot-2")
    return args