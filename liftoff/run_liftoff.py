from liftoff import write_new_gff, liftover_types, extract_features, polish, align_features, lift_features
import argparse
from pyfaidx import Fasta, Faidx
import subprocess
import os






def main(arglist=None):
    args = parse_args(arglist)
    run_all_liftoff_steps(args)


def run_all_liftoff_steps(args):
    lifted_feature_list = {}
    unmapped_features = []

    # ---------------------------------------------
    # build database
    # ---------------------------------------------
    if args.chroms:
        ref_chroms, target_chroms = parse_chrm_files(args.chroms)
    else:
        ref_chroms, target_chroms = [args.reference], [args.target]
    parent_features_to_lift = get_parent_features_to_lift(args.f)

    all_features = parent_features_to_lift + ['gene_pc', 'gene_pseudo']
    feature_hierarchy, feature_db, ref_parent_order = extract_features.extract_features_to_lift(ref_chroms, 'chrm_by_chrm', all_features, args)

    if not args.prot_prioritize:
        # ---------------------------------------------
        # regular mode
        # ---------------------------------------------
        log('Aligning all feature types.')
        liftover_types.lift_original_annotation(ref_chroms, target_chroms, lifted_feature_list, args, unmapped_features, feature_db, feature_hierarchy, ref_parent_order, all_features)
        if len(unmapped_features) > 0 and target_chroms[0] != args.target:
            log("Mapping unaligned features against all.")
            unmapped_features = liftover_types.map_unmapped_genes_agaisnt_all(unmapped_features, [args.reference], [args.target], lifted_feature_list, feature_db, feature_hierarchy, ref_parent_order, args, all_features)
        if args.copies:
            log("Mapping extra copies against all.")
            liftover_types.map_extra_copies([args.reference], [args.target], lifted_feature_list, feature_hierarchy, feature_db, ref_parent_order, args, all_features)
    else:
        # ---------------------------------------------
        # liftover protein-coding genes
        # ---------------------------------------------
        log('Aligning protein-coding features.')
        liftover_types.lift_original_annotation(ref_chroms, target_chroms, lifted_feature_list, args, unmapped_features, feature_db, feature_hierarchy, ref_parent_order, ['gene_pc'])
        if len(unmapped_features) > 0 and target_chroms[0] != args.target:
            log("Mapping unaligned features against all.")
            unmapped_features = liftover_types.map_unmapped_genes_agaisnt_all(unmapped_features, [args.reference], [args.target], lifted_feature_list, feature_db, feature_hierarchy, ref_parent_order, args, ['gene_pc'])
        if args.copies:
            log("Mapping extra copies against all.")
            liftover_types.map_extra_copies([args.reference], [args.target], lifted_feature_list, feature_hierarchy, feature_db, ref_parent_order, args, ['gene_pc'])
        # ---------------------------------------------
        # liftover pseudogenes
        # ---------------------------------------------
        log('Aligning pseudogene features.')
        liftover_types.lift_original_annotation(ref_chroms, target_chroms, lifted_feature_list, args, unmapped_features, feature_db, feature_hierarchy, ref_parent_order, ['gene_pseudo'])
        if len(unmapped_features) > 0 and target_chroms[0] != args.target:
            log("Mapping unaligned features against all.")
            unmapped_features = liftover_types.map_unmapped_genes_agaisnt_all(unmapped_features, [args.reference], [args.target], lifted_feature_list, feature_db, feature_hierarchy, ref_parent_order, args, ['gene_pseudo'])
        if args.copies:
            log("Mapping extra copies against all.")
            liftover_types.map_extra_copies([args.reference], [args.target], lifted_feature_list, feature_hierarchy, feature_db, ref_parent_order, args, ['gene_pseudo'])
        # ---------------------------------------------
        # align other types and unplaced genes
        # ---------------------------------------------
        log('Aligning other feature types.')
        liftover_types.lift_original_annotation(ref_chroms, target_chroms, lifted_feature_list, args, unmapped_features, feature_db, feature_hierarchy, ref_parent_order, parent_features_to_lift)
        if len(unmapped_features) > 0 and target_chroms[0] != args.target:
            log("Mapping unaligned features against all.")
            unmapped_features = liftover_types.map_unmapped_genes_agaisnt_all(unmapped_features, [args.reference], [args.target], lifted_feature_list, feature_db, feature_hierarchy, ref_parent_order, args, parent_features_to_lift)
        if args.copies:
            log("Mapping extra copies against all.")
            liftover_types.map_extra_copies([args.reference], [args.target], lifted_feature_list, feature_hierarchy, feature_db, ref_parent_order, args, parent_features_to_lift)

    if args.unplaced and args.chroms:
        log("Mapping unplaced genes")
        ref_chroms, target_chroms = parse_chrm_files(args.unplaced)[0], [args.target]
        liftover_types.map_unplaced_genes(unmapped_features, ref_chroms, target_chroms, lifted_feature_list, feature_db, feature_hierarchy, ref_parent_order, args, all_features)

    # ---------------------------------------------
    # output and post-processing
    # ---------------------------------------------
    if args.cds or args.polish:
        check_cds(lifted_feature_list, feature_hierarchy, args)
    write_new_gff.write_new_gff(lifted_feature_list, args, feature_db)

    if args.polish:
        log("Polishing annotations")
        find_and_polish_broken_cds(args, lifted_feature_list, feature_hierarchy, ref_chroms, target_chroms, unmapped_features, feature_db, ref_parent_order, ['gene_pc'])
        if args.o == 'stdout':
            write_new_gff.write_new_gff(lifted_feature_list, args, feature_db)
        else:
            old_o = args.o; args.o += '_polished'
            write_new_gff.write_new_gff(lifted_feature_list, args, feature_db)
            args.o = old_o

    if args.fix_orfs:
        run_fix_orfs(args)

    with open(args.u, 'w') as file:
        for f in unmapped_features:
            file.write(f.id + "\n")
    
    log('Liftoff complete.')


def log(msg):
    print(f'[Info]: {msg}')


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
        '--prot_prioritize', default=True, action='store_true', 
        help='heuristics that prioritizes protein-coding genes during lift-over. [default=True]'
    )
    pcgrp.add_argument(
        '--prot_S', default=0.9, type=float, metavar='FLOAT',
        help='protein S-Score. When --prot_prioritize is True, the -s score for lifting over proteins. [default=0.9]',
    )
    pcgrp.add_argument(
        '--fix_orfs', default=False, action='store_true',
        help='fix open-reading-frames of CDSes. Similar to --polish and could '
             'be used with --polish but with more stringent checks and fixes for '
             'in-frame stop codons and protein truncations. Requires --prot. [default=False]'
    )
    pcgrp.add_argument(
        '--prot', metavar='FILE',
        help='protein sequence file. Required when --fix_orfs is True. Recommand '
             'to use "hg38" for human genome.'
    )

    regiongrp = parser.add_argument_group('Special genomic regions')
    regiongrp.add_argument(
        '--chrY_separate', default=False, action='store_true', 
        help='separate annotation for chromosome Y. Requires --annot_2. Applies ' 
             'protein prioritization heuristics separately to chrY when '
             '--prot_prioritize is enabled. [default=False]'
    )
    regiongrp.add_argument(
        '--rDNA_separate', default=False, action='store_true',
        help='separate annotation for rDNA arrays. Requires --annot_2. [default=False]',
    )
    regiongrp.add_argument(
        '--annot_2', metavar='FILE',
        help='secondary annotation file. Required when either (or both) of previous two options is enabled.'
    )
    regiongrp.add_argument(
        '--short_gene', action='store_true',
        help='handles short regions like VDJ regions that are poorly managed by '
             'default minimap2 parameters. Uses minimap2 -sr for genes between '
             '30 and 100 base pairs, and Bowtie2 for genes shorter than 30 bps. '
             '[default=False]'
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
        parser.error("-unplaced must be used with -chroms")\

    if args.chrY_separate and not args.annot_2:
        parser.error("--chrY-separate must be used with --annot-2")
    if args.rDNA_separate and not args.annot_2:
        parser.error("--rDNA_separate must be used with --annot-2")
    if args.fix_orfs and not args.prot:
        parser.error("--fix_orfs must be used with args.prot")
    return args


def parse_chrm_files(chroms_file):
    with open(chroms_file, 'r') as chroms:
        ref_chroms, target_chroms = [], []
        for line in chroms.readlines():
            kv_pair = line.rstrip().split(",")
            ref_chroms.append(kv_pair[0])
            assert(len(kv_pair) == 2)
            target_chroms.append(kv_pair[1])
    return ref_chroms, target_chroms


def get_parent_features_to_lift(feature_types_file):
    feature_types = ["gene"]
    if feature_types_file:
        with open(feature_types_file, 'r') as file:
            for line in file.readlines():
                if line.rstrip() not in feature_types:
                    feature_types.append(line.rstrip())
    return feature_types


def find_and_polish_broken_cds(args, lifted_feature_list,feature_hierarchy, ref_chroms, target_chroms,
                               unmapped_features, feature_db, ref_parent_order, feature_types):
    args.subcommand = "polish"
    polish_lifted_features = {}
    ref_fa, target_fa = Fasta(args.reference), Fasta(args.target)
    for target_feature in lifted_feature_list:
        aligned_segments_new = {}
        if polish.polish_annotations(lifted_feature_list, ref_fa, target_fa, args, feature_hierarchy, target_feature):
            aligned_segments = align_features.align_features_to_target(ref_chroms, target_chroms, args, feature_hierarchy, "chrm_by_chrm", unmapped_features, feature_types)
            aligned_segments_new[target_feature] = list(aligned_segments.values())[0]
            for seg in aligned_segments_new[target_feature]:
                seg.query_name = target_feature
            args.d = 100000000
            s_var = args.prot_S if args.prot_prioritize else args.s
            lift_features.lift_all_features(aligned_segments_new, args.a, feature_db, feature_hierarchy,
                                            unmapped_features, polish_lifted_features, s_var, None, args,
                                            ref_parent_order, feature_types)

    check_cds(polish_lifted_features, feature_hierarchy, args)
    for feature in polish_lifted_features:
        original_feature = lifted_feature_list[feature][0]
        polished_feature = polish_lifted_features[feature][0]
        replace = False
        if 'valid_ORFs' not in polished_feature.attributes or int(polished_feature.attributes['valid_ORFs'][0]) > \
                int(original_feature.attributes['valid_ORFs'][0]):
            replace = True
        elif polished_feature.attributes['valid_ORFs'][0] == original_feature.attributes['valid_ORFs'][0]:
            if polished_feature.attributes['sequence_ID'][0] > original_feature.attributes['sequence_ID'][0]:
                replace = True
            elif polished_feature.attributes['coverage'][0] > original_feature.attributes['coverage'][0]:
                replace = True
        if replace:
            lifted_feature_list[feature] = polish_lifted_features[feature]


def check_cds(feature_list, feature_hierarchy, args):
    ref_faidx, target_faidx = Fasta(args.reference), Fasta(args.target)
    for target_feature in feature_list:
        target_sub_features = polish.get_sub_features(feature_list, target_feature)
        ref_sub_features = polish.get_sub_features(feature_hierarchy.children, target_sub_features[0].id)
        polish.find_and_check_cds(target_sub_features, ref_sub_features, ref_faidx,
                                                             target_faidx, feature_list[target_feature])


def run_fix_orfs(args):
    log("Run fix_orfs.sh")
    log("Check CDSes in reference.")
    cmd = ['/ccb/sw/packages/fix_orfs/fix_orfs.sh', '-g', args.reference, '-p', args.prot, '-a', args.g, '-t', args.p, '-o', os.path.join(args.dir, 'p1'), '-d']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, _ = proc.communicate()
    print(stdout.decode())

    log("Fix protein coding annotations.")
    broken_fn = os.path.join(args.dir, 'p1.broken.txt')
    if args.o != 'stdout' and args.polish:
        args.o += "_polished"
    cmd = ['/ccb/sw/packages/fix_orfs/fix_orfs.sh', '-g', args.target, '-p', args.prot, '-a', args.o, '-t', args.p, '-o', os.path.join(args.dir, 'final'), '-b', broken_fn]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, _ = proc.communicate()
    print(stdout.decode())
    log(f"Fixed annotation in {args.dir}/final.adjusted_cds.gff")


if __name__ == "__main__":
    main()
