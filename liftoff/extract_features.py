import gffutils
from pyfaidx import Fasta
from liftoff import liftoff_utils, feature_hierarchy, new_feature
from liftoff.liftoff_utils import LiftoverType, find_parent_order
import os
import sys
import numpy as np
import ujson as json


def create_hierarchical_database(args, parents_to_lift):
    if not os.path.exists(args.dir):
        os.mkdir(args.dir)
    gffutils.constants.ignore_url_escape_characters = True
    feature_db = build_database(args.db, args.g, not args.infer_transcripts, not args.infer_genes, not args.no_prot_prior)
    feature_hierarchy, parent_order = seperate_parents_and_children(feature_db, parents_to_lift)
    return feature_hierarchy, feature_db, parent_order


def transform (f):
    if f.featuretype == 'gene':
        if 'gene_type' in f.attributes and f['gene_type'][0] == 'protein_coding':
            f.featuretype = 'gene_pc'
        elif 'gene_biotype' in f.attributes and f['gene_biotype'][0] == 'protein_coding':
            f.featuretype = 'gene_pc'
    return f 


def build_database(db, gff_file, disable_transcripts, disable_genes, protein_priority):
    feature_db = None
    if db:
        if protein_priority and not db.endswith("_curated_db"):
            print('[Warning]: Improper suffix; may not be able to extract gene type information.')
        print(f"[Info]: Extracting features from {db}")
        feature_db = gffutils.FeatureDB(db)
        if protein_priority:
            pc     = feature_db.count_features_of_type(featuretype='gene_pc')
            other  = feature_db.count_features_of_type(featuretype='gene')
            if pc == 0:
                sys.exit(f'[Fatal]: Extracted {pc} protein-coding genes.')
            else:
                print(f'[Info]: Extracted {pc} protein-coding genes.\n'
                      f'[Info]: Extracted {other} other genes.')
        else:
            all_g  = feature_db.count_features_of_type(featuretype='gene')
            print(f'[Info]: Extracted {all_g} genes.')
    else:
        try:
            print(f"[Info]: Extracting features from {gff_file}")
            if protein_priority:
                suffix = "_curated_db"
                func = transform 
            else:
                suffix = "_db"
                func = None
            feature_db = gffutils.create_db(
                gff_file, gff_file + suffix, 
                merge_strategy="create_unique", force=True,
                disable_infer_transcripts=disable_transcripts,
                disable_infer_genes=disable_genes, 
                verbose=True, transform=func)
            print('[Info]: Database build succeeded.')
            if protein_priority:
                pc     = feature_db.count_features_of_type(featuretype='gene_pc')
                other  = feature_db.count_features_of_type(featuretype='gene')
                print(f'[Info]: Extracted {pc} protein-coding genes.\n'
                      f'[Info]: Extracted {other} other genes.')
            else:
                all_g  = feature_db.count_features_of_type(featuretype='gene')
                print(f'[Info]: Extracted {all_g} genes.')
        except Exception as e:
            print(f'[Error]: Exception occurred while creating database: {e}')
            find_problem_line(gff_file)
    return feature_db


def find_problem_line(gff_file):
    f = open(gff_file, 'r')
    lines = f.readlines()
    for i in range (len(lines)):
        line = lines[i]
        if line[0] != "#":
            try:
                gffutils.create_db(line, ":memory:", from_string=True, force=True)
            except:
                sys.exit("[Fatal]: Incorrect GFF/GTF syntax on line " + str(i + 1))


def seperate_parents_and_children(feature_db, parent_types_to_lift):
    c = feature_db.conn.cursor()
    relations = [list(feature) for feature in c.execute(
                ''' SELECT * FROM relations
                    JOIN features AS a ON a.id = relations.parent
                    JOIN features AS b ON b.id = relations.child
                ''') if feature[0] != feature[1]]
    all_ids   = [list(feature)[0] for feature in c.execute(
                '''SELECT * FROM features''')]
    all_parent_ids   = [relation[0] for relation in relations]
    all_children_ids = [relation[1] for relation in relations]
    lowest_children = np.setdiff1d(all_ids, all_parent_ids)
    highest_parents = np.setdiff1d(all_ids, all_children_ids)
    intermediates = set(all_children_ids).intersection(set(all_parent_ids))
    parent_dict, child_dict, intermediate_dict = dict(), dict(), dict()
    add_parents(parent_dict, child_dict, highest_parents, parent_types_to_lift, feature_db)
    add_children(parent_dict, child_dict, lowest_children, feature_db)
    add_intermediates(intermediates, intermediate_dict, feature_db)
    parent_order          = find_parent_order([parent for parent in list(parent_dict.values()) if parent])
    ref_feature_hierarchy = feature_hierarchy.feature_hierarchy(parent_dict, intermediate_dict, child_dict)
    return ref_feature_hierarchy, parent_order


def add_parents(parent_dict, child_dict, highest_parents, parent_types_to_lift, feature_db):
    c = feature_db.conn.cursor()
    cond = ', '.join(f'"{w}"' for w in highest_parents)
    for result in c.execute(f'''SELECT * FROM features WHERE id IN ({cond})'''):
        feature_tup = tuple(result)
        parent = new_feature.new_feature(feature_tup[0], feature_tup[3], feature_tup[1], feature_tup[2],feature_tup[7],
                                         feature_tup[4], feature_tup[5], feature_tup[8], json.loads(feature_tup[9]))
        if parent.featuretype in parent_types_to_lift or parent.featuretype in ['gene_pc', 'gene_pseudo']:
            parent_dict[parent.id] = parent
            child_dict[parent.id] = []


def add_children(parent_dict, child_dict, lowest_children, feature_db):
    c = feature_db.conn.cursor()
    cond = ', '.join(f'"{w}"' for w in lowest_children)
    c.execute(f'''SELECT * FROM relations JOIN features ON features.id  = relations.child WHERE relations.child IN ({cond})''')
    results = c.fetchall()
    added_children_ids = []
    for result in results:
        feature_tup = tuple(result)
        parent = feature_tup[0]
        if parent in parent_dict:
            child = new_feature.new_feature(feature_tup[3], feature_tup[6], feature_tup[4], feature_tup[5], feature_tup[10], 
                                            feature_tup[7], feature_tup[8], feature_tup[11], json.loads(feature_tup[12]))
            if child.featuretype != "intron":
                if "Parent" not in child.attributes:
                    add_parent_tag(child, feature_db)
                child_dict[parent].append(child)
                added_children_ids.append(child.id)
    single_level_features = np.setdiff1d(lowest_children, added_children_ids)
    for feature in single_level_features:
        if feature in parent_dict:
            child_dict[feature] = [parent_dict[feature]]


def add_parent_tag(feature, feature_db):
    parent_id = ""
    parents = [parent for parent in feature_db.parents(feature.id, level=1) if feature.id != parent.id]
    if len(parents) > 0:
        parent_id = parents[0].id
    else:
        parents = [parent for parent in feature_db.parents(feature.id) if feature.id != parent.id]
        if len(parents) > 0:
            parent_id = parents[0].id
    feature.attributes["Parent"] = [parent_id]



def add_intermediates(intermediate_ids, intermediate_dict, feature_db):
    c = feature_db.conn.cursor()
    cond = ', '.join(f'"{w}"' for w in intermediate_ids)
    for result in c.execute(f'''SELECT * FROM features WHERE id IN ({cond})'''):
        feature_tup = tuple(result)
        intermediate_feature = new_feature.new_feature(feature_tup[0], feature_tup[3], feature_tup[1], feature_tup[2], feature_tup[7],
                                                       feature_tup[4], feature_tup[5], feature_tup[8], json.loads(feature_tup[9]))
        intermediate_dict[intermediate_feature.id] = intermediate_feature
        if "Parent" not in intermediate_feature.attributes:
            add_parent_tag(intermediate_feature, feature_db)


def get_gene_sequences(parent_dict, ref_chroms, args, liftover_type):
    fai = Fasta(args.reference)
    if liftover_type == LiftoverType.UNPLACED:
        open(args.dir + "/unplaced_genes.fa", 'w')
    for chrom in ref_chroms:
        parent_feature = sorted(list(parent_dict.values()), key=lambda x: x.seqid)
        if len(parent_feature) == 0:
            sys.exit("[Fatal]: GFF does not contain any gene features. Use -f to provide a list of other feature types to lift over.")
        fasta_hdl = get_fasta_hdl(chrom, args.reference, liftover_type, args.dir)
        write_gene_sequences_to_file(chrom, args.reference, fai, parent_feature, fasta_hdl, args)
        fasta_hdl.close()



def get_fasta_hdl(chrom_name, reference_fasta_name, ltype, inter_files):
    if chrom_name == reference_fasta_name and ltype in [LiftoverType.ONE2ONE, LiftoverType.COPIES]:
        fasta_out_name = "reference_all"
    elif ltype == LiftoverType.UNMAPPED:
        fasta_out_name = "unmapped_to_expected_chrom"
    elif ltype == LiftoverType.UNPLACED:
        fasta_out_name = "unplaced"
    else:
        fasta_out_name = chrom_name
    open_mode = 'a' if ltype == LiftoverType.UNPLACED else 'w'
    return open(f'{inter_files}/{fasta_out_name}_gene.fa', open_mode)


def write_gene_sequences_to_file(chrom_name, reference_fasta_name, reference_fasta_idx, parents, fasta_out, args):
    if chrom_name == reference_fasta_name:
        current_chrom = parents[0].seqid
    else:
        current_chrom = chrom_name
    chrom_seq = reference_fasta_idx[current_chrom][:].seq
    for parent in parents:
        if parent.seqid != current_chrom and chrom_name == reference_fasta_name:
            current_chrom = parent.seqid
            chrom_seq = reference_fasta_idx[current_chrom][:].seq
        if parent.seqid == chrom_name or chrom_name == reference_fasta_name:
            gene_length = parent.end - parent.start + 1
            parent.start = round(max(1, parent.start - args.flank * gene_length))
            parent.end = round(min(parent.end + args.flank * gene_length, len(chrom_seq)))
            parent_seq = chrom_seq[parent.start - 1: parent.end]
            fasta_out.write(">" + parent.id + "\n" + str(parent_seq) + "\n")
