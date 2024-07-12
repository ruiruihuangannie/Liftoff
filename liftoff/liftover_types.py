from liftoff import fix_overlapping_features, lift_features, align_features, extract_features
from liftoff.liftoff_utils import LiftoverType, parse_chrm_files, convert_id_to_original

def liftover(obj, unmapped_features, ltype):
    ref_tgt = {
        LiftoverType.ONE2ONE:   tuple([obj.ref_chroms, obj.tgt_chroms]),
        LiftoverType.UNMAPPED:  tuple([[obj.args.reference], [obj.args.target]]),
        LiftoverType.COPIES:    tuple([[obj.args.reference], [obj.args.target]]),
        LiftoverType.UNPLACED:  None if not obj.args.unplaced else tuple([parse_chrm_files(obj.args.unplaced)[0], [obj.args.target]])
    }

    if ltype in [LiftoverType.ONE2ONE, LiftoverType.COPIES]:
        extract_features.get_gene_sequences(
            obj.feature_hierarchy.parents, 
            ref_tgt[ltype][0], 
            obj.args, 
            ltype,
        )
    elif ltype == LiftoverType.UNMAPPED:
        extract_features.get_gene_sequences(
            get_unmapped_genes(obj.unmapped_features), 
            ref_tgt[ltype][0], 
            obj.args, 
            ltype
        )
    elif ltype == LiftoverType.UNPLACED:
        extract_features.get_gene_sequences(
            get_unplaced_seq(obj.args.reference, obj.feature_hierarchy),
            ref_tgt[ltype][0], 
            obj.args, 
            ltype
        )
        
    aligned_segments = align_features.align_features_to_target(
        obj.args,
        obj.feature_hierarchy,
        unmapped_features,
        ref_tgt[ltype][0],
        ref_tgt[ltype][1],
        ltype
    )
    
    feature_locations = None
    lift_features.lift_all_features(
        obj,
        aligned_segments, 
        feature_locations,
        unmapped_features,
        ltype
    )
    fix_overlapping_features.fix_incorrectly_overlapping_features(
        obj,
        obj.lifted_features,
        obj.lifted_features,
        aligned_segments,
        unmapped_features,
        ltype
    )
    for entry in obj.lifted_features:
        obj.settled_features.add(convert_id_to_original(entry))


def get_unmapped_genes(unmapped_features):
    unmapped_dict = {}
    for feature in unmapped_features:
        unmapped_dict[feature.id] = feature
    return unmapped_dict


def get_unplaced_seq(ref_chroms, feature_hierarchy):
    unplaced_dict = {}
    for feature_name in feature_hierarchy.parents:
        feature = feature_hierarchy.parents[feature_name]
        if feature.seqid in ref_chroms:
            unplaced_dict[feature.id] = feature
    return unplaced_dict

