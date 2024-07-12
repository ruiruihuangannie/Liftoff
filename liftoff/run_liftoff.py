from liftoff import parser, write_new_gff, liftover_types, extract_features, polish, align_features, lift_features
from liftoff.liftoff_utils import LiftoverType, parse_chrm_files
import subprocess
import os


def main(arglist=None):
    args = parser.parse_args(arglist)

    # redirect minimap2 outupts to log file
    log_fn = os.path.join(args.dir, 'minimap.log')
    if os.path.isfile(log_fn): os.remove(log_fn)
    
    # build database
    if args.chroms:
        ref_chroms, tgt_chroms = parse_chrm_files(args.chroms)
    else:
        ref_chroms, tgt_chroms = [args.reference], [args.target]
    
    parent_features_to_lift = get_parent_features_to_lift(args.f)
    all_features = parent_features_to_lift if args.no_prot_prior else parent_features_to_lift + ['gene_pc']

    feature_hierarchy, feature_db, ref_parent_order = extract_features.create_hierarchical_database(
        args,
        all_features
    )
    
    Liftoff(
        ref_chroms=ref_chroms,
        tgt_chroms=tgt_chroms,
        features=all_features,
        feature_hierarchy=feature_hierarchy,
        feature_db=feature_db,
        ref_parent_order=ref_parent_order,
        args=args
    ).run_liftoff()


class Liftoff:
    def __init__(
        self,
        *,
        ref_chroms,
        tgt_chroms,
        features,
        feature_hierarchy,
        feature_db,
        ref_parent_order,
        args
    ):
        self.ref_chroms = ref_chroms
        self.tgt_chroms = tgt_chroms
        self.features  = features
        self.feature_hierarchy = feature_hierarchy
        self.feature_db = feature_db
        self.ref_parent_order = ref_parent_order
        self.lifted_features = dict()
        self.unmapped_features = list()
        self.settled_features = set()
        self.args = args
    
    def run_liftoff(self):
        self.lift_single_match_annotation()
        self.lift_extra_copies_annotation()
        self.lift_unplaced_seq_annotation()
        self.write_unmapped_features2file()
        self.post_process_and_output2file()
        log('Liftoff complete.')
    
    def lift_single_match_annotation(self):
        if self.args.chroms:
            log('Aligning features chrom-by-chrom.')
            liftover_types.liftover(self, self.unmapped_features, LiftoverType.ONE2ONE)
            if len(self.unmapped_features) > 0:
                log("Mapping unaligned features against all.")
                unmapped_features = list()
                liftover_types.liftover(self, unmapped_features, LiftoverType.UNMAPPED)
                self.unmapped_features = unmapped_features
        else:
            log('Aligning features against all.')
            liftover_types.liftover(self, self.unmapped_features, LiftoverType.ONE2ONE)
    
    def lift_extra_copies_annotation(self):
        if self.args.copies:
            log("Mapping extra copies against all.")
            unmapped_features = list()
            liftover_types.liftover(self, unmapped_features, LiftoverType.COPIES)
            post_len = len(self.unmapped_features)


    def lift_unplaced_seq_annotation(self):
        if self.args.unplaced and self.args.chroms:
            log("Mapping unplaced genes")
            # ref_chroms, target_chroms = parse_chrm_files(args.unplaced)[0], [args.target]
            liftover_types.liftover(self, self.unmapped_features, LiftoverType.UNPLACED)
    
    def post_process_and_output2file(self):
        # output and post-processing
        if not self.args.cds and not self.args.polish:
            write_new_gff.write_new_gff(self)
        else:
            polish.check_cds(self.lifted_features, self.feature_hierarchy, self.args)
            write_new_gff.write_new_gff(self)
        if self.args.polish:
            log("Polishing annotations")
            polish.find_and_polish_broken_cds(self, self.args, LiftoverType.ONE2ONE)
            if self.args.o == 'stdout':
                write_new_gff.write_new_gff(self)
            else:
                self.args.o += '_polished'
                write_new_gff.write_new_gff(self)
    
    def write_unmapped_features2file(self):
        with open(self.args.u, 'w') as file:
            for f in self.unmapped_features:
                file.write(f.id + "\n")


def log(msg, debug_mode=True):
    if debug_mode:
        print(f'[Info]: {msg}')


def get_parent_features_to_lift(feature_types_file):
    feature_types = ["gene"]
    if feature_types_file:
        with open(feature_types_file, 'r') as file:
            for line in file.readlines():
                if line.rstrip() not in feature_types:
                    feature_types.append(line.rstrip())
    return feature_types


if __name__ == "__main__":
    main()