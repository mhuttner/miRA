#include <stdio.h>
#include "cluster.h"
#include "Lfold/Lfold.h"
#include "fasta.h"
#include "util.h"
#ifndef VFOLD_H
#define VFOLD_H

struct structure_info {
  char *structure_string;
  size_t n;
  int start;
  double mfe;
  double pvalue;
  double mean;
  double sd;
  int external_loop_count;
  double paired_fraction;
  int stem_start;
  int stem_end;
  int stem_start_with_mismatch;
  int stem_end_with_mismatch;
  int is_valid;
};

struct foldable_sequence {
  struct cluster *c;
  char *seq;
  size_t n;
  struct structure_info *structure;
};

struct sequence_list {
  struct foldable_sequence **sequences;
  size_t n;
};

int vfold(int argc, char **argv);
int map_clusters(struct sequence_list **seq_list, struct cluster_list *c_list,
                 struct genome_sequence *seq_table);
int fold_sequences(struct sequence_list *seq_list,
                   struct configuration_params *config);
int write_json_result(struct sequence_list *seq_list, char *filename);
int calculate_mfe_distribution(struct foldable_sequence *fs,
                               int permutation_count);
int find_optimal_structure(struct structure_list *s_list,
                           struct foldable_sequence *fs,
                           struct configuration_params *config);
int check_folding_constraints(struct foldable_sequence *fs,
                              struct configuration_params *config);
int check_pvalue(struct foldable_sequence *fs,
                 struct configuration_params *config);
int write_foldable_sequence(FILE *fp, struct foldable_sequence *fs);
int reverse_complement(struct foldable_sequence *s);
int get_rand_int(int max);
int fisher_yates_shuffle(char *seq, int n);
int free_structure_info(struct structure_info *info);
int free_sequence_list(struct sequence_list *seq_list);
int free_foldable_sequence(struct foldable_sequence *s);

#endif