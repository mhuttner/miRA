#include "cluster.h"
#include "fasta.h"
#ifndef VFOLD_H
#define VFOLD_H

struct foldable_sequence {
  struct cluster *c;
  char *seq;
  size_t n;
  char *structure;
  size_t st_n;
  int st_start;
};

struct sequence_list {
  struct foldable_sequence **sequences;
  size_t n;
};

int vfold(int argc, char **argv);
int map_clusters(struct sequence_list **seq_list, struct cluster_list *c_list,
                 struct genome_sequence *seq_table);
int fold_sequences(struct sequence_list *seq_list);
int reverse_complement(struct foldable_sequence *s);
int free_sequence_list(struct sequence_list *seq_list);
int free_foldable_sequence(struct foldable_sequence *s);

#endif