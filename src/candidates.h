#include <stddef.h>
#include "defs.h"
#include "vfold.h"

#ifndef COVERAGE_H
#define COVERAGE_H

struct micro_rna_candidate {
  u64 id;
  char *chrom;
  char strand;
  u64 start;
  u64 end;
  char *sequence;
  char *structure;
  double mfe;
  double pvalue;
  double mean;
  double sd;
};

struct candidate_list {
  struct micro_rna_candidate **candidates;
  size_t n;
  size_t capacity;
};
int create_empty_candidate_list(struct candidate_list **cand_list);
int add_candidate_to_list(struct candidate_list *cand_list,
                          struct micro_rna_candidate *cand);
int convert_seq_list_to_cand_list(struct candidate_list **cand_list,
                                  struct sequence_list *seq_list);
int create_micro_rna_candidate(struct micro_rna_candidate **cand,
                               struct foldable_sequence *fs);
int write_candidate_file(struct candidate_list *cand_list, char *filename);
int write_candidate_line(FILE *fp, struct micro_rna_candidate *cand);
int read_candidate_file(struct candidate_list **cand_list, char *filename);
int parse_candidate_line(struct micro_rna_candidate **cand, char *line);
int free_candidate_list(struct candidate_list *cand_list);
int free_micro_rna_candidate(struct micro_rna_candidate *cand);

#endif