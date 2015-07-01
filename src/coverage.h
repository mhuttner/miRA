

#ifndef COVERAGE_H
#define COVERAGE_H

#include "defs.h"
#include "uthash.h"
#include "parse_sam.h"
#include "candidates.h"
#include "util.h"
#include "reads.h"

struct chrom_coverage {
  char name[1024];
  long length;
  u32 *coverage_plus;
  u32 *coverage_minus;
  UT_hash_handle hh;
};
enum strand_arm { FIVE_PRIME = 5, THREE_PRIME = 3 };

struct candidate_subsequence_list {
  struct candidate_subsequence **mature_sequences;
  size_t capacity;
  size_t n;
};

struct candidate_subsequence {
  u32 start;
  u32 end;
  u64 coverage;
  double paired_fraction;
  i8 is_valid;
  i8 is_star;
  i8 is_artificial;
  enum strand_arm arm;
  struct unique_read_list *reads;
  struct candidate_subsequence *matching_sequence;
};

struct extended_candidate {
  struct micro_rna_candidate *cand;
  struct candidate_subsequence_list *possible_micro_rnas;
  struct candidate_subsequence *mature_micro_rna;
  struct candidate_subsequence *star_micro_rna;
  u32 total_reads;
  double total_read_percent;
  i8 is_valid;
};

struct extended_candidate_list {
  struct extended_candidate **candidates;
  size_t n;
};

int coverage(int argc, char **argv);
int coverage_main(struct configuration_params *config, char *executable_file,
                  char *mira_file, char *sam_file, char *output_path);
int create_coverage_table(struct chrom_coverage **table, struct sam_file *sam);
int coverage_test_candidates(struct extended_candidate_list *ecand_list,
                             struct chrom_coverage **coverage_table,
                             struct sam_file *sam,
                             struct configuration_params *config);
int find_mature_micro_rnas(struct extended_candidate *ecand,
                           struct chrom_coverage *chrom_cov,
                           struct configuration_params *config);
int get_coverage_in_range(u64 *result, u32 *cov_list, u64 start, u64 end);
int extend_all_candidates(struct extended_candidate_list **ecand_list,
                          struct candidate_list *cand_list);

int create_extended_candidate(struct extended_candidate **ecand,
                              struct micro_rna_candidate *cand);
int
create_candidate_subseqence_list(struct candidate_subsequence_list **css_list);
int
append_candidate_subseqence_list(struct candidate_subsequence_list *css_list,
                                 struct candidate_subsequence *mature);
int
free_candidate_subsequence_list(struct candidate_subsequence_list *css_list);
int create_candidate_subseqence(struct candidate_subsequence **cs, u32 start,
                                u32 end, u64 coverage, double paired_fraction);
int free_extended_candidate_list(struct extended_candidate_list *ecand_list);
int free_extended_candidate(struct extended_candidate *ecand);

int free_candidate_subsequence(struct candidate_subsequence *cs);
int free_coverage_table(struct chrom_coverage **table);
#endif