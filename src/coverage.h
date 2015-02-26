#include "defs.h"
#include "uthash.h"
#include "parse_sam.h"
#include "candidates.h"
#include "util.h"

#ifndef COVERAGE_H
#define COVERAGE_H

struct chrom_coverage {
  char name[1024];
  long length;
  u32 *coverage_plus;
  u32 *coverage_minus;
  UT_hash_handle hh;
};
enum strand_arm { FIVE_PRIME = 5, THREE_PRIME = 3 };

struct candidate_subsequence {
  u32 start;
  u32 end;
  u64 coverage;
  double paired_fraction;
  int is_valid;
  enum strand_arm arm;
};

struct extended_candidate {
  struct micro_rna_candidate *cand;
  struct candidate_subsequence *mature_micro_rna;
  struct candidate_subsequence *star_micro_rna;
  int is_valid;
};

struct extended_candidate_list {
  struct extended_candidate **candidates;
  size_t n;
};

int coverage(int argc, char **argv);
int create_coverage_table(struct chrom_coverage **table, struct sam_file *sam);
int coverage_test_candidates(struct extended_candidate_list *ecand_list,
                             struct chrom_coverage **coverage_table,
                             struct configuration_params *config);
int find_mature_micro_rna(struct extended_candidate *ecand,
                          struct chrom_coverage *chrom_cov,
                          struct configuration_params *config);
int find_star_micro_rna(struct extended_candidate *ecand,
                        struct chrom_coverage *chrom_cov,
                        struct configuration_params *config);
int find_matching_bracket_index(u32 *result, u32 target, char *structure,
                                size_t n);
int get_coverage_in_range(u64 *result, u32 *cov_list, u64 start, u64 end);
int extend_all_candidates(struct extended_candidate_list **ecand_list,
                          struct candidate_list *cand_list);
int validate_mature_micro_rna(struct extended_candidate *ecand,
                              struct configuration_params *config);
int determine_mature_arm(struct extended_candidate *ecand);
int check_for_valid_folding(struct extended_candidate *ecand);
int check_for_terminal_mismatches(struct extended_candidate *ecand);
int check_for_adjacent_mismatches(struct extended_candidate *ecand);
int create_extended_candidate(struct extended_candidate **ecand,
                              struct micro_rna_candidate *cand);
int create_candidate_subseqence(struct candidate_subsequence **cs, u32 start,
                                u32 end, u64 coverage, double paired_fraction);
int free_extended_candidate_list(struct extended_candidate_list *ecand_list);
int free_extended_candidate(struct extended_candidate *ecand);
int free_candidate_subsequence(struct candidate_subsequence *cs);
int free_coverage_table(struct chrom_coverage **table);
#endif