

#ifndef MIRNA_VALIDATION_H
#define MIRNA_VALIDATION_H

#include "defs.h"
#include "candidates.h"
#include "util.h"

/* forward declaration of struct in coverage.h */
struct extended_candidate;
struct candidate_subsequence;
struct chrom_coverage;

int validate_mature_micro_rnas(struct extended_candidate *ecand,
                               struct configuration_params *config);
int mirna_coverage_compare(const void *a, const void *b);
int mirna_total_coverage_compare(const void *a, const void *b);
int filter_mature_micro_rnas(struct extended_candidate *ecand,
                             struct chrom_coverage *chrom_cov,
                             struct configuration_params *config);
int determine_mature_arm(struct extended_candidate *ecand);
int check_for_valid_folding(struct extended_candidate *ecand);
int check_for_terminal_mismatches(struct extended_candidate *ecand);
int check_for_adjacent_mismatches(struct extended_candidate *ecand);
int find_star_micro_rna(struct candidate_subsequence **result,
                        struct micro_rna_candidate *cand,
                        struct candidate_subsequence *mature_mirna,
                        struct chrom_coverage *chrom_cov,
                        struct configuration_params *config);
int find_matching_bracket_index(u32 *result, u32 target, char *structure,
                                size_t n);
#endif