#include "mirna_validation.h"
#include "errors.h"
#include "util.h"
#include "coverage.h"
#include "structure_evaluation.h"
#include <stdlib.h>

int validate_mature_micro_rnas(struct extended_candidate *ecand,
                               struct configuration_params *config) {
  determine_mature_arm(ecand);
  check_for_valid_folding(ecand);
  if (config->allow_two_terminal_mismatches == 0) {
    check_for_terminal_mismatches(ecand);
  }
  if (config->allow_three_mismatches == 0) {
    check_for_adjacent_mismatches(ecand);
  }

  return E_SUCCESS;
}
int mirna_coverage_compare(const void *a, const void *b) {
  struct candidate_subsequence *csa = *(struct candidate_subsequence **)a;
  struct candidate_subsequence *csb = *(struct candidate_subsequence **)b;
  return (i64)csb->coverage - (i64)csa->coverage;
}
int mirna_total_coverage_compare(const void *a, const void *b) {
  struct candidate_subsequence *csa = *(struct candidate_subsequence **)a;
  struct candidate_subsequence *csb = *(struct candidate_subsequence **)b;
  return ((i64)csb->coverage + (i64)csb->matching_sequence->coverage) -
         ((i64)csa->coverage + (i64)csa->matching_sequence->coverage);
}

int filter_mature_micro_rnas(struct extended_candidate *ecand,
                             struct chrom_coverage *chrom_cov,
                             struct configuration_params *config) {
  struct candidate_subsequence_list *css_list = ecand->possible_micro_rnas;
  struct candidate_subsequence *star_css = NULL;
  struct candidate_subsequence *css = NULL;
  struct candidate_subsequence *css_tmp = NULL;
  size_t j = 0;
  for (size_t i = 0; i < css_list->n; i++) {
    css = css_list->mature_sequences[i];
    if (css->is_valid) {
      css_list->mature_sequences[j] = css_list->mature_sequences[i];
      j++;
    } else {
      free_candidate_subsequence(css_list->mature_sequences[i]);
      css_list->mature_sequences[i] = NULL;
    }
  }
  css_list->n = j;

  qsort(css_list->mature_sequences, css_list->n,
        sizeof(struct candidate_subsequence *), mirna_coverage_compare);

  for (size_t i = 0; i < css_list->n; i++) {
    css = css_list->mature_sequences[i];
    css->matching_sequence = NULL;
    if (css->is_star) {
      continue;
    }
    find_star_micro_rna(&star_css, ecand->cand, css, chrom_cov, config);
    if (star_css == NULL) {
      css->is_star = 1;
      continue;
    }
    for (size_t j = i + 1; j < css_list->n; j++) {
      css_tmp = css_list->mature_sequences[j];
      if (css_tmp->is_star) {
        continue;
      }
      int min_offset = config->min_dicer_offset;
      int max_offset = config->max_dicer_offset;
      i32 start_offset = abs((i32)css_tmp->start - (i32)star_css->start);
      i32 end_offset = abs((i32)css_tmp->end - (i32)star_css->end);
      if (min_offset <= start_offset && start_offset < max_offset &&
          min_offset <= end_offset && end_offset < max_offset) {
        css_tmp->is_star = 1;
        css->matching_sequence = css_tmp;
        css->matching_sequence->matching_sequence = css;
        break;
      }
    }
    if (css->matching_sequence == NULL) {
      css->matching_sequence = star_css;
      css->matching_sequence->matching_sequence = css;
    } else {
      free_candidate_subsequence(star_css);
    }
  }
  j = 0;
  for (size_t i = 0; i < css_list->n; i++) {
    css = css_list->mature_sequences[i];
    if (!(css->is_star)) {
      css_list->mature_sequences[j] = css_list->mature_sequences[i];
      j++;
    } else {
      css_list->mature_sequences[i] = NULL;
    }
  }
  css_list->n = j;
  qsort(css_list->mature_sequences, css_list->n,
        sizeof(struct candidate_subsequence *), mirna_total_coverage_compare);
  return E_SUCCESS;
}

int determine_mature_arm(struct extended_candidate *ecand) {
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence_list *css_list = ecand->possible_micro_rnas;
  struct candidate_subsequence *mature_mirna = NULL;
  for (size_t i = 0; i < css_list->n; i++) {
    mature_mirna = css_list->mature_sequences[i];

    u32 cand_center = (cand->end - cand->start + 1) / 2;
    u32 mature_center = (mature_mirna->end + mature_mirna->start + 1) / 2;
    if (mature_center > cand_center) {
      if (cand->strand == '+') {
        mature_mirna->arm = THREE_PRIME;
      } else {
        mature_mirna->arm = FIVE_PRIME;
      }
    } else {
      if (cand->strand == '+') {
        mature_mirna->arm = FIVE_PRIME;
      } else {
        mature_mirna->arm = THREE_PRIME;
      }
    }
  }
  return E_SUCCESS;
}

int check_for_valid_folding(struct extended_candidate *ecand) {
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence_list *css_list = ecand->possible_micro_rnas;
  struct candidate_subsequence *mature_mirna = NULL;
  for (size_t i = 0; i < css_list->n; i++) {
    mature_mirna = css_list->mature_sequences[i];
    if (mature_mirna == NULL) {
      return E_NO_MATURE_MI_RNA;
    }
    int loop_count = 0;
    get_loop_count_for_substructure(&loop_count, cand->structure,
                                    mature_mirna->start, mature_mirna->end);
    if (loop_count > 0) {
      mature_mirna->is_valid = 0;
    }
  }
  return E_SUCCESS;
}
int check_for_terminal_mismatches(struct extended_candidate *ecand) {
  const char MISMATCH_CHAR = '.';
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence_list *css_list = ecand->possible_micro_rnas;
  struct candidate_subsequence *mature_mirna = NULL;
  for (size_t i = 0; i < css_list->n; i++) {
    mature_mirna = css_list->mature_sequences[i];
    if (cand->structure[mature_mirna->start] == MISMATCH_CHAR &&
        cand->structure[mature_mirna->start + 1] == MISMATCH_CHAR) {
      mature_mirna->is_valid = 0;
    }
    if (cand->structure[mature_mirna->end - 2] == MISMATCH_CHAR &&
        cand->structure[mature_mirna->end - 1] == MISMATCH_CHAR) {
      mature_mirna->is_valid = 0;
    }
  }
  return E_SUCCESS;
}
int check_for_adjacent_mismatches(struct extended_candidate *ecand) {
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence_list *css_list = ecand->possible_micro_rnas;
  struct candidate_subsequence *mature_mirna = NULL;
  const char MISMATCH_CHAR = '.';

  for (size_t i = 0; i < css_list->n; i++) {
    mature_mirna = css_list->mature_sequences[i];
    int mismatch_count = 0;
    for (size_t i = mature_mirna->start; i < mature_mirna->end; i++) {
      if (cand->structure[i] == MISMATCH_CHAR) {
        mismatch_count++;
        if (mismatch_count >= 4) {
          mature_mirna->is_valid = 0;
          break;
        }
      } else {
        mismatch_count = 0;
      }
    }
  }
  return E_SUCCESS;
}

int find_star_micro_rna(struct candidate_subsequence **result,
                        struct micro_rna_candidate *cand,
                        struct candidate_subsequence *mature_mirna,
                        struct chrom_coverage *chrom_cov,
                        struct configuration_params *config) {
  const char MISMATCH_CHAR = '.';
  const int DICER_OFFSET = 2;

  u32 bracket_start = mature_mirna->start;
  while (cand->structure[bracket_start] == MISMATCH_CHAR) {
    bracket_start++;
    if (bracket_start > mature_mirna->end) {
      /* should never happen */
      return E_UNKNOWN;
    }
  }
  u32 bracket_end = mature_mirna->end - 1;
  while (cand->structure[bracket_end] == MISMATCH_CHAR) {
    bracket_end--;
    if (bracket_end < mature_mirna->start) {
      /* should never happen */
      return E_UNKNOWN;
    }
  }
  u32 star_start = 0;
  u32 star_end = 0;
  size_t n = cand->end - cand->start;
  find_matching_bracket_index(&star_end, bracket_start, cand->structure, n);
  find_matching_bracket_index(&star_start, bracket_end, cand->structure, n);

  star_start += DICER_OFFSET;
  star_end += DICER_OFFSET + 1;
  if (star_end > n) {
    star_end = n;
  }
  size_t l = star_end - star_start;
  if (l < config->min_duplex_length || l >= config->max_duplex_length) {
    return E_NO_STAR_MI_RNA_FOUND;
  }

  u32 *cov_list = chrom_cov->coverage_plus;
  if (cand->strand == '-') {
    cov_list = chrom_cov->coverage_minus;
  }

  u64 coverage = 0;
  get_coverage_in_range(&coverage, cov_list, star_start + cand->start,
                        star_end + cand->start);

  struct candidate_subsequence *star_micro_rna = NULL;
  create_candidate_subseqence(&star_micro_rna, star_start, star_end, coverage,
                              -1.0);
  star_micro_rna->is_star = 1;
  star_micro_rna->is_artificial = 1;
  *result = star_micro_rna;
  return E_SUCCESS;
}

int find_matching_bracket_index(u32 *result, u32 target, char *structure,
                                size_t n) {
  const char OPENING_BRACKET = '(';
  const char CLOSING_BRACKET = ')';
  int dir = 0;
  char initial_bracket;
  char inverse_bracket;
  if (structure[target] == OPENING_BRACKET) {
    initial_bracket = OPENING_BRACKET;
    inverse_bracket = CLOSING_BRACKET;
    dir = 1;
  }
  if (structure[target] == CLOSING_BRACKET) {
    initial_bracket = CLOSING_BRACKET;
    inverse_bracket = OPENING_BRACKET;
    dir = -1;
  }
  if (dir == 0) {
    return E_STRUCTURE_IS_INVALID;
  }
  u32 open_count = 0;
  u32 matching_index = 0;
  for (i64 i = target; i < n && i >= 0; i += dir) {
    char c = structure[i];
    if (c == initial_bracket) {
      open_count++;
    }
    if (c == inverse_bracket) {
      open_count--;
      if (open_count == 0) {
        matching_index = i;
        break;
      }
    }
  }
  if (open_count > 0) {
    return E_STRUCTURE_IS_INVALID;
  }
  *result = matching_index;

  return E_SUCCESS;
}
