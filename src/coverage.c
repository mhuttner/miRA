#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "coverage.h"
#include "parse_sam.h"
#include "errors.h"
#include "util.h"
#include "candidates.h"
#include "uthash.h"
#include "defs.h"
#include "structure_evaluation.h"
#include "reporting.h"

static int print_help();

int coverage(int argc, char **argv) {
  char *config_file = NULL;
  int c;
  int log_level = LOG_LEVEL_BASIC;
  char default_output_file[] = "output.tmp";
  char *output_file = default_output_file;

  while ((c = getopt(argc, argv, "c:o:hvq")) != -1) {
    switch (c) {
    case 'c':
      config_file = optarg;
      break;
    case 'o':
      output_file = optarg;
      break;
    case 'h':
      print_help();
      return E_SUCCESS;
    case 'v':
      log_level = LOG_LEVEL_VERBOSE;
    case 'q':
      log_level = LOG_LEVEL_QUIET;
    default:
      break;
    }
  }
  if (optind + 2 > argc) { /* missing input file */
    printf("No Input Files specified\n\n");
    print_help();
    return E_NO_FILE_SPECIFIED;
  }
  struct configuration_params *config = NULL;
  initialize_configuration(&config, config_file);
  log_configuration(config);
  int err;
  struct sam_file *sam = NULL;
  struct candidate_list *c_list = NULL;
  struct extended_candidate_list *ec_list = NULL;
  struct chrom_coverage *cov_table = NULL;

  err = read_candidate_file(&c_list, argv[optind]);
  if (err) {
    goto error;
  }
  err = parse_sam(&sam, argv[optind + 1]);
  if (err) {
    goto error;
  }
  err = create_coverage_table(&cov_table, sam);
  if (err) {
    goto error;
  }
  err = extend_all_candidates(&ec_list, c_list);
  if (err) {
    c_list = NULL;
    goto error;
  }
  err = coverage_test_candidates(ec_list, &cov_table, config);
  if (err) {
    goto error;
  }
  err = report_valid_candiates(ec_list);
  if (err) {
    goto error;
  }
  return E_SUCCESS;

error:
  if (ec_list != NULL) {
    free_extended_candidate_list(ec_list);
  }
  if (cov_table != NULL) {
    free_coverage_table(&cov_table);
  }
  if (sam != NULL) {
    free_sam(sam);
  }
  if (c_list != NULL) {
    free_candidate_list(c_list);
  }
  print_error(err);
  return err;
}

static int print_help() {
  printf("Description:\n"
         "    Coverage based verification on micro RNA candidates \n"
         "Usage: miRA coverage <input miRA file> <input SAM file>\n");
  return E_SUCCESS;
}

int create_coverage_table(struct chrom_coverage **table, struct sam_file *sam) {
  struct chrom_coverage *chrom_cov = NULL;
  struct sq_header *header = NULL;
  struct sam_entry *entry = NULL;
  for (size_t i = 0; i < sam->header_n; i++) {
    header = sam->headers[i];
    chrom_cov = (struct chrom_coverage *)malloc(sizeof(struct chrom_coverage));
    strncpy(chrom_cov->name, header->sn, 1024);
    chrom_cov->name[1023] = 0;
    chrom_cov->length = header->ln;
    chrom_cov->coverage_plus = (u32 *)malloc(chrom_cov->length * sizeof(u32));
    chrom_cov->coverage_minus = (u32 *)malloc(chrom_cov->length * sizeof(u32));
    for (u32 i = 0; i < chrom_cov->length; i++) {
      chrom_cov->coverage_plus[i] = 0;
      chrom_cov->coverage_minus[i] = 0;
    }
    HASH_ADD_STR(*table, name, chrom_cov);
  }
  u32 start;
  u32 stop;
  for (size_t i = 0; i < sam->n; i++) {
    entry = sam->entries[i];
    if (strcmp(chrom_cov->name, entry->rname) != 0) {
      HASH_FIND_STR(*table, entry->rname, chrom_cov);
    }
    if (chrom_cov == NULL) {
      return E_CHROMOSOME_NOT_FOUND;
    }

    start = entry->pos - 1;
    stop = start + strlen(entry->seq);
    u32 *cov_list = chrom_cov->coverage_plus;
    if (entry->flag & REV_COMPLM) {
      cov_list = chrom_cov->coverage_minus;
    }

    for (u32 i = start; i < stop; i++) {
      cov_list[i]++;
    }
  }

  return E_SUCCESS;
}
int coverage_test_candidates(struct extended_candidate_list *cand_list,
                             struct chrom_coverage **coverage_table,
                             struct configuration_params *config) {
  struct extended_candidate *ecand = NULL;
  struct micro_rna_candidate *cand = NULL;
  struct chrom_coverage *chrom_cov = NULL;

  int err = 0;
  for (size_t i = 0; i < cand_list->n; i++) {
    ecand = cand_list->candidates[i];
    ecand->is_valid = 0;
    cand = ecand->cand;
    printf("Cluster %lld\n", cand->id);
    HASH_FIND_STR(*coverage_table, cand->chrom, chrom_cov);
    if (chrom_cov == NULL) {
      return E_CHROMOSOME_NOT_FOUND;
    }
    err = find_mature_micro_rna(ecand, chrom_cov, config);
    if (err != E_SUCCESS || ecand->mature_micro_rna == NULL) {
      continue;
    }
    validate_mature_micro_rna(ecand, config);
    if (ecand->mature_micro_rna->is_valid == 0) {
      continue;
    }
    err = find_star_micro_rna(ecand, chrom_cov, config);
    if (err != E_SUCCESS || ecand->star_micro_rna == NULL) {
      continue;
    }
    ecand->is_valid = 1;
  }
  return E_SUCCESS;
}

int find_mature_micro_rna(struct extended_candidate *ecand,
                          struct chrom_coverage *chrom_cov,
                          struct configuration_params *config) {
  u32 *cov_list = chrom_cov->coverage_plus;
  struct micro_rna_candidate *cand = ecand->cand;
  if (cand->strand == '-') {
    cov_list = chrom_cov->coverage_minus;
  }

  size_t n = cand->end - cand->start;

  size_t neg_spike_count = 0;
  size_t pos_spike_count = 0;
  u64 total_coverage = 0;
  u64 best_start = 0;
  u64 best_end = 0;
  u64 best_coverage = 0;
  double best_paired_fraction = 0.0;
  struct candidate_subsequence *mature_micro_rna = NULL;

  get_coverage_in_range(&total_coverage, cov_list, cand->start, cand->end);
  int err = 0;
  u64 *pos_cov_spikes = (u64 *)malloc(n * sizeof(u64));
  u64 *neg_cov_spikes = (u64 *)malloc(n * sizeof(u64));
  if (pos_cov_spikes == NULL || neg_cov_spikes == NULL) {
    err = E_MALLOC_FAIL;
    goto error;
  }

  for (size_t i = 0; i < (n - 1); i++) {
    u64 global_index = cand->start + i;
    double cov_local = (double)cov_list[global_index] / total_coverage;
    double cov_next = (double)cov_list[global_index + 1] / total_coverage;
    if (cov_next - cov_local > config->min_coverage) {
      pos_cov_spikes[pos_spike_count] = i;
      pos_spike_count++;
    }
    if (cov_local - cov_next > config->min_coverage) {
      neg_cov_spikes[neg_spike_count] = i;
      neg_spike_count++;
    }
  }

  for (size_t i = 0; i < pos_spike_count; i++) {
    for (size_t j = 0; j < neg_spike_count; j++) {
      u64 start = pos_cov_spikes[i] + 1;
      u64 end = neg_cov_spikes[j];
      u64 length = end - start;
      if (length <= 0) {
        continue;
      }
      if (length < config->min_mature_strand_length) {
        continue;
      }
      if (length >= config->max_mature_strand_length) {
        continue;
      }
      double paired_fraction = 0.0;
      get_paired_fraction_for_substructure(&paired_fraction, cand->structure,
                                           start, end);
      if (paired_fraction < config->min_paired_fraction) {
        continue;
      }
      u64 segment_coverage = 0;
      get_coverage_in_range(&segment_coverage, cov_list, cand->start + start,
                            cand->start + end);
      if (segment_coverage > best_coverage) {
        best_coverage = segment_coverage;
        best_start = start;
        best_end = end;
        best_paired_fraction = paired_fraction;
      }
    }
  }
  free(pos_cov_spikes);
  free(neg_cov_spikes);
  pos_cov_spikes = NULL;
  neg_cov_spikes = NULL;
  if (best_coverage == 0) {
    err = E_NO_MATURE_MI_RNA_FOUND;
    goto error;
  }
  err = create_candidate_subseqence(&mature_micro_rna, best_start, best_end,
                                    best_coverage, best_paired_fraction);
  if (err != E_SUCCESS) {
    goto error;
  }
  ecand->mature_micro_rna = mature_micro_rna;
  return E_SUCCESS;
error:
  if (pos_cov_spikes != NULL) {
    free(pos_cov_spikes);
  }
  if (neg_cov_spikes != NULL) {
    free(pos_cov_spikes);
  }
  return err;
}

int find_star_micro_rna(struct extended_candidate *ecand,
                        struct chrom_coverage *chrom_cov,
                        struct configuration_params *config) {
  const char MISMATCH_CHAR = '.';
  const int DICER_OFFSET = 2;

  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence *mature_mirna = ecand->mature_micro_rna;

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
  if (l < config->min_mature_strand_length ||
      l >= config->max_mature_strand_length) {
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
  ecand->star_micro_rna = star_micro_rna;

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

int get_coverage_in_range(u64 *result, u32 *cov_list, u64 start, u64 end) {
  u64 total_coverage = 0;
  for (size_t i = start; i < end; i++) {
    total_coverage += cov_list[i];
  }
  *result = total_coverage;
  return E_SUCCESS;
}
int extend_all_candidates(struct extended_candidate_list **ecand_list,
                          struct candidate_list *cand_list) {
  struct extended_candidate_list *ecand_list_tmp =
      (struct extended_candidate_list *)malloc(
          sizeof(struct extended_candidate_list));
  if (ecand_list_tmp == NULL) {
    free_candidate_list(cand_list);
    return E_MALLOC_FAIL;
  }
  ecand_list_tmp->candidates = (struct extended_candidate **)malloc(
      cand_list->n * sizeof(struct extended_candidate *));
  if (ecand_list_tmp->candidates == NULL) {
    free_candidate_list(cand_list);
    free(ecand_list_tmp);
    return E_MALLOC_FAIL;
  }
  ecand_list_tmp->n = 0;
  struct extended_candidate *ecand = NULL;
  int err = 0;
  for (size_t i = 0; i < cand_list->n; i++) {
    err = create_extended_candidate(&ecand, cand_list->candidates[i]);
    if (err != E_SUCCESS) {
      free_extended_candidate_list(ecand_list_tmp);
      for (size_t j = i; j < cand_list->n; j++) {
        free_micro_rna_candidate(cand_list->candidates[i]);
      }
      free(cand_list->candidates);
      free(cand_list);
      return err;
    }
    ecand_list_tmp->candidates[i] = ecand;
    ecand_list_tmp->n++;
  }
  *ecand_list = ecand_list_tmp;
  free(cand_list->candidates);
  free(cand_list);
  return E_SUCCESS;
}

int validate_mature_micro_rna(struct extended_candidate *ecand,
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
int determine_mature_arm(struct extended_candidate *ecand) {
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence *mature_mirna = ecand->mature_micro_rna;

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

  return E_SUCCESS;
}

int check_for_valid_folding(struct extended_candidate *ecand) {
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence *mature_mirna = ecand->mature_micro_rna;
  if (mature_mirna == NULL) {
    return E_NO_MATURE_MI_RNA;
  }
  int loop_count = 0;
  get_loop_count_for_substructure(&loop_count, cand->structure,
                                  mature_mirna->start, mature_mirna->end);
  if (loop_count > 0) {
    mature_mirna->is_valid = 0;
  }
  return E_SUCCESS;
}
int check_for_terminal_mismatches(struct extended_candidate *ecand) {
  const char MISMATCH_CHAR = '.';
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence *mature_mirna = ecand->mature_micro_rna;
  if (cand->structure[mature_mirna->start] == MISMATCH_CHAR &&
      cand->structure[mature_mirna->start + 1] == MISMATCH_CHAR) {
    mature_mirna->is_valid = 0;
  }
  if (cand->structure[mature_mirna->end - 2] == MISMATCH_CHAR &&
      cand->structure[mature_mirna->end - 1] == MISMATCH_CHAR) {
    mature_mirna->is_valid = 0;
  }
  return E_SUCCESS;
}
int check_for_adjacent_mismatches(struct extended_candidate *ecand) {
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence *mature_mirna = ecand->mature_micro_rna;
  const char MISMATCH_CHAR = '.';
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
  return E_SUCCESS;
}

int create_extended_candidate(struct extended_candidate **ecand,
                              struct micro_rna_candidate *cand) {
  struct extended_candidate *ecand_tmp =
      (struct extended_candidate *)malloc(sizeof(struct extended_candidate));
  if (ecand_tmp == NULL) {
    return E_MALLOC_FAIL;
  }
  ecand_tmp->cand = cand;
  ecand_tmp->mature_micro_rna = NULL;
  ecand_tmp->star_micro_rna = NULL;
  *ecand = ecand_tmp;
  return E_SUCCESS;
}

int create_candidate_subseqence(struct candidate_subsequence **cs, u32 start,
                                u32 end, u64 coverage, double paired_fraction) {
  struct candidate_subsequence *cs_tmp = (struct candidate_subsequence *)malloc(
      sizeof(struct candidate_subsequence));
  if (cs_tmp == NULL) {
    return E_MALLOC_FAIL;
  }
  cs_tmp->start = start;
  cs_tmp->end = end;
  cs_tmp->coverage = coverage;
  cs_tmp->paired_fraction = paired_fraction;
  cs_tmp->is_valid = 1;
  *cs = cs_tmp;
  return E_SUCCESS;
}
int free_extended_candidate_list(struct extended_candidate_list *ecand_list) {
  for (size_t i = 0; i < ecand_list->n; i++) {
    free_extended_candidate(ecand_list->candidates[i]);
  }
  free(ecand_list->candidates);
  free(ecand_list);
  return E_SUCCESS;
}

int free_extended_candidate(struct extended_candidate *ecand) {
  free_micro_rna_candidate(ecand->cand);
  free_candidate_subsequence(ecand->mature_micro_rna);
  free_candidate_subsequence(ecand->star_micro_rna);
  free(ecand);
  return E_SUCCESS;
}

int free_candidate_subsequence(struct candidate_subsequence *cs) {
  if (cs != NULL) {
    free(cs);
  }
  return E_SUCCESS;
}

int free_coverage_table(struct chrom_coverage **table) {
  struct chrom_coverage *cov = NULL;
  struct chrom_coverage *tmp = NULL;

  HASH_ITER(hh, *table, cov, tmp) {
    HASH_DEL(*table, cov);
    free(cov->coverage_plus);
    free(cov->coverage_minus);
    free(cov);
  }
  return E_SUCCESS;
}
