#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "coverage.h"
#include "parse_sam.h"
#include "errors.h"
#include "util.h"
#include "candidates.h"
#include "uthash.h"
#include "defs.h"
#include "reporting.h"
#include "reads.h"
#include "structure_evaluation.h"
#include "mirna_validation.h"

static int print_help();

int coverage(int argc, char **argv) {
  char *config_file = NULL;
  int c;
  int log_level = LOG_LEVEL_BASIC;

  while ((c = getopt(argc, argv, "c:o:hvq")) != -1) {
    switch (c) {
    case 'c':
      config_file = optarg;
      break;
    case 'h':
      print_help();
      return E_SUCCESS;
    case 'v':
      log_level = LOG_LEVEL_VERBOSE;
      break;
    case 'q':
      log_level = LOG_LEVEL_QUIET;
      break;
    default:
      break;
    }
  }
  if (optind + 3 > argc) { /* missing input file */
    printf("No Input Files specified\n\n");
    print_help();
    return E_NO_FILE_SPECIFIED;
  }
  struct configuration_params *config = NULL;
  initialize_configuration(&config, config_file);
  if (log_level != LOG_LEVEL_BASIC) {
    config->log_level = log_level;
  }
  log_configuration(config);
  int err = coverage_main(config, argv[-1], argv[optind], argv[optind + 1],
                          argv[optind + 2], NULL);
  free(config);
  return err;
}

static int print_help() {
  printf("Description:\n"
         "    Coverage based verification on micro RNA candidates \n"
         "Usage: miRA coverage [-c config file] [-o output file] [-q] [-v] \n"
         "    [-h] <input miRA file> <input SAM file> <output directory>\n");
  return E_SUCCESS;
}

int coverage_main(struct configuration_params *config, char *executable_file,
                  char *mira_file, char *sam_file, char *output_path,
                  char *selected_crom) {
  int err;
  struct sam_file *sam = NULL;
  struct candidate_list *c_list = NULL;
  struct extended_candidate_list *ec_list = NULL;
  struct chrom_coverage *cov_table = NULL;
  log_basic_timestamp(config->log_level, "Coverage based verification...\n");

  for (size_t i = strlen(executable_file); i > 0; i--) {
    if (executable_file[i] == '/') {
      executable_file[i] = 0;
      break;
    }
  }
  log_verbose_timestamp(config->log_level, "\tReading candidate file...\n");
  err = read_candidate_file(&c_list, mira_file);
  if (err) {
    goto error;
  }
  log_verbose_timestamp(config->log_level, "\tParsing sam file...\n");
  err = parse_sam(&sam, sam_file, selected_crom);
  if (err) {
    goto error;
  }
  log_verbose_timestamp(config->log_level, "\tCreating coverage table...\n");
  err = create_coverage_table(&cov_table, sam);
  if (err) {
    goto error;
  }
  log_verbose_timestamp(config->log_level, "\tExtending candidates...\n");
  err = extend_all_candidates(&ec_list, c_list);
  if (err) {
    c_list = NULL;
    goto error;
  }
  log_verbose_timestamp(config->log_level,
                        "\tCoverage testing candidates...\n");
  err = coverage_test_candidates(ec_list, &cov_table, sam, config);
  if (err) {
    goto error;
  }
  log_verbose_timestamp(config->log_level, "\tFreeing Sam file...\n");

  free_sam(sam);
  log_verbose_timestamp(config->log_level, "\tAll OK\n");
  log_basic_timestamp(config->log_level,
                      "Coverage based verification completed\n");
  log_basic_timestamp(config->log_level, "Generating reports...\n");
  err = report_valid_candiates(ec_list, &cov_table, executable_file,
                               output_path, config);
  if (err) {
    goto error;
  }
  log_basic_timestamp(config->log_level, "Generating reports completed\n");
  free_coverage_table(&cov_table);
  free_extended_candidate_list(ec_list);
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
                             struct sam_file *sam,
                             struct configuration_params *config) {
  struct extended_candidate *ecand = NULL;
  struct micro_rna_candidate *cand = NULL;
  struct chrom_coverage *chrom_cov = NULL;

  int err = 0;
  for (size_t i = 0; i < cand_list->n; i++) {
    ecand = cand_list->candidates[i];
    ecand->is_valid = 0;
    cand = ecand->cand;
    HASH_FIND_STR(*coverage_table, cand->chrom, chrom_cov);
    if (chrom_cov == NULL) {
      return E_CHROMOSOME_NOT_FOUND;
    }

    err = find_mature_micro_rnas(ecand, chrom_cov, config);
    if (err != E_SUCCESS) {
      continue;
    }

    validate_mature_micro_rnas(ecand, config);

    err = filter_mature_micro_rnas(ecand, chrom_cov, config);
    if (err != E_SUCCESS) {
      continue;
    }

    err = count_unique_reads(ecand, sam);
    if (err != E_SUCCESS) {
      continue;
    }
    ecand->is_valid = 1;
  }

  return E_SUCCESS;
}

int find_mature_micro_rnas(struct extended_candidate *ecand,
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

  struct candidate_subsequence *mature_micro_rna = NULL;

  get_coverage_in_range(&total_coverage, cov_list, cand->start, cand->end);
  int err = 0;

  err = create_candidate_subseqence_list(&(ecand->possible_micro_rnas));
  u64 *pos_cov_spikes = (u64 *)malloc(n * sizeof(u64));
  u64 *neg_cov_spikes = (u64 *)malloc(n * sizeof(u64));
  if (pos_cov_spikes == NULL || neg_cov_spikes == NULL || err) {
    err = E_MALLOC_FAIL;
    goto error;
  }

  for (size_t i = 0; i < (n - 1); i++) {
    u64 global_index = cand->start + i;
    u32 cov_local = cov_list[global_index];
    u32 cov_next = cov_list[global_index + 1];
    if (cov_next > cov_local) {
      pos_cov_spikes[pos_spike_count] = i;
      pos_spike_count++;
    }
    if (cov_local > cov_next) {
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
      if (length < config->min_duplex_length) {
        continue;
      }
      if (length >= config->max_duplex_length) {
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
      err = create_candidate_subseqence(&mature_micro_rna, start, end,
                                        segment_coverage, paired_fraction);
      if (err) {
        goto error;
      }
      err = append_candidate_subseqence_list(ecand->possible_micro_rnas,
                                             mature_micro_rna);
      if (err) {
        goto error;
      }
    }
  }
  free(pos_cov_spikes);
  free(neg_cov_spikes);
  pos_cov_spikes = NULL;
  neg_cov_spikes = NULL;
  return E_SUCCESS;
error:
  if (ecand->possible_micro_rnas != NULL) {
    free_candidate_subsequence_list(ecand->possible_micro_rnas);
  }
  if (pos_cov_spikes != NULL) {
    free(pos_cov_spikes);
  }
  if (neg_cov_spikes != NULL) {
    free(pos_cov_spikes);
  }
  return err;
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

int create_extended_candidate(struct extended_candidate **ecand,
                              struct micro_rna_candidate *cand) {
  struct extended_candidate *ecand_tmp =
      (struct extended_candidate *)malloc(sizeof(struct extended_candidate));
  if (ecand_tmp == NULL) {
    return E_MALLOC_FAIL;
  }
  ecand_tmp->cand = cand;
  ecand_tmp->possible_micro_rnas = NULL;
  *ecand = ecand_tmp;
  return E_SUCCESS;
}

int create_candidate_subseqence_list(
    struct candidate_subsequence_list **css_list) {
  const int INITIAL_CAPACITY = 32;
  struct candidate_subsequence_list *tmp_list =
      (struct candidate_subsequence_list *)malloc(
          sizeof(struct candidate_subsequence_list));
  if (tmp_list == NULL) {
    return E_MALLOC_FAIL;
  }
  tmp_list->capacity = INITIAL_CAPACITY;
  tmp_list->n = 0;
  tmp_list->mature_sequences = (struct candidate_subsequence **)malloc(
      tmp_list->capacity * sizeof(struct candidate_subsequence *));
  if (tmp_list->mature_sequences == NULL) {
    free(tmp_list);
    return E_MALLOC_FAIL;
  }
  *css_list = tmp_list;
  return E_SUCCESS;
}
int append_candidate_subseqence_list(
    struct candidate_subsequence_list *css_list,
    struct candidate_subsequence *mature) {
  if (css_list->n >= css_list->capacity) {
    css_list->capacity *= 2;
    struct candidate_subsequence **tmp =
        (struct candidate_subsequence **)realloc(
            css_list->mature_sequences,
            css_list->capacity * sizeof(struct candidate_subsequence *));
    if (tmp == NULL) {
      return E_REALLOC_FAIL;
    }
    css_list->mature_sequences = tmp;
  }
  css_list->mature_sequences[css_list->n] = mature;
  css_list->n++;
  return E_SUCCESS;
}
int free_candidate_subsequence_list(
    struct candidate_subsequence_list *css_list) {
  for (size_t i = 0; i < css_list->n; i++) {
    free_candidate_subsequence(css_list->mature_sequences[i]);
  }
  free(css_list->mature_sequences);
  free(css_list);
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
  cs_tmp->is_star = 0;
  cs_tmp->is_artificial = 0;
  cs_tmp->reads = NULL;
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
  free_candidate_subsequence_list(ecand->possible_micro_rnas);
  free(ecand);
  return E_SUCCESS;
}

int free_candidate_subsequence(struct candidate_subsequence *cs) {
  if (cs->reads != NULL) {
    free_unique_read_list(cs->reads);
  }
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
