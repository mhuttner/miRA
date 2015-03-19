#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "float.h"
#include "vfold.h"
#include "errors.h"
#include "fasta.h"
#include "cluster.h"
#include "defs.h"
#include "bed.h"
#include "util.h"
#include "Lfold/Lfold.h"
#include "Lfold/fold.h"
#include "structure_evaluation.h"
#include "candidates.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static int print_help();

int vfold(int argc, char *argv[]) {
  srand(time(NULL));
  int c;
  int log_level = LOG_LEVEL_BASIC;
  char default_output_file[] = "output";
  char *output_file = default_output_file;
  char *config_file = NULL;

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
  err = vfold_main(config, argv[optind], argv[optind + 1], output_file);
  free(config);
  return err;
}

int vfold_main(struct configuration_params *config, char *bed_file,
               char *fasta_file, char *output_file) {

#ifdef _OPENMP
  omp_set_num_threads(config->openmp_thread_count);
#endif

  struct cluster_list *c_list = NULL;
  struct genome_sequence *seq_table = NULL;
  struct sequence_list *seq_list = NULL;
  struct candidate_list *cand_list = NULL;
  char *json_output_file = NULL;
  char *mira_output_file = NULL;
  int err;
  err = read_bed_file(&c_list, bed_file);
  if (err) {
    goto bed_read_err;
  }
  err = read_fasta_file(&seq_table, fasta_file);
  if (err) {
    goto fasta_read_err;
  }
  err = map_clusters(&seq_list, c_list, seq_table);
  if (err) {
    goto map_error;
  }
  /*clusters freed by map_clusters */
  free_sequence_table(seq_table);

  err = fold_sequences(seq_list, config);
  if (err) {
    goto fold_error;
  }
  json_output_file = (char *)malloc((strlen(output_file) + 6) * sizeof(char));
  mira_output_file = (char *)malloc((strlen(output_file) + 6) * sizeof(char));
  if (json_output_file == NULL || mira_output_file == NULL) {
    err = E_MALLOC_FAIL;
    goto fold_error;
  }
  sprintf(json_output_file, "%s.json", output_file);
  sprintf(mira_output_file, "%s", output_file);
  err = write_json_result(seq_list, json_output_file);
  if (err) {
    goto write_error;
  }
  err = convert_seq_list_to_cand_list(&cand_list, seq_list);
  if (err) {
    goto convert_error;
  }
  err = write_candidate_file(cand_list, mira_output_file);
  if (err) {
    goto convert_error;
  }
  free(json_output_file);
  free(mira_output_file);

  free_sequence_list(seq_list);

  return E_SUCCESS;
convert_error:
write_error:
  free(json_output_file);
  free(mira_output_file);
fold_error:
  free_sequence_list(seq_list);
  print_error(err);
  return err;
map_error:
  free_sequence_table(seq_table);
  print_error(err);
  return err;
fasta_read_err:
  free_clusters(c_list);
bed_read_err:
  print_error(err);
  return err;
}

static int print_help() {
  printf("Description:\n"
         "    fold tries to fold rna sequences and calculates secondary\n"
         "    structure information \n"
         "Usage: miRA fold [-c config file] [-o output file] [-q] [-v] [-h] \n"
         "    <input BED file> <input FASTA file>\n"
         "\n"
         "Options:\n"
         "-c <config file>\n"
         "    pass a configuration file to the programm, containing parameters "

         );
  return E_SUCCESS;
}

int map_clusters(struct sequence_list **seq_list, struct cluster_list *c_list,
                 struct genome_sequence *seq_table) {
  struct sequence_list *tmp_seq_list =
      (struct sequence_list *)malloc(sizeof(struct sequence_list));
  if (tmp_seq_list == NULL) {
    return E_MALLOC_FAIL;
  }
  size_t n = c_list->n;
  tmp_seq_list->n = 0;
  tmp_seq_list->sequences = (struct foldable_sequence **)malloc(
      n * sizeof(struct foldable_sequence *));
  if (tmp_seq_list->sequences == NULL) {
    free(tmp_seq_list);
    return E_MALLOC_FAIL;
  }
  struct foldable_sequence *fs = NULL;
  struct cluster *c = NULL;
  struct genome_sequence *gs = NULL;
  for (size_t i = 0; i < n; i++) {
    c = c_list->clusters[i];
    HASH_FIND_STR(seq_table, c->chrom, gs);
    if (gs == NULL) {
      free_sequence_list(tmp_seq_list);
      return E_NEEDED_SEQUENCE_NOT_FOUND;
    }
    if (c->flank_end > gs->n + 1) {
      free_sequence_list(tmp_seq_list);
      return E_INVALID_FASTA_SEQUENCE_LENGTH;
    }
    fs = (struct foldable_sequence *)malloc(sizeof(struct foldable_sequence));
    fs->c = c;

    size_t l = c->flank_end - c->flank_start - 1;
    fs->seq = (char *)malloc((l + 1) * sizeof(char));
    if (fs->seq == NULL) {
      free_sequence_list(tmp_seq_list);
      return E_MALLOC_FAIL;
    }

    memcpy(fs->seq, gs->data + c->flank_start, l);
    fs->seq[l] = 0;
    fs->n = l + 1;

    fs->structure = NULL;
    if (c->strand == '-') {
      int err = reverse_complement(fs);
      if (err) {
        free_sequence_list(tmp_seq_list);
        return err;
      }
    }
    tmp_seq_list->n++;
    tmp_seq_list->sequences[i] = fs;
  }

  free(c_list->clusters);
  free(c_list);
  *seq_list = tmp_seq_list;

  return E_SUCCESS;
}
int fold_sequences(struct sequence_list *seq_list,
                   struct configuration_params *config) {
  struct foldable_sequence *fs = NULL;
  struct structure_list *s_list = NULL;
  size_t progress_count = 0;

  log_basic(config->log_level, "Initializing folding...\n");
#ifdef _OPENMP
#pragma omp parallel for private(fs, s_list) schedule(dynamic)
#endif
  for (size_t i = 0; i < seq_list->n; i++) {

#pragma omp critical
    {
      progress_count++;
      log_basic(config->log_level, "Folding sequence %5ld \\%5ld ... \n",
                progress_count, seq_list->n);
    }
    fs = seq_list->sequences[i];
    int max_length = fs->n;
    if (config->max_precursor_length > 0 &&
        config->max_precursor_length < max_length) {
      max_length = config->max_precursor_length;
    }
    Lfold(&s_list, fs->seq, max_length);
    find_optimal_structure(s_list, fs, config);
    free_structure_list(s_list);
    if (fs->structure == NULL) {
      continue;
    }
    evaluate_structure(fs->structure);

    check_folding_constraints(fs, config);
    if (fs->structure->is_valid == 0) {
      continue;
    }
    calculate_mfe_distribution(fs, config->permutation_count);
    check_pvalue(fs, config);
    if (fs->structure->is_valid == 0) {
      continue;
    }
    if (config->log_level == LOG_LEVEL_VERBOSE) {
      write_foldable_sequence(NULL, fs);
    }
  }
  log_basic(config->log_level, "Folding completed successfully.\n");
  return E_SUCCESS;
};

int write_json_result(struct sequence_list *seq_list, char *filename) {
  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    return E_UNKNOWN_FILE_IO_ERROR;
  }
  fprintf(fp, "{\n");
  for (size_t i = 0; i < seq_list->n - 1; i++) {
    if (seq_list->sequences[i]->structure == NULL) {
      continue;
    }
    if (seq_list->sequences[i]->structure->is_valid == 0) {
      continue;
    }
    write_foldable_sequence(fp, seq_list->sequences[i]);
    fprintf(fp, ",\n");
  }
  write_foldable_sequence(fp, seq_list->sequences[seq_list->n - 1]);
  fprintf(fp, "}");
  fclose(fp);
  return E_SUCCESS;
}

int calculate_mfe_distribution(struct foldable_sequence *fs,
                               int permutation_count) {
  if (fs->structure == NULL) {
    return E_NO_STRUCTURE;
  }
  char *seq_copy = (char *)malloc((fs->n) * sizeof(char));
  if (seq_copy == NULL) {
    return E_MALLOC_FAIL;
  }
  double *mfe_list = (double *)malloc(permutation_count * sizeof(double));
  if (mfe_list == NULL) {
    free(seq_copy);
    return E_MALLOC_FAIL;
  }
  char *tmp = (char *)malloc((fs->n) * sizeof(char));
  if (tmp == NULL) {
    free(seq_copy);
    free(mfe_list);
    return E_MALLOC_FAIL;
  }

  memcpy(seq_copy, fs->seq, fs->n);
  seq_copy[fs->n - 1] = 0;

  for (int i = 0; i < permutation_count; i++) {
    fisher_yates_shuffle(seq_copy, fs->n - 1);
    mfe_list[i] = fold(seq_copy, tmp) / fs->n;
  }

  struct structure_info *si = fs->structure;
  double mfe_mean = mean(mfe_list, permutation_count);
  double mfe_sd = sd(mfe_list, permutation_count, mfe_mean);
  double mfe_pvalue = pvalue(mfe_mean, mfe_sd, si->mfe);
  si->mean = mfe_mean;
  si->sd = mfe_sd;
  si->pvalue = mfe_pvalue;
  free(mfe_list);
  free(tmp);
  free(seq_copy);

  return E_SUCCESS;
}

int find_optimal_structure(struct structure_list *s_list,
                           struct foldable_sequence *fs,
                           struct configuration_params *config) {
  struct secondary_structure *ss = NULL;
  struct secondary_structure *best_ss = NULL;
  double min_mfe = DBL_MAX;
  struct cluster *c = fs->c;
  u64 local_core_start = c->start - c->flank_start;
  /* offset du to bed file format */
  u64 local_core_end = c->end - c->flank_start - 1;
  for (size_t i = 0; i < s_list->n; i++) {
    ss = s_list->structures[i];
    /*the string should always be null terminated, if not the library would have
     * failed already */
    size_t l = strlen(ss->structure_string);
    if (ss->start > local_core_start + 1) {
      continue;
    }
    if (ss->start + l - 1 < local_core_end + 1) {
      continue;
    }
    double mfe_per_bp = ss->mfe / (double)l;
    if (mfe_per_bp < min_mfe) {
      best_ss = ss;
      min_mfe = mfe_per_bp;
    }
  }
  if (best_ss == NULL) {
    fs->structure = NULL;
    return E_NO_OPTIMAL_STRUCTURE_FOUND;
  }
  size_t n = strlen(best_ss->structure_string);
  fs->structure =
      (struct structure_info *)malloc(sizeof(struct structure_info));
  if (fs->structure == NULL) {
    return E_MALLOC_FAIL;
  }

  char *st_string = (char *)malloc((n + 1) * sizeof(char));
  if (st_string == NULL) {
    free(fs->structure);
    return E_MALLOC_FAIL;
  }

  memcpy(st_string, best_ss->structure_string, n);
  st_string[n] = 0;
  fs->structure->structure_string = st_string;
  fs->structure->n = n;
  fs->structure->start = best_ss->start;
  fs->structure->mfe = min_mfe;

  return E_SUCCESS;
}

int check_folding_constraints(struct foldable_sequence *fs,
                              struct configuration_params *config) {
  struct structure_info *si = fs->structure;
  if (si->n < config->min_precursor_length) {
    goto invalid_structure;
  }
  if (si->external_loop_count >= config->max_hairpin_count) {
    goto invalid_structure;
  }
  if (abs(si->stem_end_with_mismatch - si->stem_start_with_mismatch) + 1 <
      config->min_double_strand_length) {
    goto invalid_structure;
  }
  if (si->mfe >= config->min_mfe_per_nt) {
    goto invalid_structure;
  }
  si->is_valid = 1;
  return E_SUCCESS;

invalid_structure:
  si->is_valid = 0;
  return E_STRUCTURE_IS_INVALID;
}
int check_pvalue(struct foldable_sequence *fs,
                 struct configuration_params *config) {
  if (fs->structure->pvalue > config->max_pvalue) {
    fs->structure->is_valid = 0;
    return E_STRUCTURE_IS_INVALID;
  }
  return E_SUCCESS;
}

int write_foldable_sequence(FILE *fp, struct foldable_sequence *fs) {
  if (fp == NULL) {
    fp = stdout;
  }
  struct cluster *c = fs->c;
  struct structure_info *si = fs->structure;
  if (si == NULL) {
    return E_NO_STRUCTURE;
  }

  long structure_start;
  long structure_end;
  if (c->strand == '+') {
    /*0 based */
    structure_start = si->start - 1 + c->flank_start;
    structure_end = structure_start + si->n;
  } else {
    /*0 based */
    structure_end = c->flank_end - (si->start - 1);
    structure_start = structure_end - si->n;
  }
  char *structure_sequence = (char *)malloc((si->n + 1) * sizeof(char));
  if (structure_sequence == NULL) {
    return E_MALLOC_FAIL;
  }
  memcpy(structure_sequence, fs->seq + si->start - 1, si->n);
  structure_sequence[si->n] = 0;

  fprintf(fp, "\"Cluster_%lld_%s\":{\n", c->id,
          (c->strand == '-') ? "minus" : "plus");
  fprintf(fp, "\t\"chromosome\":\"%s\",\n", c->chrom);
  fprintf(fp, "\t\"structure_start\":%ld,\n", structure_start);
  fprintf(fp, "\t\"structure_end\":%ld,\n", structure_end);
  fprintf(fp, "\t\"structure_length\":%ld,\n", si->n);
  fprintf(fp, "\t\"strand\":\"%c\",\n", c->strand);
  fprintf(fp, "\t\"sequence\":\"%s\",\n", structure_sequence);
  fprintf(fp, "\t\"structure\":\"%s\",\n", si->structure_string);
  fprintf(fp, "\t\"mfe\":\"%7.5lf kcal/mol/bp\",\n", si->mfe);
  fprintf(fp, "\t\"external_loops\":%d,\n", si->external_loop_count);
  fprintf(fp, "\t\"longest_stem_0mm\":\"%d (nt %d-%d)\",\n",
          si->stem_end - si->stem_start + 1, si->stem_start, si->stem_end);
  fprintf(fp, "\t\"longest_stem_1mm\":\"%d (nt %d-%d)\",\n",
          si->stem_end_with_mismatch - si->stem_start_with_mismatch + 1,
          si->stem_start_with_mismatch, si->stem_end_with_mismatch);
  fprintf(fp, "\t\"paired_fraction\":%7.5lf,\n", si->paired_fraction);
  fprintf(fp, "\t\"mfe_mean\":%7.5lf,\n", si->mean);
  fprintf(fp, "\t\"mfe_sd\":%7.5lf,\n", si->sd);
  fprintf(fp, "\t\"mfe_precise\":%8lf,\n", si->mfe);
  fprintf(fp, "\t\"z_value\":%7.5lf,\n", (si->mfe - si->mean) / si->sd);
  fprintf(fp, "\t\"p_value\":%7.5le,\n", si->pvalue);
  fprintf(fp, "}");
  return E_SUCCESS;
}

int reverse_complement(struct foldable_sequence *s) {
  char *reversed = NULL;
  int err = reverse_complement_sequence_string(&reversed, s->seq, s->n);
  if (err) {
    return err;
  }
  free(s->seq);
  s->seq = reversed;
  return E_SUCCESS;
}

/* Generates a uniformly distributed random integer in the range [0,max) */
int get_rand_int(int max) {
  /* make sure the random the possible amount of random numbers is a mutiple of
   * max to get a uniform distribution */
  int limit = RAND_MAX - RAND_MAX % max;
  int random_number;
  do {
    random_number = rand();
  } while (random_number >= limit);
  return random_number % max;
}

int fisher_yates_shuffle(char *seq, int n) {
  int j;
  char tmp;
  for (int i = n - 1; i > 0; i--) {
    j = get_rand_int(i + 1);
    tmp = seq[j];
    seq[j] = seq[i];
    seq[i] = tmp;
  }
  return E_SUCCESS;
}

int free_structure_info(struct structure_info *info) {
  free(info->structure_string);
  free(info);
  return E_SUCCESS;
}

int free_sequence_list(struct sequence_list *seq_list) {
  for (size_t i = 0; i < seq_list->n; i++) {
    free_foldable_sequence(seq_list->sequences[i]);
  }
  free(seq_list);
  return E_SUCCESS;
}

int free_foldable_sequence(struct foldable_sequence *s) {
  if (s == NULL) {
    return E_SUCCESS;
  }
  free_cluster(s->c);
  free(s->seq);
  if (s->structure != NULL) {
    free_structure_info(s->structure);
  }
  free(s);
  return E_SUCCESS;
}
