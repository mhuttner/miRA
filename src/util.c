
#include "util.h"
#include "errors.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <stddef.h>

static void set_default_config(struct configuration_params *config) {
  config->log_level = LOG_LEVEL_BASIC;
  config->openmp_thread_count = 1;

  config->cluster_gap_size = 10;
  config->cluster_min_reads = 10;
  config->cluster_window_size = 200;
  config->cluster_max_length = 2000;

  config->max_precursor_length = 0;
  config->min_precursor_length = 50;
  config->min_mfe_per_nt = -0.4;
  config->max_hairpin_count = 4;
  config->min_double_strand_length = 20;
  config->permutation_count = 100;
  config->max_pvalue = 0.01;

  config->min_coverage = 0.0001;
  config->min_paired_fraction = 0.67;
  config->min_mature_strand_length = 18;
  config->max_mature_strand_length = 30;
  config->allow_three_mismatches = 0;
  config->allow_two_terminal_mismatches = 0;
  config->create_coverage_plots = 1;
  config->create_structure_images = 1;
  config->create_coverage_images = 1;
  config->cleanup_auxiliary_files = 1;
}

static int parse_config_file(struct configuration_params *config,
                             char *config_file) {
  const int MAXLINELENGTH = 1024;
  const char *integer_tokens[] = {
      "log_level",                     "openmp_thread_count",
      "cluster_gap_size",              "cluster_min_reads",
      "cluster_window_size",           "cluster_max_length",
      "max_precursor_length",          "min_precursor_length",
      "max_hairpin_count",             "min_double_strand_length",
      "permutation_count",             "min_mature_strand_length",
      "max_mature_strand_length",      "allow_three_mismatches",
      "allow_two_terminal_mismatches", "create_coverage_plots",
      "create_structure_images",       "create_coverage_images",
      "cleanup_auxiliary_files"};
  int integer_token_offsets[] = {
      (int)offsetof(struct configuration_params, log_level),
      (int)offsetof(struct configuration_params, openmp_thread_count),
      (int)offsetof(struct configuration_params, cluster_gap_size),
      (int)offsetof(struct configuration_params, cluster_min_reads),
      (int)offsetof(struct configuration_params, cluster_window_size),
      (int)offsetof(struct configuration_params, cluster_max_length),
      (int)offsetof(struct configuration_params, max_precursor_length),
      (int)offsetof(struct configuration_params, min_precursor_length),
      (int)offsetof(struct configuration_params, max_hairpin_count),
      (int)offsetof(struct configuration_params, min_double_strand_length),
      (int)offsetof(struct configuration_params, permutation_count),
      (int)offsetof(struct configuration_params, min_mature_strand_length),
      (int)offsetof(struct configuration_params, max_mature_strand_length),
      (int)offsetof(struct configuration_params, allow_three_mismatches),
      (int)offsetof(struct configuration_params, allow_two_terminal_mismatches),
      (int)offsetof(struct configuration_params, create_coverage_plots),
      (int)offsetof(struct configuration_params, create_structure_images),
      (int)offsetof(struct configuration_params, create_coverage_images),
      (int)offsetof(struct configuration_params, cleanup_auxiliary_files)};
  const int integer_token_count = 19;
  const char *double_tokens[] = {"min_mfe_per_nt", "max_pvalue", "min_coverage",
                                 "min_paired_fraction"};
  int double_token_offsets[] = {
      (int)offsetof(struct configuration_params, min_mfe_per_nt),
      (int)offsetof(struct configuration_params, max_pvalue),
      (int)offsetof(struct configuration_params, min_coverage),
      (int)offsetof(struct configuration_params, min_paired_fraction)};
  const int double_token_count = 4;

  const char COMMENT_CHAR = '#';
  FILE *fp = fopen(config_file, "r");
  char line[MAXLINELENGTH];
  if (fp == NULL) {
    return E_FILE_NOT_FOUND;
  }
  while (fgets(line, sizeof(line), fp) != NULL) {
    char *end = strchr(line, COMMENT_CHAR);
    if (end == NULL) {
      end = strchr(line, '\n');
      if (end == NULL) {
        strchr(line, '\0');
      }
      if (end == NULL) {
        end = line + MAXLINELENGTH - 1;
      }
      /* ignore everthing after comment */
      *end = 0;
      for (int i = 0; i < integer_token_count; i++) {
        char *match = strstr(line, integer_tokens[i]);
        if (match == NULL) {
          continue;
        }
        char *eq = strchr(match, '=');
        if (eq == NULL) {
          break;
        }

        char *tmp = NULL;
        long value = strtol(eq + 1, &tmp, 10);
        if (tmp == eq + 1 || tmp == NULL) {
          break;
        }
        int *target = (int *)((long)config + integer_token_offsets[i]);
        *target = (int)value;
        break;
      }
      for (int i = 0; i < double_token_count; i++) {
        char *match = strstr(line, double_tokens[i]);
        if (match == NULL) {
          continue;
        }
        char *eq = strchr(match, '=');
        if (eq == NULL) {
          break;
        }
        char *tmp = NULL;
        double value = strtod(eq + 1, &tmp);
        if (tmp == eq + 1 || tmp == NULL) {
          break;
        }
        double *target = (double *)((long)config + double_token_offsets[i]);
        *target = (double)value;
        break;
      }
    }
  }
  fclose(fp);
  return E_SUCCESS;
}

int initialize_configuration(struct configuration_params **config,
                             char *config_file) {
  struct configuration_params *tmp_config =
      (struct configuration_params *)malloc(
          sizeof(struct configuration_params));
  if (tmp_config == NULL) {
    return E_MALLOC_FAIL;
  }
  set_default_config(tmp_config);
  if (config_file != NULL) {
    parse_config_file(tmp_config, config_file);
  }
  *config = tmp_config;
  return E_SUCCESS;
}

int reverse_complement_sequence_string(char **result, char *seq, size_t n) {
  const char *pairs[] = {"AT", "GC", "UA", "YR", "SS",
                         "WW", "KM", "BV", "DH", "NN"};
  const int pairs_n = 10;
  char *tmp = (char *)malloc(n * sizeof(char));
  if (tmp == NULL) {
    return E_MALLOC_FAIL;
  }
  size_t k = n - 1;
  tmp[k--] = 0;
  for (size_t i = 0; i < n; i++) {
    for (int j = 0; j < pairs_n; j++) {
      if (seq[i] == pairs[j][0]) {
        tmp[k--] = pairs[j][1];
        break;
      }
      if (seq[i] == pairs[j][1]) {
        tmp[k--] = pairs[j][0];
        break;
      }
    }
  }
  *result = tmp;
  return E_SUCCESS;
};

int create_file_path(char **file_path, const char *path, const char *filename) {
  size_t path_n = strnlen(path, 1024);
  size_t file_n = strnlen(filename, 1024);
  if (path_n >= 1024 || file_n >= 1024) {
    return E_FILE_DESCRIPTOR_TOO_LONG;
  }
  char *file_path_tmp = (char *)malloc((path_n + file_n + 2) * sizeof(char));
  if (file_path_tmp == NULL) {
    return E_MALLOC_FAIL;
  }
  if (path[path_n - 1] == '/') {
    sprintf(file_path_tmp, "%s%s", path, filename);
  } else {
    sprintf(file_path_tmp, "%s/%s", path, filename);
  }
  *file_path = file_path_tmp;
  return E_SUCCESS;
}

void log_configuration(struct configuration_params *config) {
  log_basic(config->log_level, "Configuartion Parameters:\n");
  log_basic(config->log_level, "    log_level %d\n", config->log_level);
  log_basic(config->log_level, "    openmp_thread_count %d\n",
            config->openmp_thread_count);
  log_basic(config->log_level, "    cluster_gap_size %d\n",
            config->cluster_gap_size);
  log_basic(config->log_level, "    cluster_min_reads %d\n",
            config->cluster_min_reads);
  log_basic(config->log_level, "    cluster_window_size %d\n",
            config->cluster_window_size);
  log_basic(config->log_level, "    cluster_max_length %d\n",
            config->cluster_max_length);
  log_basic(config->log_level, "    max_precursor_length %d\n",
            config->max_precursor_length);
  log_basic(config->log_level, "    min_precursor_length %d\n",
            config->min_precursor_length);
  log_basic(config->log_level, "    min_mfe_per_nt %lf\n",
            config->min_mfe_per_nt);
  log_basic(config->log_level, "    max_hairpin_count %d\n",
            config->max_hairpin_count);
  log_basic(config->log_level, "    min_double_strand_length %d\n",
            config->min_double_strand_length);
  log_basic(config->log_level, "    permutation_count %d\n",
            config->permutation_count);
  log_basic(config->log_level, "    max_pvalue %lf\n", config->max_pvalue);
  log_basic(config->log_level, "    min_coverage %lf\n", config->min_coverage);
  log_basic(config->log_level, "    min_paired_fraction %lf\n",
            config->min_paired_fraction);
  log_basic(config->log_level, "    min_mature_strand_length %d\n",
            config->min_mature_strand_length);
  log_basic(config->log_level, "    max_mature_strand_length %d\n",
            config->max_mature_strand_length);
  log_basic(config->log_level, "    allow_three_mismatches %d\n",
            config->allow_three_mismatches);
  log_basic(config->log_level, "    allow_two_terminal_mismatches %d\n",
            config->allow_two_terminal_mismatches);
  log_basic(config->log_level, "    create_coverage_plots %d\n",
            config->create_coverage_plots);
  log_basic(config->log_level, "    create_structure_images %d\n",
            config->create_structure_images);
  log_basic(config->log_level, "    create_coverage_images %d\n",
            config->create_coverage_images);
  log_basic(config->log_level, "    cleanup_auxiliary_files %d\n",
            config->cleanup_auxiliary_files);
  log_basic(config->log_level, "\n\n");
}

void log_basic(int loglevel, const char *msg, ...) {
  if (loglevel < LOG_LEVEL_BASIC) {
    return;
  }
  va_list fmtargs;
  va_start(fmtargs, msg);
  vprintf(msg, fmtargs);
  va_end(fmtargs);
  fflush(stdout);
}

void log_verbose(int loglevel, const char *msg, ...) {
  if (loglevel < LOG_LEVEL_VERBOSE) {
    return;
  }
  va_list fmtargs;
  va_start(fmtargs, msg);
  vprintf(msg, fmtargs);
  va_end(fmtargs);
  fflush(stdout);
}

double mean(double *list, int n) {
  double sum = 0.0;
  if (n <= 0) {
    return sum;
  }
  for (int i = 0; i < n; i++) {
    sum += list[i];
  }
  return sum / (double)n;
}

double sd(double *list, int n, double mean) {
  double sd_sum = 0.0;
  if (n <= 1) {
    return sd_sum;
  }
  for (int i = 0; i < n; i++) {
    double tmp = list[i] - mean;
    sd_sum += tmp * tmp;
  }
  return sqrt(sd_sum / (double)(n - 1));
}

double pvalue(double mean, double sd, double value) {
  return erfc(fabs(value - mean) / (sqrt(2) * sd)) / 2.0;
}