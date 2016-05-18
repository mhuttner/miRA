
#include "util.h"
#include "errors.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <stddef.h>
#include <time.h>
#include <string.h>

int create_text_buffer(struct text_buffer **buffer) {
  const size_t INITIAL_SIZE = 4096;
  struct text_buffer *tmp_buffer =
      (struct text_buffer *)malloc(sizeof(struct text_buffer));
  if (tmp_buffer == NULL) {
    return E_MALLOC_FAIL;
  }
  tmp_buffer->capacity = INITIAL_SIZE;
  tmp_buffer->data = (char *)malloc(tmp_buffer->capacity * sizeof(char));
  if (tmp_buffer->data == NULL) {
    free(tmp_buffer);
    return E_MALLOC_FAIL;
  }
  tmp_buffer->free_space = tmp_buffer->capacity;
  tmp_buffer->start = tmp_buffer->data;

  *buffer = tmp_buffer;
  return E_SUCCESS;
}
void print_to_text_buffer(struct text_buffer *buffer, const char *msg, ...) {
  // temporary disabled because of memory error (see issue #1)
  return;
  if (buffer == NULL) {
    return;
  }
  if (buffer->free_space <= 0) {
    return;
  }
  int n;
  va_list fmtargs;
  va_start(fmtargs, msg);
  n = vsnprintf(buffer->data, buffer->free_space, msg, fmtargs);
  if (n > 0) {
    if (buffer->free_space < n) {
      buffer->capacity *= 2;
      char *tmp_buffer =
          (char *)realloc(buffer->start, buffer->capacity * sizeof(char));
      if (tmp_buffer == NULL) {
        buffer->free_space = 0;
      } else {
        buffer->free_space += buffer->capacity;
        long offset = buffer->data - buffer->start;
        buffer->start = tmp_buffer;
        buffer->data = buffer->start + offset;
      }
    } else {
      buffer->data += n;
      buffer->free_space -= n;
    }
  }
  va_end(fmtargs);
}
int free_text_buffer(struct text_buffer *buffer) {
  free(buffer->start);
  buffer->data = NULL;
  free(buffer);
  return E_SUCCESS;
}

static void set_default_config(struct configuration_params *config) {
  config->log_level = LOG_LEVEL_BASIC;
  config->openmp_thread_count = 1;

  config->cluster_gap_size = 10;
  config->cluster_min_reads = 10;
  config->cluster_flank_size = 200;
  config->cluster_max_length = 2000;

  config->max_precursor_length = 0;
  config->min_precursor_length = 50;
  config->max_mfe_per_nt = -0.4;
  config->max_hairpin_count = 4;
  config->min_double_strand_length = 20;
  config->permutation_count = 100;
  config->max_pvalue = 0.01;

  config->min_dicer_offset = 0;
  config->max_dicer_offset = 3;

  config->min_coverage = 0.0001;
  config->min_paired_fraction = 0.67;
  config->min_duplex_length = 18;
  config->max_duplex_length = 30;
  config->allow_three_mismatches = 0;
  config->allow_two_terminal_mismatches = 0;
  config->create_coverage_plots = 1;
  config->create_structure_plots = 1;
  config->create_structure_coverage_plots = 1;
  config->cleanup_auxiliary_files = 1;
}

static void set_algae_config(struct configuration_params *config) {
  config->log_level = LOG_LEVEL_BASIC;
  config->openmp_thread_count = 10;

  config->cluster_gap_size = 10;
  config->cluster_min_reads = 10;
  config->cluster_flank_size = 200;
  config->cluster_max_length = 2000;

  config->max_precursor_length = 0;
  config->min_precursor_length = 0;
  config->max_mfe_per_nt = -0.3;
  config->max_hairpin_count = 4;
  config->min_double_strand_length = 18;
  config->permutation_count = 100;
  config->max_pvalue = 0.01;

  config->min_coverage = 0.01;
  config->min_paired_fraction = 0.55;
  config->min_duplex_length = 18;
  config->max_duplex_length = 30;
  config->allow_three_mismatches = 1;
  config->allow_two_terminal_mismatches = 0;
}
static void set_plant_config(struct configuration_params *config) {
  config->log_level = LOG_LEVEL_BASIC;
  config->openmp_thread_count = 10;

  config->cluster_gap_size = 10;
  config->cluster_min_reads = 10;
  config->cluster_flank_size = 200;
  config->cluster_max_length = 2000;

  config->max_precursor_length = 0;
  config->min_precursor_length = 0;
  config->max_mfe_per_nt = -0.2;
  config->max_hairpin_count = 4;
  config->min_double_strand_length = 18;
  config->permutation_count = 100;
  config->max_pvalue = 0.01;

  config->min_coverage = 0.01;
  config->min_paired_fraction = 0.55;
  config->min_duplex_length = 18;
  config->max_duplex_length = 30;
  config->allow_three_mismatches = 1;
  config->allow_two_terminal_mismatches = 0;
}
static void set_animal_config(struct configuration_params *config) {
  config->log_level = LOG_LEVEL_BASIC;
  config->openmp_thread_count = 10;

  config->cluster_gap_size = 10;
  config->cluster_min_reads = 10;
  config->cluster_flank_size = 200;
  config->cluster_max_length = 2000;

  config->max_precursor_length = 200;
  config->min_precursor_length = 0;
  config->max_mfe_per_nt = -0.2;
  config->max_hairpin_count = 2;
  config->min_double_strand_length = 17;
  config->permutation_count = 100;
  config->max_pvalue = 0.01;

  config->min_coverage = 0.01;
  config->min_paired_fraction = 0.50;
  config->min_duplex_length = 18;
  config->max_duplex_length = 30;
  config->allow_three_mismatches = 1;
  config->allow_two_terminal_mismatches = 0;
}

static int parse_config_file(struct configuration_params *config,
                             char *config_file) {
  const int MAXLINELENGTH = 1024;
  const char *integer_tokens[] = {
      "log_level", "openmp_thread_count", "cluster_gap_size",
      "cluster_min_reads", "cluster_flank_size", "cluster_max_length",
      "max_precursor_length", "min_precursor_length", "max_hairpin_count",
      "min_double_strand_length", "permutation_count", "min_duplex_length",
      "max_duplex_length", "allow_three_mismatches",
      "allow_two_terminal_mismatches", "min_dicer_offset", "max_dicer_offset",
      "create_coverage_plots", "create_structure_plots",
      "create_structure_coverage_plots", "cleanup_auxiliary_files"};
  int integer_token_offsets[] = {
      (int)offsetof(struct configuration_params, log_level),
      (int)offsetof(struct configuration_params, openmp_thread_count),
      (int)offsetof(struct configuration_params, cluster_gap_size),
      (int)offsetof(struct configuration_params, cluster_min_reads),
      (int)offsetof(struct configuration_params, cluster_flank_size),
      (int)offsetof(struct configuration_params, cluster_max_length),
      (int)offsetof(struct configuration_params, max_precursor_length),
      (int)offsetof(struct configuration_params, min_precursor_length),
      (int)offsetof(struct configuration_params, max_hairpin_count),
      (int)offsetof(struct configuration_params, min_double_strand_length),
      (int)offsetof(struct configuration_params, permutation_count),
      (int)offsetof(struct configuration_params, min_duplex_length),
      (int)offsetof(struct configuration_params, max_duplex_length),
      (int)offsetof(struct configuration_params, allow_three_mismatches),
      (int)offsetof(struct configuration_params, allow_two_terminal_mismatches),
      (int)offsetof(struct configuration_params, min_dicer_offset),
      (int)offsetof(struct configuration_params, max_dicer_offset),
      (int)offsetof(struct configuration_params, create_coverage_plots),
      (int)offsetof(struct configuration_params, create_structure_plots),
      (int)offsetof(struct configuration_params,
                    create_structure_coverage_plots),
      (int)offsetof(struct configuration_params, cleanup_auxiliary_files)};
  const int integer_token_count = 21;
  const char *double_tokens[] = {"max_mfe_per_nt", "max_pvalue", "min_coverage",
                                 "min_paired_fraction"};
  int double_token_offsets[] = {
      (int)offsetof(struct configuration_params, max_mfe_per_nt),
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
    if (strncasecmp(config_file, "animal", 7) == 0) {
      set_animal_config(tmp_config);
      *config = tmp_config;
      return E_SUCCESS;
    }
    if (strncasecmp(config_file, "plant", 6) == 0) {
      set_plant_config(tmp_config);
      *config = tmp_config;
      return E_SUCCESS;
    }
    if (strncasecmp(config_file, "algae", 6) == 0) {
      set_algae_config(tmp_config);
      *config = tmp_config;
      return E_SUCCESS;
    }
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
  printf("log lev %d\n", config->log_level);
  log_basic(config->log_level, "Configuartion Parameters:\n");
  log_basic(config->log_level, "    log_level %d\n", config->log_level);
  log_basic(config->log_level, "    openmp_thread_count %d\n",
            config->openmp_thread_count);
  log_basic(config->log_level, "    cluster_gap_size %d\n",
            config->cluster_gap_size);
  log_basic(config->log_level, "    cluster_min_reads %d\n",
            config->cluster_min_reads);
  log_basic(config->log_level, "    cluster_flank_size %d\n",
            config->cluster_flank_size);
  log_basic(config->log_level, "    cluster_max_length %d\n",
            config->cluster_max_length);
  log_basic(config->log_level, "    max_precursor_length %d\n",
            config->max_precursor_length);
  log_basic(config->log_level, "    min_precursor_length %d\n",
            config->min_precursor_length);
  log_basic(config->log_level, "    max_mfe_per_nt %lf\n",
            config->max_mfe_per_nt);
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
  log_basic(config->log_level, "    min_duplex_length %d\n",
            config->min_duplex_length);
  log_basic(config->log_level, "    max_duplex_length %d\n",
            config->max_duplex_length);
  log_basic(config->log_level, "    allow_three_mismatches %d\n",
            config->allow_three_mismatches);
  log_basic(config->log_level, "    allow_two_terminal_mismatches %d\n",
            config->allow_two_terminal_mismatches);
  log_basic(config->log_level, "    min_dicer_offset %d\n",
            config->min_dicer_offset);
  log_basic(config->log_level, "    max_dicer_offset %d\n",
            config->max_dicer_offset);
  log_basic(config->log_level, "    create_coverage_plots %d\n",
            config->create_coverage_plots);
  log_basic(config->log_level, "    create_structure_plots %d\n",
            config->create_structure_plots);
  log_basic(config->log_level, "    create_structure_coverage_plots %d\n",
            config->create_structure_coverage_plots);
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
void log_basic_timestamp(int loglevel, const char *msg, ...) {
  if (loglevel < LOG_LEVEL_BASIC) {
    return;
  }
  time_t ltime;
  struct tm result;
  char stime[32];
  ltime = time(NULL);
  localtime_r(&ltime, &result);
  asctime_r(&result, stime);
  /* remove newline */
  stime[strlen(stime) - 1] = '\0';
  printf("%s --- ", stime);

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
void log_verbose_timestamp(int loglevel, const char *msg, ...) {
  if (loglevel < LOG_LEVEL_VERBOSE) {
    return;
  }
  time_t ltime;
  struct tm result;
  char stime[32];
  ltime = time(NULL);
  localtime_r(&ltime, &result);
  asctime_r(&result, stime);
  /* remove newline */
  stime[strlen(stime) - 1] = '\0';
  printf("%s --- ", stime);
  va_list fmtargs;
  va_start(fmtargs, msg);
  vprintf(msg, fmtargs);
  va_end(fmtargs);
  fflush(stdout);
}

void log_to_file(struct configuration_params *config, const char *msg, ...) {}

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