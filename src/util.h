

#include <stdlib.h>

#ifndef UTIL_H
#define UTIL_H

enum log_level {
  LOG_LEVEL_QUIET = 0,
  LOG_LEVEL_BASIC = 1,
  LOG_LEVEL_VERBOSE = 2,
};

struct configuration_params {
  int log_level;
  int openmp_thread_count;
  int cluster_gap_size;
  int cluster_min_reads;
  int cluster_flank_size;
  int cluster_max_length;

  int max_precursor_length;
  int min_precursor_length;
  double max_mfe_per_nt;
  int max_hairpin_count;
  int min_double_strand_length;
  int permutation_count;
  double max_pvalue;
  double min_coverage;
  double min_paired_fraction;
  int min_duplex_length;
  int max_duplex_length;
  int allow_three_mismatches;
  int allow_two_terminal_mismatches;

  int create_coverage_plots;
  int create_structure_plots;
  int create_structure_coverage_plots;
  int cleanup_auxiliary_files;
};

struct text_buffer {
  char *data;
  char *start;
  size_t capacity;
  size_t free_space;
};

int create_text_buffer(struct text_buffer **buffer);
void print_to_text_buffer(struct text_buffer *buffer, const char *msg, ...);
int free_text_buffer(struct text_buffer *buffer);

int initialize_configuration(struct configuration_params **config,
                             char *config_file);

int reverse_complement_sequence_string(char **result, char *seq, size_t n);
int create_file_path(char **file_path, const char *path, const char *filename);

void log_configuration(struct configuration_params *config);
void log_basic(int loglevel, const char *msg, ...);
void log_basic_timestamp(int loglevel, const char *msg, ...);
void log_verbose(int loglevel, const char *msg, ...);

double mean(double *list, int n);
double sd(double *list, int n, double mean);
double pvalue(double mean, double sd, double value);

#endif