#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "cluster.h"
#include "vfold.h"
#include "coverage.h"
#include "util.h"
#include "errors.h"
#include "full.h"
#include "reporting.h"

static int print_help();

int full(int argc, char **argv) {
  char *config_file = NULL;
  int c;
  int log_level = LOG_LEVEL_BASIC;

  while ((c = getopt(argc, argv, "c:hvq")) != -1) {
    switch (c) {
    case 'c':
      config_file = optarg;
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
  if (optind + 3 > argc) { /* missing input file(s) */
    printf("Not enough Input Files specified\n\n");
    print_help();
    return E_NO_FILE_SPECIFIED;
  }
  struct configuration_params *config = NULL;
  initialize_configuration(&config, config_file);
  log_configuration(config);

  int err;
  char *sam_file = argv[optind];
  char *fasta_file = argv[optind + 1];
  char *output_path = argv[optind + 2];

  err = create_directory_if_ne(output_path);
  if (err) {
    return err;
  }

  char bed_filename[] = "cluster_contigs.bed";
  char *bed_file_path = NULL;
  create_file_path(&bed_file_path, output_path, bed_filename);
  err = cluster_main(config, sam_file, bed_file_path);
  if (err) {
    free(bed_file_path);
    return err;
  }
  char mira_filename[] = "fold_candidates.miRA";
  char *mira_file_path = NULL;
  create_file_path(&mira_file_path, output_path, mira_filename);
  err = vfold_main(config, bed_file_path, fasta_file, mira_file_path);
  if (err) {
    free(bed_file_path);
    free(mira_file_path);
    return err;
  }

  err = coverage_main(config, mira_file_path, sam_file, output_path);
  free(bed_file_path);
  free(mira_file_path);
  free(config);
  if (err) {
    return err;
  }
  log_basic(config->log_level,
            "All steps completed successfully. Exiting... \n");
  return E_SUCCESS;
}

static int print_help() {
  printf("Description:\n"
         "    Coverage based verification on micro RNA candidates \n"
         "Usage: miRA full <input SAM file> <input FASTA file> <output "
         "directory>\n");
  return E_SUCCESS;
}