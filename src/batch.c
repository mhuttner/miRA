#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "cluster.h"
#include "vfold.h"
#include "coverage.h"
#include "util.h"
#include "errors.h"
#include "batch.h"
#include "reporting.h"

static int print_help();

int batch(int argc, char **argv) {
  char *config_file = NULL;
  int c;
  int log_level = LOG_LEVEL_BASIC;

  while ((c = getopt(argc, argv, "c:s:hvq")) != -1) {
    switch (c) {
    case 's':
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
  if (optind + 3 > argc) { /* missing input file(s) */
    printf("Not enough Input Files specified\n\n");
    print_help();
    return E_NO_FILE_SPECIFIED;
  }
  struct configuration_params *config = NULL;
  initialize_configuration(&config, config_file);
  if (log_level != LOG_LEVEL_BASIC) {
    config->log_level = log_level;
  }
  log_configuration(config);

  int err;
  char *sam_file = argv[optind];
  char *fasta_file = argv[optind + 1];
  char *output_path = argv[optind + 2];
  char *mira_bin = argv[-1];

  err = create_directory_if_ne(output_path);
  if (err) {
    return err;
  }

  struct sam_file *sam_headers = NULL;
  err = parse_sam_headers(&sam_headers, sam_file);
  if (err) {
    return err;
  }
  for (size_t i = 0; i < sam_headers->header_n; ++i) {
    char *selected_chrom = sam_headers->headers[i]->sn;
    log_basic(
        config->log_level,
        "##############################################################\n");
    log_basic(
        config->log_level,
        "#                        BATCH %3ld                           #\n",
        i + 1);
    log_basic(
        config->log_level,
        "##############################################################\n");
    log_basic(config->log_level, "Chromosome: %s \n\n", selected_chrom);

    char *full_output_path = NULL;
    create_file_path(&full_output_path, output_path, selected_chrom);
    batch_main(config, mira_bin, sam_file, fasta_file, full_output_path,
               selected_chrom);
    free(full_output_path);
  }

  free_sam(sam_headers);
  free(config);
  return E_SUCCESS;
}

int batch_main(struct configuration_params *config, char *mira_bin,
               char *sam_file, char *fasta_file, char *output_path,
               char *selected_chrom) {
  int err;

  err = create_directory_if_ne(output_path);
  if (err) {
    return err;
  }

  char bed_filename[] = "cluster_contigs.bed";
  char *bed_file_path = NULL;
  create_file_path(&bed_file_path, output_path, bed_filename);

  err = cluster_main(config, sam_file, bed_file_path, selected_chrom);
  if (err) {
    free(bed_file_path);
    return err;
  }
  char mira_filename[] = "fold_candidates.miRA";
  char *mira_file_path = NULL;
  create_file_path(&mira_file_path, output_path, mira_filename);
  err = vfold_main(config, bed_file_path, fasta_file, mira_file_path,
                   selected_chrom);
  if (err) {
    free(bed_file_path);
    free(mira_file_path);
    return err;
  }
  err = coverage_main(config, mira_bin, mira_file_path, sam_file, output_path,
                      selected_chrom);

  free(bed_file_path);
  free(mira_file_path);
  if (err) {
    return err;
  }
  log_basic(config->log_level, "All steps completed successfully.\n");
  return E_SUCCESS;
}

static int print_help() {
  printf("Description:\n"
         "    Coverage based verification on micro RNA candidates \n"
         "    Batches all files based on the chromosome (rname) and \n"
         "    runs a full miRA analysis for each chromosome separately\n"
         "Usage: miRA batch <input SAM file> <input FASTA file> <output "
         "directory>\n");
  return E_SUCCESS;
}