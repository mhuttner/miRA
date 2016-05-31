#ifndef BATCH_H
#define BATCH_H

int batch(int argc, char **argv);
int batch_main(struct configuration_params *config, char *mira_bin,
               char *sam_file, char *fasta_file, char *output_path,
               char *selected_chrom);

#endif