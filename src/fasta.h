

#ifndef FASTA_H
#define FASTA_H

#include <stddef.h>
#include "uthash.h"

struct genome_sequence {
  char chrom[1024];
  char *data;
  size_t n;
  size_t capacity;
  UT_hash_handle hh;
};

int read_fasta_file(struct genome_sequence **sequence_table, char *filename,
                    char *selected_crom);
int create_genome_sequence(struct genome_sequence **seq);
int strip_newlines(char **dest, size_t *nsize, char *buffer, size_t size);
int append_to_genome_sequence(struct genome_sequence *seq, char *data,
                              size_t size);
int free_sequence_table(struct genome_sequence *table);
#endif