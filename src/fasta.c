#include "fasta.h"
#include "errors.h"
#include "fasta.h"
#include <stdlib.h>
#include <stdio.h>

int read_fasta_file(struct genome_sequence **sequence_table, char *filename) {
  static const char name_marker = '>';
  static const int MAXLINELENGHT = 2048;
  struct genome_sequence *current_seq = NULL;
  int err, l;
  char *sep;
  char line[MAXLINELENGHT];
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    return E_FILE_NOT_FOUND;
  }
  while (fgets(line, sizeof(line), fp) != NULL) {
    if (line[0] == name_marker) {
      if (current_seq != NULL) {
        HASH_ADD_STR(*sequence_table, chrom, current_seq);
      }
      create_genome_sequence(&current_seq);
      sep = strchr(line, ' ');
      if (sep == NULL) {
        sep = strchr(line, '\n');
      }
      if (sep == NULL) {
        sep = strchr(line, '\0');
      }
      if (sep == NULL) {
        err = E_UNKNOWN_FILE_IO_ERROR;
        goto parse_error;
      }
      l = sep - line - 1;
      memcpy(current_seq->chrom, line + 1, l);
      current_seq->chrom[l] = 0;
      continue;
    }
    sep = strchr(line, '\n');
    if (sep == NULL) {
      sep = strchr(line, '\0');
    }
    if (sep == NULL) {
      err = E_UNKNOWN_FILE_IO_ERROR;
      goto parse_error;
    }
    err = append_to_genome_sequence(current_seq, line, sep - line);
    if (err) {
      goto parse_error;
    }
  }
  if (*sequence_table == NULL) {
    fclose(fp);
    return E_UNKNOWN_FILE_IO_ERROR;
  }
  HASH_ADD_STR(*sequence_table, chrom, current_seq);
  fclose(fp);
  return E_SUCCESS;

parse_error:
  free_sequence_table(*sequence_table);
  fclose(fp);
  return err;
}

int create_genome_sequence(struct genome_sequence **seq) {
  static const int STARTINGSIZE = 4096;
  struct genome_sequence *s =
      (struct genome_sequence *)malloc(sizeof(struct genome_sequence));
  if (s == NULL) {
    return E_MALLOC_FAIL;
  }
  s->data = (char *)malloc(STARTINGSIZE * sizeof(char));
  if (s->data == NULL) {
    free(s);
    return E_MALLOC_FAIL;
  }
  s->capacity = STARTINGSIZE;
  s->n = 0;
  *seq = s;
  return E_SUCCESS;
}
int strip_newlines(char **dest, size_t *nsize, char *buffer, size_t size) {
  char *new_buffer = (char *)malloc(size * sizeof(char));
  size_t new_size = 0;
  for (size_t i = 0; i < size; i++) {
    if (buffer[i] != '\n') {
      new_buffer[new_size++] = buffer[i];
    }
  }
  char *tmp = (char *)realloc(new_buffer, new_size * sizeof(char));
  if (tmp == NULL) {
    free(new_buffer);
    return E_REALLOC_FAIL;
  }
  *dest = tmp;
  *nsize = new_size;
  return E_SUCCESS;
}

int append_to_genome_sequence(struct genome_sequence *seq, char *data,
                              size_t size) {
  if (seq == NULL)
    return E_SUCCESS;
  while (seq->capacity - seq->n < size) {
    seq->capacity *= 2;
    char *tmp = (char *)realloc(seq->data, seq->capacity * sizeof(char));
    if (tmp == NULL) {
      return E_REALLOC_FAIL;
    }
    seq->data = tmp;
  }
  memcpy(seq->data + seq->n, data, size);
  seq->n += size;
  return E_SUCCESS;
}

int free_sequence_table(struct genome_sequence *table) {
  struct genome_sequence *s = NULL;
  struct genome_sequence *tmp = NULL;
  HASH_ITER(hh, table, s, tmp) {
    HASH_DEL(table, s);
    free(s->data);
    free(s);
  }
  return E_SUCCESS;
}