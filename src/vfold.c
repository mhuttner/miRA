#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>
#include "vfold.h"
#include "errors.h"
#include "fasta.h"
#include "cluster.h"
#include "bed.h"
#include "Lfold/lFold.h"

static int print_help();

int vfold(int argc, char *argv[]) {
  int c;
  char default_output_file[] = "contigs.bed";
  char *output_file = default_output_file;

  while ((c = getopt(argc, argv, "o:h")) != -1) {
    switch (c) {
    case 'o':
      output_file = optarg;
      break;
    case 'h':
      print_help();
      return E_SUCCESS;
    default:
      break;
    }
  }
  if (optind + 2 > argc) { /* missing input file */
    printf("No Input Files specified\n\n");
    print_help();
    return E_NO_FILE_SPECIFIED;
  }
  struct cluster_list *c_list = NULL;
  struct genome_sequence *seq_table = NULL;
  struct sequence_list *seq_list = NULL;
  int err;

  err = read_bed_file(&c_list, argv[optind]);
  if (err) {
    goto bed_read_err;
  }

  err = read_fasta_file(&seq_table, argv[optind + 1]);
  if (err) {
    goto fasta_read_err;
  }
  err = map_clusters(&seq_list, c_list, seq_table);
  if (err) {
    goto map_error;
  }
  /*clusters freed by map_clusters */
  free_sequence_table(seq_table);

  fold_sequences(seq_list);

  return E_SUCCESS;
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
         "Usage: miRA fold <input BED file> <input FASTA file>\n");
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
int fold_sequences(struct sequence_list *seq_list) {
  struct foldable_sequence *fs = NULL;

  for (size_t i = 92; i < seq_list->n; i++) {
    fs = seq_list->sequences[i];
    Lfold(fs->seq, NULL, fs->n);
  }

  return E_SUCCESS;
};

int reverse_complement(struct foldable_sequence *s) {
  const char *pairs[] = {"AT", "GC", "UA", "YR", "SS",
                         "WW", "KM", "BV", "DH", "NN"};
  const int pairs_n = 10;
  char *tmp = (char *)malloc(s->n * sizeof(char));
  if (tmp == NULL) {
    return E_MALLOC_FAIL;
  }
  size_t k = s->n - 1;
  tmp[k--] = 0;
  for (size_t i = 0; i < s->n; i++) {
    for (int j = 0; j < pairs_n; j++) {
      if (s->seq[i] == pairs[j][0]) {
        tmp[k--] = pairs[j][1];
        break;
      }
      if (s->seq[i] == pairs[j][1]) {
        tmp[k--] = pairs[j][0];
        break;
      }
    }
  }
  free(s->seq);
  s->seq = tmp;

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
  free(s);
  return E_SUCCESS;
}
