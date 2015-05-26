#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "cluster.h"
#include "parse_sam.h"
#include "bed.h"
#include "errors.h"
#include "string.h"
#include "uthash.h"
#include "util.h"

static int print_help();

int cluster(int argc, char **argv) {
  int c;
  char *config_file = NULL;
  int log_level = LOG_LEVEL_BASIC;
  char default_output_file[] = "contigs.bed";
  char *output_file = default_output_file;

  while ((c = getopt(argc, argv, "c:o:h")) != -1) {
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
      break;
    case 'q':
      log_level = LOG_LEVEL_QUIET;
      break;
    default:
      break;
    }
  }
  if (optind >= argc) { /* missing input file */
    printf("No SAM File specified\n\n");
    print_help();
    return E_NO_FILE_SPECIFIED;
  }
  struct configuration_params *config = NULL;
  initialize_configuration(&config, config_file);
  config->log_level = log_level;
  log_configuration(config);
  int err;
  err = cluster_main(config, argv[optind], output_file);
  free(config);
  return err;
}

static int print_help() {
  printf(
      "Description:\n"
      "    cluster generates a list of main expression contigs based on\n"
      "    alignment data.\n"
      "Usage: miRA cluster [-c config file] [-o output file] [-q] [-v] [-h]\n "
      "    <input SAM file> \n");
  return E_SUCCESS;
}

int cluster_main(struct configuration_params *config, char *sam_file,
                 char *output_file) {
  struct cluster_list *list = NULL;
  struct chrom_info *chromosome_table = NULL;
  log_basic_timestamp(config->log_level, "Clustering reads...\n");

  parse_clusters(&chromosome_table, &list, sam_file);

  int err;

  sort_clusters(list, compare_strand_chrom_start);
  err = merge_clusters(list, 0);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  err = filter_clusters(list, config->cluster_min_reads);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  err = merge_clusters(list, config->cluster_gap_size);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  err = extend_clusters(list, &chromosome_table, config->cluster_flank_size);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  err = merge_extended_clusters(list, config->cluster_max_length);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  err = filter_extended_clusters(list, config->cluster_max_length);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  sort_clusters(list, compare_chrom_flank);

  err = write_bed_file(output_file, list);
  if (err) {
    goto error_clusters;
  }

  free_chromosome_table(&chromosome_table);
  free_clusters(list);
  log_basic_timestamp(config->log_level,
                      "Clustering completed successfully.\n");

  return E_SUCCESS;

error_clusters:
  print_error(err);
  free_clusters(list);
  free_chromosome_table(&chromosome_table);
  return err;
}

int parse_clusters(struct chrom_info **table, struct cluster_list **list,
                   char *file) {
  static const int MAXLINELENGHT = 2048;
  static const int STARTINGSIZE = 1024;
  struct cluster_list *tmp_list = NULL;
  int err = create_clusters(&tmp_list, STARTINGSIZE);
  if (err != E_SUCCESS) {
    return err;
  }

  FILE *fp = fopen(file, "r");
  if (fp == NULL) {
    free(tmp_list->clusters);
    free(tmp_list);
    return E_FILE_NOT_FOUND;
  }
  char line[MAXLINELENGHT];
  struct sam_entry *tmp_entry = NULL;
  struct sq_header *tmp_header = NULL;
  struct chrom_info *info = NULL;
  struct cluster *c = NULL;
  while (fgets(line, sizeof(line), fp) != NULL) {
    int result = parse_line(&tmp_entry, line);
    if (result == E_SAM_HEADER_LINE) {
      int err = parse_header(&tmp_header, line);
      if (err != E_SUCCESS) {
        continue;
      }
      info = (struct chrom_info *)malloc(sizeof(struct chrom_info));
      if (info == NULL) {
        free_sam_header(tmp_header);
        continue;
      }
      strncpy(info->name, tmp_header->sn, 1024);
      info->length = tmp_header->ln;
      HASH_ADD_STR(*table, name, info);
      free_sam_header(tmp_header);
      continue;
    }
    if (result != E_SUCCESS) {
      continue;
    }

    if (tmp_list->n == tmp_list->capacity) {
      tmp_list->capacity *= 2;
      struct cluster **tmp = (struct cluster **)realloc(
          tmp_list->clusters, tmp_list->capacity * sizeof(struct cluster *));
      if (tmp == NULL) {
        free_sam_entry(tmp_entry);
        free(tmp_list->clusters);
        free(tmp_list);
        fclose(fp);
        return E_MALLOC_FAIL;
      }
      tmp_list->clusters = tmp;
    }
    c = (struct cluster *)malloc(sizeof(struct cluster));
    if (c == NULL) {
      continue;
    }
    sam_to_cluster(c, tmp_entry, tmp_list->n);
    free_sam_entry(tmp_entry);
    tmp_list->clusters[tmp_list->n] = c;
    tmp_list->n++;
  }
  fclose(fp);
  *list = tmp_list;

  return E_SUCCESS;
}

int create_clusters(struct cluster_list **list, size_t n) {
  struct cluster_list *tmp_list =
      (struct cluster_list *)malloc(sizeof(struct cluster_list));
  if (tmp_list == NULL) {
    return E_MALLOC_FAIL;
  }
  tmp_list->capacity = n;
  tmp_list->clusters =
      (struct cluster **)malloc(tmp_list->capacity * sizeof(struct cluster *));
  if (tmp_list->clusters == NULL) {
    free(tmp_list);
    return E_MALLOC_FAIL;
  }
  tmp_list->n = 0;
  *list = tmp_list;
  return E_SUCCESS;
}

int sort_clusters(struct cluster_list *list,
                  int (*comparison_func)(const void *c1, const void *c2)) {
  qsort(list->clusters, list->n, sizeof(struct cluster *), comparison_func);
  return E_SUCCESS;
}

int compare_strand_chrom_start(const void *c1, const void *c2) {
  struct cluster *cl1 = *(struct cluster **)c1;
  struct cluster *cl2 = *(struct cluster **)c2;
  int diff;
  diff = cl1->strand - cl2->strand;
  if (diff != 0)
    return diff;
  diff = strcmp(cl1->chrom, cl2->chrom);
  if (diff != 0)
    return diff;
  diff = cl1->start - cl2->start;
  return diff;
}
int compare_chrom_flank(const void *c1, const void *c2) {
  struct cluster *cl1 = *(struct cluster **)c1;
  struct cluster *cl2 = *(struct cluster **)c2;
  int diff;
  diff = strcmp(cl1->chrom, cl2->chrom);
  if (diff != 0)
    return diff;
  diff = cl1->flank_start - cl2->flank_start;
  return diff;
}
int compare_strand_chrom_flank(const void *c1, const void *c2) {
  struct cluster *cl1 = *(struct cluster **)c1;
  struct cluster *cl2 = *(struct cluster **)c2;
  int diff;
  diff = cl1->strand - cl2->strand;
  if (diff != 0)
    return diff;
  diff = strcmp(cl1->chrom, cl2->chrom);
  if (diff != 0)
    return diff;
  diff = cl1->flank_start - cl2->flank_start;
  return diff;
}

int merge_clusters(struct cluster_list *list, int max_gap) {
  if (list->n == 0) {
    return E_SUCCESS;
  }
  struct cluster **new_clusters =
      (struct cluster **)malloc(list->n * sizeof(struct cluster *));
  if (new_clusters == NULL) {
    return E_MALLOC_FAIL;
  }

  new_clusters[0] = list->clusters[0];
  struct cluster **top = new_clusters;
  for (size_t i = 1; i < list->n; i++) {
    struct cluster *c = list->clusters[i];
    if (c->strand != (*top)->strand || strcmp(c->chrom, (*top)->chrom) != 0 ||
        c->start > (*top)->end + max_gap) {
      top++;
      *top = c;
    } else {
      if (c->end > (*top)->end)
        (*top)->end = c->end;
      (*top)->readcount += c->readcount;
      free_cluster(c);
    }
  }
  size_t new_size = top + 1 - new_clusters;
  struct cluster **tmp = (struct cluster **)realloc(
      new_clusters, new_size * sizeof(struct cluster *));
  if (tmp == NULL) {
    free(new_clusters);
    return E_REALLOC_FAIL;
  }
  free(list->clusters);
  list->clusters = tmp;
  list->n = new_size;
  list->capacity = new_size;
  return E_SUCCESS;
}
int filter_clusters(struct cluster_list *list, int minreads) {
  struct cluster **new_clusters =
      (struct cluster **)malloc(list->n * sizeof(struct cluster *));
  if (new_clusters == NULL) {
    return E_MALLOC_FAIL;
  }
  struct cluster **top = new_clusters;
  for (size_t i = 0; i < list->n; i++) {
    struct cluster *c = list->clusters[i];
    if (c->readcount >= minreads) {
      *top = c;
      top++;
    } else {
      free_cluster(c);
    }
  }
  size_t new_size = top - new_clusters;
  struct cluster **tmp = (struct cluster **)realloc(
      new_clusters, new_size * sizeof(struct cluster *));
  if (tmp == NULL) {
    free(new_clusters);
    return E_REALLOC_FAIL;
  }
  free(list->clusters);
  list->clusters = tmp;
  list->n = new_size;
  list->capacity = new_size;
  return E_SUCCESS;
}

int extend_clusters(struct cluster_list *list, struct chrom_info **table,
                    int window) {
  struct cluster *c;
  struct chrom_info *info;
  int start;
  int end;

  for (size_t i = 0; i < list->n; i++) {
    c = list->clusters[i];
    HASH_FIND_STR(*table, c->chrom, info);
    if (info == NULL) {
      return E_CHROMOSOME_NOT_FOUND;
    }
    start = c->start - window;
    end = c->end + window;

    c->flank_start = (start > 0) ? start : 0;
    c->flank_end = (end < info->length) ? end : info->length;
  }
  return E_SUCCESS;
}
int merge_extended_clusters(struct cluster_list *list, int max_length) {
  if (list->n == 0) {
    return E_SUCCESS;
  }
  struct cluster **new_clusters =
      (struct cluster **)malloc(list->n * sizeof(struct cluster *));
  if (new_clusters == NULL) {
    return E_MALLOC_FAIL;
  }
  new_clusters[0] = list->clusters[0];
  struct cluster **top = new_clusters;
  int total_readcount = (*top)->readcount;
  for (size_t i = 1; i < list->n; i++) {
    struct cluster *c = list->clusters[i];
    if (c->strand != (*top)->strand || strcmp(c->chrom, (*top)->chrom) != 0 ||
        c->flank_start > (*top)->flank_end ||
        c->flank_end - (*top)->flank_start > max_length) {
      (*top)->readcount = total_readcount;
      top++;
      *top = c;
      total_readcount = (*top)->readcount;
    } else {
      total_readcount += c->readcount;
      if (c->readcount >= (*top)->readcount) {
        (*top)->start = c->start;
        (*top)->end = c->end;
        (*top)->readcount = c->readcount;
      }
      if (c->flank_end > (*top)->flank_end) {
        (*top)->flank_end = c->flank_end;
      }
      (*top)->readcount += c->readcount;
      free_cluster(c);
    }
  }
  size_t new_size = top + 1 - new_clusters;
  struct cluster **tmp = (struct cluster **)realloc(
      new_clusters, new_size * sizeof(struct cluster *));
  if (tmp == NULL) {
    free(new_clusters);
    return E_REALLOC_FAIL;
  }
  free(list->clusters);
  list->clusters = tmp;
  list->n = new_size;
  list->capacity = new_size;
  return E_SUCCESS;
}
int filter_extended_clusters(struct cluster_list *list, int max_length) {
  struct cluster **new_clusters =
      (struct cluster **)malloc(list->n * sizeof(struct cluster *));
  if (new_clusters == NULL) {
    return E_MALLOC_FAIL;
  }
  struct cluster **top = new_clusters;
  for (size_t i = 0; i < list->n; i++) {
    struct cluster *c = list->clusters[i];
    if (c->flank_end - c->flank_start < max_length) {
      *top = c;
      top++;
    } else {
      free_cluster(c);
    }
  }
  size_t new_size = top - new_clusters;
  struct cluster **tmp = (struct cluster **)realloc(
      new_clusters, new_size * sizeof(struct cluster *));
  if (tmp == NULL) {
    free(new_clusters);
    return E_REALLOC_FAIL;
  }
  free(list->clusters);
  list->clusters = tmp;
  list->n = new_size;
  list->capacity = new_size;
  return E_SUCCESS;
}

int sam_to_cluster(struct cluster *cluster, struct sam_entry *entry, long id) {
  const char POSITIVE_STRAND_SYMBOL = '+';
  const char NEGATIVE_STRAND_SYMBOL = '-';

  cluster->id = id;
  if (entry->flag & REV_COMPLM) {
    cluster->strand = NEGATIVE_STRAND_SYMBOL;
  } else {
    cluster->strand = POSITIVE_STRAND_SYMBOL;
  }
  int n = strlen(entry->rname) + 1;
  cluster->chrom = (char *)malloc(n * sizeof(char));
  if (cluster->chrom == NULL)
    return E_MALLOC_FAIL;
  strcpy(cluster->chrom, entry->rname);

  cluster->start = entry->pos;
  cluster->end = entry->pos + strlen(entry->seq);
  cluster->readcount = 1;

  return E_SUCCESS;
}

int free_chromosome_table(struct chrom_info **table) {
  struct chrom_info *info = NULL;
  struct chrom_info *tmp = NULL;

  HASH_ITER(hh, *table, info, tmp) {
    HASH_DEL(*table, info);
    free(info);
  }
  return E_SUCCESS;
}

int free_clusters(struct cluster_list *list) {
  for (size_t i = 0; i < list->n; i++) {
    free_cluster(list->clusters[i]);
  }
  free(list->clusters);
  free(list);
  return E_SUCCESS;
}
int free_cluster(struct cluster *c) {
  free(c->chrom);
  free(c);
  return E_SUCCESS;
}
