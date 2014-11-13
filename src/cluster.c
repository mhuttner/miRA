#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "cluster.h"
#include "parse_sam.h"
#include "errors.h"
#include "string.h"
#include "uthash.h"

int cluster(int argc, char **argv) {
  int c;
  int gap_size = 10;
  int min_reads = 10;
  int window_size = 200;
  int max_length = 2000;
  int tmp;

  while ((c = getopt(argc, argv, "g:w:r:l:o:h")) != -1) {
    switch (c) {
    case 'g':
      tmp = atoi(optarg);
      if (tmp != 0) {
        gap_size = tmp;
      } else {
        print_error(E_INVALID_ARGUMENT);
        return E_INVALID_ARGUMENT;
      }
      break;
    case 'w':
      tmp = atoi(optarg);
      if (tmp != 0) {
        window_size = tmp;
      } else {
        print_error(E_INVALID_ARGUMENT);
        return E_INVALID_ARGUMENT;
      }
      break;
    case 'r':
      tmp = atoi(optarg);
      if (tmp != 0) {
        min_reads = tmp;
      } else {
        print_error(E_INVALID_ARGUMENT);
        return E_INVALID_ARGUMENT;
      }
      break;
    case 'l':
      tmp = atoi(optarg);
      if (tmp != 0) {
        max_length = tmp;
      } else {
        print_error(E_INVALID_ARGUMENT);
        return E_INVALID_ARGUMENT;
      }
      break;
    case 'o':
      break;
    case 'h':
      print_help();
      return E_SUCCESS;
    default:
      break;
    }
  }
  if (optind >= argc) { /* missing input file */
    printf("No SAM File specified\n\n");
    print_help();
    return E_NO_FILE_SPECIFIED;
  }
  struct sam_file *sam = NULL;
  int err = parse_sam(&sam, argv[optind]);
  if (err != E_SUCCESS) {
    print_error(err);
    return E_UNKNOWN;
  }

  struct cluster_list *list = NULL;
  struct chrom_info *chromosome_table = NULL;

  create_chromosome_table(&chromosome_table, sam);

  err = create_clusters(&list, sam);
  if (err != E_SUCCESS) {
    print_error(err);
    free_sam(sam);
    return err;
  }
  free_sam(sam);

  sort_clusters(list, compare_strand_chrom_start);
  err = merge_clusters(list, 0);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  err = filter_clusters(list, min_reads);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }

  err = merge_clusters(list, gap_size);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  err = extend_clusters(list, &chromosome_table, window_size);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  err = merge_extended_clusters(list, max_length);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  err = filter_extended_clusters(list, max_length);
  if (err != E_SUCCESS) {
    goto error_clusters;
  }
  sort_clusters(list, compare_chrom_flank);

  print_bed_file(NULL, list);

  free_chromosome_table(&chromosome_table);
  free_clusters(list);

  return E_SUCCESS;

error_clusters:
  print_error(err);
  free_clusters(list);
  free_chromosome_table(&chromosome_table);
  return err;
}

int print_help() {
  printf("Description:\n"
         "    cluster generates a list of main expression contigs based on\n"
         "    alignment data.\n"
         "Usage: miRA cluster [-g gapsize] [-w windowsize] [-f minreadnumber]\n"
         "                    [-l maxlength] [-o outputfile] [-h] "
         "<input SAM file> \n");
  return E_SUCCESS;
}

int create_clusters(struct cluster_list **list, struct sam_file *sam) {
  struct cluster_list *tmp_list =
      (struct cluster_list *)malloc(sizeof(struct cluster_list));
  if (tmp_list == NULL) {
    return E_MALLOC_FAIL;
  }
  tmp_list->capacity = sam->n;
  tmp_list->clusters =
      (struct cluster *)malloc(tmp_list->capacity * sizeof(struct cluster));
  if (tmp_list->clusters == NULL) {
    free(tmp_list);
    return E_MALLOC_FAIL;
  }
  size_t n = 0;
  for (size_t i = 0; i < sam->n; i++) {
    int result = sam_to_cluster(tmp_list->clusters + n, sam->entries + i, n);
    if (result != E_SUCCESS)
      continue;
    n++;
  }
  tmp_list->n = n;
  *list = tmp_list;
  return E_SUCCESS;
}

int sort_clusters(struct cluster_list *list,
                  int (*comparison_func)(const void *c1, const void *c2)) {
  qsort(list->clusters, list->n, sizeof(struct cluster), comparison_func);
  return E_SUCCESS;
}

int compare_strand_chrom_start(const void *c1, const void *c2) {
  struct cluster *cl1 = (struct cluster *)c1;
  struct cluster *cl2 = (struct cluster *)c2;
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
  struct cluster *cl1 = (struct cluster *)c1;
  struct cluster *cl2 = (struct cluster *)c2;
  int diff;
  diff = strcmp(cl1->chrom, cl2->chrom);
  if (diff != 0)
    return diff;
  diff = cl1->flank_start - cl2->flank_start;
  return diff;
}
int compare_strand_chrom_flank(const void *c1, const void *c2) {
  struct cluster *cl1 = (struct cluster *)c1;
  struct cluster *cl2 = (struct cluster *)c2;
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
  struct cluster *new_clusters =
      (struct cluster *)malloc(list->n * sizeof(struct cluster));
  if (new_clusters == NULL) {
    return E_MALLOC_FAIL;
  }
  new_clusters[0] = list->clusters[0];
  struct cluster *top = new_clusters;
  for (size_t i = 1; i < list->n; i++) {
    struct cluster *c = list->clusters + i;
    if (c->strand != top->strand || strcmp(c->chrom, top->chrom) != 0 ||
        c->start > top->end + max_gap) {
      top++;
      *top = *c;
    } else {
      if (c->end > top->end)
        top->end = c->end;
      top->readcount += c->readcount;
      free_cluster(c);
    }
  }
  size_t new_size = top + 1 - new_clusters;
  struct cluster *tmp = (struct cluster *)realloc(
      new_clusters, new_size * sizeof(struct cluster));
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
  struct cluster *new_clusters =
      (struct cluster *)malloc(list->n * sizeof(struct cluster));
  if (new_clusters == NULL) {
    return E_MALLOC_FAIL;
  }
  struct cluster *top = new_clusters;
  for (size_t i = 0; i < list->n; i++) {
    struct cluster *c = list->clusters + i;
    if (c->readcount >= minreads) {
      *top = *c;
      top++;
    } else {
      free_cluster(c);
    }
  }
  size_t new_size = top - new_clusters;
  struct cluster *tmp = (struct cluster *)realloc(
      new_clusters, new_size * sizeof(struct cluster));
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
    c = list->clusters + i;
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
  struct cluster *new_clusters =
      (struct cluster *)malloc(list->n * sizeof(struct cluster));
  if (new_clusters == NULL) {
    return E_MALLOC_FAIL;
  }
  new_clusters[0] = list->clusters[0];
  struct cluster *top = new_clusters;
  int total_readcount = 0;
  for (size_t i = 1; i < list->n; i++) {
    struct cluster *c = list->clusters + i;
    if (c->strand != top->strand || strcmp(c->chrom, top->chrom) != 0 ||
        c->flank_start > top->flank_end ||
        c->flank_end - top->flank_start > max_length) {
      top->readcount = total_readcount;
      top++;
      *top = *c;
      total_readcount = top->readcount;
    } else {
      total_readcount += c->readcount;
      if (c->readcount >= top->readcount) {
        top->start = c->start;
        top->end = c->end;
        top->readcount = c->readcount;
      }
      if (c->flank_end > top->flank_end) {
        top->flank_end = c->flank_end;
      }
      top->readcount += c->readcount;
      free_cluster(c);
    }
  }
  size_t new_size = top + 1 - new_clusters;
  struct cluster *tmp = (struct cluster *)realloc(
      new_clusters, new_size * sizeof(struct cluster));
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
  struct cluster *new_clusters =
      (struct cluster *)malloc(list->n * sizeof(struct cluster));
  if (new_clusters == NULL) {
    return E_MALLOC_FAIL;
  }
  struct cluster *top = new_clusters;
  for (size_t i = 0; i < list->n; i++) {
    struct cluster *c = list->clusters + i;
    if (c->flank_end - c->flank_start < max_length) {
      *top = *c;
      top++;
    } else {
      free_cluster(c);
    }
  }
  size_t new_size = top - new_clusters;
  struct cluster *tmp = (struct cluster *)realloc(
      new_clusters, new_size * sizeof(struct cluster));
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

int print_bed_file(char *filename, struct cluster_list *list) {
  struct cluster *c = NULL;
  long fixed_start;
  long fixed_flank_start;
  for (size_t i = 0; i < list->n; i++) {
    c = list->clusters + i;
    /* BED format is 0 based and does not include the last base */
    fixed_start = (c->start > 0) ? c->start - 1 : 0;
    fixed_flank_start = (c->flank_start > 0) ? c->flank_start - 1 : 0;
    fprintf(stdout, "%s\t%ld\t%ld\tCluster_%ld\t%d\t%c\t%ld\t%ld\t%d\t%ld\n",
            c->chrom, fixed_flank_start, c->flank_end, i, 0, c->strand,
            fixed_start, c->end, 0, c->readcount);
  }
  return E_SUCCESS;
}

int create_chromosome_table(struct chrom_info **table, struct sam_file *sam) {
  struct chrom_info *info = NULL;
  struct sq_header *h = NULL;

  for (size_t i = 0; i < sam->header_n; i++) {
    info = (struct chrom_info *)malloc(sizeof(struct chrom_info));
    h = sam->headers + i;
    strncpy(info->name, h->sn, 1024);
    info->length = h->ln;
    HASH_ADD_STR(*table, name, info);
  }

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
    free_cluster(list->clusters + i);
  }
  free(list->clusters);
  free(list);
  return E_SUCCESS;
}
int free_cluster(struct cluster *c) {
  free(c->chrom);
  return E_SUCCESS;
}
