#include <stdlib.h>
#include "../src/parse_sam.h"
#include "../src/cluster.h"
#include "../src/errors.h"
#include "testerino.h"

int create_test_clusters(struct cluster_list **list, int n, char *strands,
                         char *chromosomes[], long *starts, long *ends,
                         long *reads, long *flank_starts, long *flank_ends) {
  struct cluster_list *tmp =
      (struct cluster_list *)malloc(sizeof(struct cluster_list));
  tmp->n = n;
  tmp->capacity = n;
  struct cluster **clusters =
      (struct cluster **)malloc(n * sizeof(struct cluster *));
  tmp->clusters = clusters;

  char default_strand = '+';
  const char *default_chromosome = "chrom_1";
  long default_start = 100;
  long default_end = 120;
  long default_reads = 5;
  long default_flank_start = 50;
  long default_flank_end = 200;
  struct cluster *c = NULL;
  for (size_t i = 0; i < tmp->n; i++) {
    c = (struct cluster *)malloc(sizeof(struct cluster));
    c->id = i;
    c->strand = default_strand;
    c->chrom = (char *)malloc(21 * sizeof(char));
    strncpy(c->chrom, default_chromosome, 20);
    c->start = default_start;
    c->end = default_end;
    c->readcount = default_reads;
    c->flank_start = default_flank_start;
    c->flank_end = default_flank_end;
    if (strands != NULL) {
      c->strand = strands[i];
    }
    if (chromosomes != NULL) {
      strncpy(c->chrom, chromosomes[i], 20);
    }
    if (starts != NULL) {
      c->start = starts[i];
    }
    if (ends != NULL) {
      c->end = ends[i];
    }
    if (reads != NULL) {
      c->readcount = reads[i];
    }
    if (flank_starts != NULL) {
      c->flank_start = flank_starts[i];
    }
    if (flank_ends != NULL) {
      c->flank_end = flank_ends[i];
    }
    tmp->clusters[i] = c;
  }
  *list = tmp;
  return E_SUCCESS;
}

void test_sort_clusters(struct test *t) {
  t_set_msg(t, "Testing sorting of clusters...");
  struct cluster_list *list = NULL;
  char strands[4] = {'-', '+', '+', '+'};
  long starts[4] = {10, 30, 15, 10};
  long ends[4] = {20, 40, 25, 20};
  create_test_clusters(&list, 4, strands, NULL, starts, ends, NULL, NULL, NULL);
  sort_clusters(list, compare_strand_chrom_start);
  t_assert_msg(t, list->clusters[0]->id == 3, "Wrong cluster order");
  t_assert_msg(t, list->clusters[1]->id == 2, "Wrong cluster order");
  t_assert_msg(t, list->clusters[2]->id == 1, "Wrong cluster order");
  t_assert_msg(t, list->clusters[3]->id == 0, "Wrong cluster order");
  free_clusters(list);
}

void test_merge_clusters(struct test *t) {
  t_set_msg(t, "Testing merging of clusters...");
  struct cluster_list *list = NULL;
  long starts[4] = {10, 100, 30, 200};
  long ends[4] = {24, 120, 44, 220};
  create_test_clusters(&list, 4, NULL, NULL, starts, ends, NULL, NULL, NULL);
  sort_clusters(list, compare_strand_chrom_start);
  int err = merge_clusters(list, 10);
  t_assert_msg(t, err == E_SUCCESS, "merging failed");
  struct cluster *c;
  for (size_t i = 0; i < list->n; i++) {
    c = list->clusters[i];
    t_log(t, "%ld %ld %ld \n", c->id, c->start, c->end);
  }
  t_assert_msg(t, list->n == 3, "Merging failed");

  free_clusters(list);
}

void test_filter_clusters(struct test *t) {
  t_set_msg(t, "Testing filtering of clusters...");
  struct cluster_list *list = NULL;

  long test_readcounts[4] = {5, 10, 20, 3};
  int err = create_test_clusters(&list, 4, NULL, NULL, NULL, NULL,
                                 test_readcounts, NULL, NULL);

  err = filter_clusters(list, 10);
  t_assert_msg(t, err == E_SUCCESS, "Filtering failed");
  t_assert_msg(t, list->n == 2, "Filtering incorrect");
  struct cluster *c;
  for (size_t i = 0; i < list->n; i++) {
    c = list->clusters[i];
    t_assert_msg(t, c->readcount >= 10,
                 "Cluster with not enough reads in result");
  }
  free_clusters(list);
}

void test_merge_extended_clusters(struct test *t) {
  t_set_msg(t, "Testing merging of extended clusters...");
  struct cluster_list *list = NULL;
  long test_starts[] = {110, 210, 710, 410};
  long test_ends[] = {120, 220, 720, 420};
  long test_readcounts[4] = {10, 20, 30, 40};
  long test_flank_starts[4] = {100, 200, 700, 400};
  long test_flank_ends[4] = {300, 400, 900, 600};
  create_test_clusters(&list, 4, NULL, NULL, test_starts, test_ends,
                       test_readcounts, test_flank_starts, test_flank_ends);
  sort_clusters(list, compare_strand_chrom_start);
  int err = merge_extended_clusters(list, 1000);
  t_assert_msg(t, err == E_SUCCESS, "Merging failed");
  struct cluster *c = NULL;
  for (size_t i = 0; i < list->n; i++) {
    c = list->clusters[i];
    t_log(t, "%ld %c %s %ld %ld %ld %ld %ld \n", c->id, c->strand, c->chrom,
          c->start, c->end, c->readcount, c->flank_start, c->flank_end);
  }
  t_assert_msg(t, list->clusters[0]->start == 410, "Merging incorrect");
  t_assert_msg(t, list->clusters[0]->flank_start == 100, "Merging incorrect");
  t_assert_msg(t, list->clusters[0]->flank_end == 600, "Merging incorrect");

  free_clusters(list);
}
