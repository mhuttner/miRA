#include <stddef.h>
#include "parse_sam.h"
#include "uthash.h"

#ifndef CLUSTER_H
#define CLUSTER_H

struct cluster {
  long id;
  char strand;
  char *chrom;
  long start;
  long end;
  long readcount;
  long flank_start;
  long flank_end;
};

struct cluster_list {
  size_t n;
  size_t capacity;
  struct cluster *clusters;
};

struct chrom_info {
  char name[1024];
  long length;
  UT_hash_handle hh;
};

int cluster(int argc, char **argv);
int print_help();

int create_clusters(struct cluster_list **list, struct sam_file *sam);

int sort_clusters(struct cluster_list *list, int compare_strand);
int compare_clusters(const void *c1, const void *c2);
int compare_extended_clusters(const void *c1, const void *c2);
int compare_extended_clusters_strand(const void *c1, const void *c2);
int merge_clusters(struct cluster_list *list, int max_gap);
int filter_clusters(struct cluster_list *list, int minreads);
int extend_clusters(struct cluster_list *list, struct chrom_info **table,
                    int window);
int merge_extended_clusters(struct cluster_list *list, int max_length);
int filter_extended_clusters(struct cluster_list *list, int max_length);
int sam_to_cluster(struct cluster *cluster, struct sam_entry *entry, long id);

int print_bed_file(char *filename, struct cluster_list *list);

int create_chromosome_table(struct chrom_info **table, struct sam_file *sam);
int free_chromosome_table(struct chrom_info **table);

int free_clusters(struct cluster_list *list);
int free_cluster(struct cluster *c);

#endif