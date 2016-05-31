#include <stddef.h>
#include "parse_sam.h"
#include "uthash.h"
#include "defs.h"
#include "util.h"

#ifndef CLUSTER_H
#define CLUSTER_H

struct cluster {
  u64 id;
  char strand;
  char *chrom;
  u64 start;
  u64 end;
  u64 readcount;
  u64 flank_start;
  u64 flank_end;
};

struct cluster_list {
  size_t n;
  size_t capacity;
  struct cluster **clusters;
};

struct chrom_info {
  char name[1024];
  long length;
  UT_hash_handle hh;
};

int cluster(int argc, char **argv);
int cluster_main(struct configuration_params *config, char *sam_file,
                 char *output_file, char *selected_crom);

int parse_clusters(struct configuration_params *config,
                   struct chrom_info **table, struct cluster_list **list,
                   char *file, char *selected_crom);
int create_clusters(struct cluster_list **list, size_t n);

int sort_clusters(struct cluster_list *list,
                  int (*comparison_func)(const void *c1, const void *c2));
int compare_strand_chrom_start(const void *c1, const void *c2);
int compare_chrom_flank(const void *c1, const void *c2);
int compare_strand_chrom_flank(const void *c1, const void *c2);
int merge_clusters(struct cluster_list *list, int max_gap);
int filter_clusters(struct cluster_list *list, int minreads);
int extend_clusters(struct cluster_list *list, struct chrom_info **table,
                    int window);
int merge_extended_clusters(struct cluster_list *list, int max_length);
int filter_extended_clusters(struct cluster_list *list, int max_length);
int sam_to_cluster(struct cluster *cluster, struct sam_entry *entry, long id);

int free_chromosome_table(struct chrom_info **table);

int free_clusters(struct cluster_list *list);
int free_cluster(struct cluster *c);

#endif