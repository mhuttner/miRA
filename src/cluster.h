#include <stddef.h>

#ifndef CLUSTER_H
#define CLUSTER_H 

struct cluster {
    long id;
    char strand;
    char* chrom;
    long start;
    long end;
    long readcount;
};

struct cluster_list {
    size_t n;
    size_t capacity;
    struct cluster* clusters;
};


int cluster(int argc,char** argv);
int print_help();

int create_clusters(struct cluster_list** index,struct sam_file* sam);
int sam_to_cluster(struct cluster* cluster,struct sam_entry* entry,long id);

int free_clusters(struct cluster_list* list);
int free_cluster(struct cluster c);

#endif