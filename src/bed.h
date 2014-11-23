#include "cluster.h"

#ifndef BED_FILE_IO_H
#define BED_FILE_IO_H

int write_bed_file(char *filename, struct cluster_list *list);
int read_bed_file(struct cluster_list **list, char *filename);
int parse_bed_line(struct cluster **result, char *line);

#endif