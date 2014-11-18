
#include <stdio.h>
#include "errors.h"
#include "cluster.h"
#include "bed_file_io.h"
#include "defs.h"

int write_bed_file(char *filename, struct cluster_list *list) {
  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    return E_FILE_WRITING_FAILED;
  }
  struct cluster *c = NULL;
  u64 fixed_start;
  u64 fixed_flank_start;
  for (size_t i = 0; i < list->n; i++) {
    c = list->clusters[i];
    /* BED format is 0 based and does not include the last base */
    fixed_start = (c->start > 0) ? c->start - 1 : 0;
    fixed_flank_start = (c->flank_start > 0) ? c->flank_start - 1 : 0;
    fprintf(fp, "%s\t%llu\t%llu\tCluster_%ld\t%d\t%c\t%llu\t%llu\t%d\t%llu\n",
            c->chrom, fixed_flank_start, c->flank_end, i, 0, c->strand,
            fixed_start, c->end, 0, c->readcount);
  }
  fclose(fp);
  return E_SUCCESS;
}

int read_bed_file(struct cluster_list **list, char *filename) {
  static const int MAXLINELENGHT = 2048;
  static const int STARTINGSIZE = 1024;
  struct cluster_list *tmp_list = NULL;
  int err = create_clusters(&tmp_list, STARTINGSIZE);
  if (err) {
    return err;
  }
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    free(tmp_list->clusters);
    free(tmp_list);
    return E_FILE_NOT_FOUND;
  }
  char line[MAXLINELENGHT];

  struct cluster *c = NULL;
  while (fgets(line, sizeof(line), fp) != NULL) {
    err = parse_bed_line(&c, line);
  }
  fclose(fp);

  return E_SUCCESS;
}

int parse_bed_line(struct cluster **result, char *line) {
  const char seperator = '\t';
  const int num_entries = 10;

  struct cluster *c = (struct cluster *)malloc(sizeof(struct cluster));
  if (c == NULL) {
    return E_MALLOC_FAIL;
  }
  char **tokens = (char **)malloc(num_entries * sizeof(char *));
  if (tokens == NULL) {
    free(c);
    return E_MALLOC_FAIL;
  }
  char *start = line;
  char *end = NULL;
  int line_done = 0;
  int current_token = 0;

  while (!line_done) {
    if (current_token >= num_entries) {
      break;
    }
    end = strchr(start, seperator);
    if (end == NULL) {
      end = strchr(start, '\n');
      if (end == NULL)
        end = strchr(start, '\0');
      line_done = 1;
    }
    long l = end - start;
    if (l == 0) {
      start += 1;
      continue;
    }
    tokens[current_token] = (char *)malloc((l + 1) * sizeof(char));
    memcpy(tokens[current_token], start, l);
    tokens[current_token][l] = 0;
    current_token++;
    start = end + 1;
  }

  if (current_token < num_entries) {
    for (int i = 0; i < current_token; i++) {
      free(tokens[i]);
    }
    free(tokens);
    return E_INVALID_BED_LINE;
  }
  /* discard unneeded tokens */
  free(tokens[4]);
  free(tokens[8]);

  char *check = NULL;

  return E_SUCCESS;
}
/*



















*/
