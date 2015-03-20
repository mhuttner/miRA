
#include <stdio.h>
#include "errors.h"
#include "cluster.h"
#include "bed.h"
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
    if (err) {
      free(tmp_list->clusters);
      free(tmp_list);
      fclose(fp);
      return err;
    }
    tmp_list->clusters[tmp_list->n] = c;
    tmp_list->n++;
    if (tmp_list->n == tmp_list->capacity) {
      tmp_list->capacity *= 2;
      struct cluster **tmp = (struct cluster **)realloc(
          tmp_list->clusters, tmp_list->capacity * sizeof(struct cluster *));
      if (tmp == NULL) {
        free(tmp_list->clusters);
        free(tmp_list);
        fclose(fp);
        return E_REALLOC_FAIL;
      }
      tmp_list->clusters = tmp;
    }
  }
  fclose(fp);
  *list = tmp_list;
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

  char *check = NULL;
  c->chrom = tokens[0];
  c->flank_start = strtol(tokens[1], &check, 10);
  if (check == tokens[1] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[1]);
  tokens[1] = NULL;

  c->flank_end = strtol(tokens[2], &check, 10);
  if (check == tokens[2] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[2]);
  tokens[2] = NULL;
  /*skip the name of the cluster and only parse the id */
  c->id = strtol(tokens[3] + 8, &check, 10);
  if (check == tokens[3] + 8 || *check != 0) {
    goto line_invalid;
  }
  free(tokens[3]);
  tokens[3] = NULL;

  free(tokens[4]);
  tokens[4] = NULL;
  c->strand = *tokens[5];
  free(tokens[5]);
  tokens[5] = NULL;
  c->start = strtol(tokens[6], &check, 10);
  if (check == tokens[6] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[6]);
  tokens[6] = NULL;
  c->end = strtol(tokens[7], &check, 10);
  if (check == tokens[7] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[7]);
  tokens[7] = NULL;
  free(tokens[8]);
  tokens[8] = NULL;
  c->readcount = strtol(tokens[9], &check, 10);
  if (check == tokens[9] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[9]);
  tokens[9] = NULL;
  *result = c;
  free(tokens);
  return E_SUCCESS;
line_invalid:
  for (int i = 0; i < num_entries; i++) {
    if (tokens[i] != NULL) {
      free(tokens[i]);
    }
  }
  free(tokens);
  return E_INVALID_BED_LINE;
}
