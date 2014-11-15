
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "parse_sam.h"
#include "errors.h"

int parse_sam(struct sam_file **sam, char *file) {
  static const int MAXLINELENGHT = 2048;
  static const int STARTINGSIZE = 1024;
  struct sam_file *data = (struct sam_file *)malloc(sizeof(struct sam_file));
  if (data == NULL) {
    return E_MALLOC_FAIL;
  }
  data->capacity = STARTINGSIZE;
  data->entries =
      (struct sam_entry **)malloc(data->capacity * sizeof(struct sam_entry*));
  if (data->entries == NULL) {
    free(data);
    return E_MALLOC_FAIL;
  }
  data->n = 0;

  data->header_cap = STARTINGSIZE;
  data->headers =
      (struct sq_header **)malloc(data->header_cap * sizeof(struct sq_header*));
  if (data->headers == NULL) {
    free(data->entries);
    free(data);
    return E_MALLOC_FAIL;
  }
  data->header_n = 0;

  FILE *fp = fopen(file, "r");
  if (fp == NULL) {
    free(data->entries);
    free(data->headers);
    free(data);
    return E_FILE_NOT_FOUND;
  }
  char line[MAXLINELENGHT];
  while (fgets(line, sizeof(line), fp) != NULL) {
    int result = parse_line(data->entries +data->n , line);
    if (result == E_SAM_HEADER_LINE) {
      int err = parse_header(data->headers +data->header_n, line);
      if (err != E_SUCCESS)
        continue;
      data->header_n++;
      if (data->header_n == data->header_cap) {
        data->header_cap *= 2;
        struct sq_header **tmph = (struct sq_header **)realloc(
            data->headers, data->header_cap * sizeof(struct sq_header*));
        if (tmph == NULL) {
          free_sam(data);
          fclose(fp);
          return E_REALLOC_FAIL;
        }
        data->headers = tmph;
      }
      continue;
    }
    if (result != E_SUCCESS)
      continue;
    data->n++;
    if (data->n == data->capacity) {
      data->capacity *= 2;
      struct sam_entry **tmp = (struct sam_entry **)realloc(
          data->entries, data->capacity * sizeof(struct sam_entry*));
      if (tmp == NULL) {
        free_sam(data);
        fclose(fp);
        return E_REALLOC_FAIL;
      }
      data->entries = tmp;
    }
  }
  fclose(fp);
  *sam = data;

  return E_SUCCESS;
}
int parse_line(struct sam_entry **entry, char *line) {
  const char seperator = '\t';
  const int num_entries = 11;
  const char header_line_marker = '@';

  struct sam_entry *e = (struct sam_entry*)malloc(sizeof(sam_entry));



  char *start = line;
  char *end = NULL;
  if (line[0] == header_line_marker)
    return E_SAM_HEADER_LINE;
  char **tokens = (char **)malloc(num_entries * sizeof(char *));
  if (tokens == NULL) {
    return E_MALLOC_FAIL;
  }
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
  int line_valid = 1;

  if (current_token < num_entries) {
    for (int i = 0; i < current_token; i++) {
      free(tokens[i]);
    }
    free(tokens);
    return E_INVALID_SAM_LINE;
  }

  char *check = NULL;

  e->qname = tokens[0];
  e->rname = tokens[2];
  e->cigar = tokens[5];
  e->rnext = tokens[6];
  e->seq = tokens[9];
  e->qual = tokens[10];

  long flag = strtol(tokens[1], &check, 10);
  if (check == tokens[1] || *check != 0) {
    line_valid = 0;
  }
  e->flag = (int)flag;
  long pos = strtol(tokens[3], &check, 10);
  if (check == tokens[3] || *check != 0) {
    line_valid = 0;
  }
  e->pos = pos;
  long mapq = strtol(tokens[4], &check, 10);
  if (check == tokens[4] || *check != 0) {
    line_valid = 0;
  }
  e->mapq = (int)mapq;
  long pnext = strtol(tokens[7], &check, 10);
  if (check == tokens[7] || *check != 0) {
    line_valid = 0;
  }
  e->pnext = pnext;
  long tlen = strtol(tokens[8], &check, 10);
  if (check == tokens[8] || *check != 0) {
    line_valid = 0;
  }
  e->tlen = tlen;

  if (line_valid) {
    free(tokens[1]);
    free(tokens[3]);
    free(tokens[4]);
    free(tokens[7]);
    free(tokens[8]);
    free(tokens);
    *entry = e;
    return E_SUCCESS;
  } else {
    for (int i = 0; i < num_entries; i++) {
      free(tokens[i]);
    }
    free(tokens);
    return E_INVALID_SAM_LINE;
  }
}

int parse_header(struct sq_header **header, char *line) {
  const char *sq_header_marker = "@SQ";
  const char *sq_sn_marker = "SN:";
  const char *sq_ln_marker = "LN:";
  const int marker_length = 3;
  const char seperator = '\t';

  struct sq_header *h = (struct sq_header*)malloc(sizeof(sq_header));

  if (strncmp(line, sq_header_marker, marker_length) != 0) {
    return E_SAM_NON_SQ_HEADER;
  }
  char *start = strchr(line, seperator) + 1;
  if (strncmp(start, sq_sn_marker, marker_length) != 0) {
    return E_INVALID_SAM_LINE;
  }
  start += marker_length;
  char *end = strchr(start, seperator);

  int l = end - start;
  char *sn = (char *)malloc((l + 1) * sizeof(char));
  if (sn == NULL) {
    return E_MALLOC_FAIL;
  }
  memcpy(sn, start, l);
  sn[l] = 0;

  start = end + 1;
  if (strncmp(start, sq_ln_marker, marker_length) != 0) {
    free(sn);
    return E_INVALID_SAM_LINE;
  }
  start += marker_length;
  end = strchr(start, seperator);
  if (end == NULL) {
    end = strchr(start, '\n');
    if (end == NULL) {
      free(sn);
      return E_INVALID_SAM_LINE;
    }
  }
  char *check = NULL;
  long ln = strtol(start, &check, 10);
  if (check != end) {
    free(sn);
    return E_INVALID_SAM_LINE;
  }
  h->sn = sn;
  h->ln = ln;
  *header = h;
  return E_SUCCESS;
}

int free_sam(struct sam_file *sam) {
  for (size_t i = 0; i < sam->n; i++) {
    free_sam_entry(sam->entries[i]);
  }
  for (size_t i = 0; i < sam->header_n; i++) {
    free_sam_header(sam->headers[i]);
  }
  free(sam->entries);
  free(sam->headers);
  free(sam);
  return E_SUCCESS;
}

int free_sam_entry(struct sam_entry *e) {
  free(e->qname);
  free(e->rname);
  free(e->cigar);
  free(e->rnext);
  free(e->seq);
  free(e->qual);
  free(e);
  return E_SUCCESS;
};

int free_sam_header(struct sam_header* h){
  free(h->sn);
  free(h);
  return E_SUCCESS;
}
