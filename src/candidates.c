#include "candidates.h"
#include "vfold.h"
#include "errors.h"
#include "cluster.h"
#include <stdlib.h>

int create_empty_candidate_list(struct candidate_list **cand_list) {
  const int INITIAL_CAPACITY = 128;
  struct candidate_list *tmp_list =
      (struct candidate_list *)malloc(sizeof(struct candidate_list));
  if (tmp_list == NULL) {
    return E_MALLOC_FAIL;
  }
  tmp_list->capacity = INITIAL_CAPACITY;
  tmp_list->n = 0;
  tmp_list->candidates = (struct micro_rna_candidate **)malloc(
      tmp_list->capacity * sizeof(struct micro_rna_candidate *));
  if (tmp_list->candidates == NULL) {
    free(tmp_list);
    return E_MALLOC_FAIL;
  }
  *cand_list = tmp_list;
  return E_SUCCESS;
}

int add_candidate_to_list(struct candidate_list *cand_list,
                          struct micro_rna_candidate *cand) {
  if (cand_list->n >= cand_list->capacity) {
    cand_list->capacity *= 2;
    struct micro_rna_candidate **tmp = (struct micro_rna_candidate **)realloc(
        cand_list->candidates,
        cand_list->capacity * sizeof(struct micro_rna_candidate *));
    if (tmp == NULL) {
      return E_REALLOC_FAIL;
    }
    cand_list->candidates = tmp;
  }
  cand_list->candidates[cand_list->n] = cand;
  cand_list->n++;
  return E_SUCCESS;
}

int convert_seq_list_to_cand_list(struct candidate_list **cand_list,
                                  struct sequence_list *seq_list) {
  struct candidate_list *tmp_cand_list = NULL;
  create_empty_candidate_list(&tmp_cand_list);
  struct foldable_sequence *fs = NULL;
  struct micro_rna_candidate *cand = NULL;
  for (size_t i = 0; i < seq_list->n; i++) {
    fs = seq_list->sequences[i];
    if (fs->structure->is_valid == 0) {
      continue;
    }
    create_micro_rna_candidate(&cand, fs);
    add_candidate_to_list(tmp_cand_list, cand);
  }
  *cand_list = tmp_cand_list;
  return E_SUCCESS;
};

int create_micro_rna_candidate(struct micro_rna_candidate **cand,
                               struct foldable_sequence *fs) {
  struct cluster *c = fs->c;
  struct structure_info *si = fs->structure;
  if (si == NULL) {
    return E_NO_STRUCTURE;
  }
  if (si->is_valid == 0) {
    return E_STRUCTURE_IS_INVALID;
  }

  long structure_start;
  long structure_end;
  if (c->strand == '+') {
    /*0 based */
    structure_start = si->start - 1 + c->flank_start;
    structure_end = structure_start + si->n;
  } else {
    /*0 based */
    structure_end = c->flank_end - (si->start - 1);
    structure_start = structure_end - si->n;
  }
  char *structure_sequence = (char *)malloc((si->n + 1) * sizeof(char));
  if (structure_sequence == NULL) {
    return E_MALLOC_FAIL;
  }
  memcpy(structure_sequence, fs->seq + si->start - 1, si->n);
  structure_sequence[si->n] = 0;
  char *structure_copy = (char *)malloc((si->n + 1) * sizeof(char));
  if (structure_copy == NULL) {
    free(structure_sequence);
    return E_MALLOC_FAIL;
  }
  memcpy(structure_copy, si->structure_string, si->n);
  structure_copy[si->n] = 0;

  size_t tmp_n = strlen(c->chrom);
  char *chrom_copy = (char *)malloc((tmp_n + 1) * sizeof(char));
  if (chrom_copy == NULL) {
    free(structure_sequence);
    free(structure_copy);
    return E_MALLOC_FAIL;
  }
  memcpy(chrom_copy, c->chrom, tmp_n);
  chrom_copy[tmp_n] = 0;

  struct micro_rna_candidate *cand_tmp =
      (struct micro_rna_candidate *)malloc(sizeof(struct micro_rna_candidate));
  if (cand_tmp == NULL) {
    free(structure_sequence);
    free(structure_copy);
    return E_MALLOC_FAIL;
  }
  cand_tmp->id = c->id;
  cand_tmp->chrom = chrom_copy;
  cand_tmp->start = structure_start;
  cand_tmp->end = structure_end;
  cand_tmp->sequence = structure_sequence;
  cand_tmp->structure = structure_copy;
  cand_tmp->mfe = si->mfe;
  cand_tmp->pvalue = si->pvalue;
  cand_tmp->mean = si->mean;
  cand_tmp->sd = si->sd;

  *cand = cand_tmp;
  return E_SUCCESS;
}
int write_candidate_file(struct candidate_list *cand_list, char *filename) {
  FILE *fp = NULL;
  if (filename != NULL) {
    fp = fopen(filename, "w");
    if (fp == NULL) {
      return E_FILE_WRITING_FAILED;
    }
  }
  struct micro_rna_candidate *cand = NULL;
  for (size_t i = 0; i < cand_list->n; i++) {
    cand = cand_list->candidates[i];
    write_candidate_line(fp, cand);
  }
  if (fp != NULL) {
    fclose(fp);
  }
  return E_SUCCESS;
}

int write_candidate_line(FILE *fp, struct micro_rna_candidate *cand) {
  if (fp == NULL) {
    fp = stdout;
  }
  fprintf(fp, "Cluster_%lld_%s\t", cand->id,
          (cand->strand == '-') ? "minus" : "plus");
  fprintf(fp, "%lld\t", cand->id);
  fprintf(fp, "%s\t", cand->chrom);
  fprintf(fp, "%c\t", cand->strand);
  fprintf(fp, "%lld\t", cand->start);
  fprintf(fp, "%lld\t", cand->end);
  fprintf(fp, "%s\t", cand->sequence);
  fprintf(fp, "%s\t", cand->structure);
  fprintf(fp, "%7.5f\t", cand->mfe);
  fprintf(fp, "%9.7e\t", cand->pvalue);
  fprintf(fp, "%7.5e\t", cand->mean);
  fprintf(fp, "%7.5e\n", cand->sd);
  return E_SUCCESS;
}
int read_candidate_file(struct candidate_list **cand_list, char *filename) {
  static const int MAXLINELENGHT = 2048;
  struct candidate_list *tmp_list = NULL;
  int err = create_empty_candidate_list(&tmp_list);
  if (err) {
    return err;
  }
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    free_candidate_list(tmp_list);
    return E_FILE_NOT_FOUND;
  }
  char line[MAXLINELENGHT];
  struct micro_rna_candidate *cand = NULL;
  while (fgets(line, sizeof(line), fp) != NULL) {
    err = parse_candidate_line(&cand, line);
    if (err) {
      free_candidate_list(tmp_list);
      return err;
    }
    add_candidate_to_list(tmp_list, cand);
  }
  fclose(fp);

  return E_SUCCESS;
}

int parse_candidate_line(struct micro_rna_candidate **cand, char *line) {
  const char seperator = '\t';
  const int num_entries = 12;

  struct micro_rna_candidate *tmp_cand =
      (struct micro_rna_candidate *)malloc(sizeof(struct micro_rna_candidate));
  if (tmp_cand == NULL) {
    return E_MALLOC_FAIL;
  }
  char **tokens = (char **)malloc(num_entries * sizeof(char *));
  if (tokens == NULL) {
    free(tmp_cand);
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
    return E_INVALID_CANDIDATE_LINE;
  }
  char *check = NULL;
  tmp_cand->id = strtol(tokens[1], &check, 10);
  if (check == tokens[1] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[1]);
  tokens[1] = NULL;

  tmp_cand->chrom = tokens[2];
  tmp_cand->strand = *tokens[3];
  free(tokens[3]);
  tokens[3] = NULL;
  tmp_cand->start = strtol(tokens[4], &check, 10);
  if (check == tokens[4] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[4]);
  tokens[4] = NULL;
  tmp_cand->end = strtol(tokens[5], &check, 10);
  if (check == tokens[5] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[5]);
  tokens[5] = NULL;
  tmp_cand->sequence = tokens[6];
  tmp_cand->structure = tokens[7];
  tmp_cand->mfe = strtod(tokens[8], &check);
  if (check == tokens[8] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[8]);
  tokens[8] = NULL;
  tmp_cand->pvalue = strtod(tokens[9], &check);
  if (check == tokens[9] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[9]);
  tokens[9] = NULL;

  tmp_cand->mean = strtod(tokens[10], &check);
  if (check == tokens[10] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[10]);
  tokens[10] = NULL;

  tmp_cand->sd = strtod(tokens[11], &check);
  if (check == tokens[11] || *check != 0) {
    goto line_invalid;
  }
  free(tokens[11]);
  tokens[11] = NULL;
  *cand = tmp_cand;
  return E_SUCCESS;
line_invalid:
  for (int i = 0; i < num_entries; i++) {
    if (tokens[i] != NULL) {
      free(tokens[i]);
    }
  }
  free(tokens);
  return E_INVALID_CANDIDATE_LINE;
}

int free_candidate_list(struct candidate_list *cand_list) {
  for (size_t i = 0; i < cand_list->n; i++) {
    free_micro_rna_candidate(cand_list->candidates[i]);
  }
  free(cand_list);
  return E_SUCCESS;
}

int free_micro_rna_candidate(struct micro_rna_candidate *cand) {
  free(cand->chrom);
  free(cand->sequence);
  free(cand->structure);
  free(cand);
  return E_SUCCESS;
}
