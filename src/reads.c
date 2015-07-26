#include "reads.h"
#include "errors.h"
#include "defs.h"
#include "coverage.h"
#include "parse_sam.h"

int count_unique_reads(struct extended_candidate *ecand, struct sam_file *sam) {

  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence_list *css_list = ecand->possible_micro_rnas;

  struct candidate_subsequence *mature_mirna = NULL;
  struct candidate_subsequence *star_mirna = NULL;
  struct sam_entry *entry = NULL;

  ecand->total_reads = 0;

  for (size_t j = 0; j < css_list->n; j++) {
    mature_mirna = css_list->mature_sequences[j];
    star_mirna = mature_mirna->matching_sequence;
    create_unique_read_list(&mature_mirna->reads);
    create_unique_read_list(&star_mirna->reads);
  }

  for (size_t i = 0; i < sam->n; i++) {
    entry = sam->entries[i];
    long entry_start = entry->pos - 1;
    char strand = (entry->flag & REV_COMPLM) ? '-' : '+';
    for (size_t j = 0; j < css_list->n; j++) {
      mature_mirna = css_list->mature_sequences[j];
      star_mirna = mature_mirna->matching_sequence;
      if (check_subsequence_match(entry, cand, mature_mirna)) {
        if (strand == '-') {
          char *reversed = NULL;
          reverse_complement_sequence_string(&reversed, entry->seq,
                                             strlen(entry->seq) + 1);
          free(entry->seq);
          entry->seq = reversed;
          strand = '+';
        }
        add_read_to_unique_read_list(mature_mirna->reads, entry_start,
                                     entry->seq);
      }
      if (check_subsequence_match(entry, cand, star_mirna)) {
        if (strand == '-') {
          char *reversed = NULL;
          reverse_complement_sequence_string(&reversed, entry->seq,
                                             strlen(entry->seq) + 1);
          free(entry->seq);
          entry->seq = reversed;
          strand = '+';
        }
        add_read_to_unique_read_list(star_mirna->reads, entry_start,
                                     entry->seq);
      }
    }
    if (entry_start >= cand->start) {
      u64 end = entry_start + strlen(entry->seq);
      if (end <= cand->end) {
        ecand->total_reads++;
      }
    }
  }
  ecand->total_read_percent = (double)ecand->total_reads / sam->n;
  return E_SUCCESS;
}

int check_subsequence_match(struct sam_entry *entry,
                            struct micro_rna_candidate *cand,
                            struct candidate_subsequence *sseq) {
  const int READCOUNT_FLANK = 30;
  size_t start = cand->start + sseq->start - READCOUNT_FLANK;
  size_t end = cand->start + sseq->end + READCOUNT_FLANK;
  long entry_start = entry->pos - 1;
  char strand = (entry->flag & REV_COMPLM) ? '-' : '+';
  if (strand != cand->strand) {
    return 0;
  }
  if (strcmp(entry->rname, cand->chrom) != 0) {
    return 0;
  }
  if (entry_start >= start) {
    u64 entry_end = entry_start + strlen(entry->seq);
    if (entry_end <= end) {
      return 1;
    }
  }
  return 0;
}

int create_unique_read(struct unique_read **read, u64 start, const char *seq) {
  struct unique_read *read_tmp =
      (struct unique_read *)malloc(sizeof(struct unique_read));
  if (read_tmp == NULL) {
    return E_MALLOC_FAIL;
  }
  read_tmp->start = start;
  size_t l = strlen(seq);
  read_tmp->seq = (char *)malloc((l + 1) * sizeof(char));
  if (read_tmp->seq == NULL) {
    free(read_tmp);
    return E_MALLOC_FAIL;
  }
  memcpy(read_tmp->seq, seq, l);
  read_tmp->seq[l] = 0;
  read_tmp->end = start + l;
  read_tmp->count = 1;
  *read = read_tmp;
  return E_SUCCESS;
}
int create_unique_read_list(struct unique_read_list **ur_list) {
  const int INITIAL_CAPACITY = 16;
  struct unique_read_list *ur_list_tmp =
      (struct unique_read_list *)malloc(sizeof(struct unique_read_list));
  if (ur_list_tmp == NULL) {
    return E_MALLOC_FAIL;
  }
  ur_list_tmp->capacity = INITIAL_CAPACITY;
  ur_list_tmp->reads = (struct unique_read **)malloc(
      ur_list_tmp->capacity * sizeof(struct unique_read *));
  if (ur_list_tmp->reads == NULL) {
    free(ur_list_tmp);
    return E_MALLOC_FAIL;
  }
  ur_list_tmp->n = 0;
  *ur_list = ur_list_tmp;
  return E_SUCCESS;
}
int append_unique_read_list(struct unique_read_list *ur_list,
                            struct unique_read *read) {
  if (ur_list->n >= ur_list->capacity) {
    ur_list->capacity *= 2;
    struct unique_read **tmp = (struct unique_read **)realloc(
        ur_list->reads, ur_list->capacity * sizeof(struct unique_read *));
    if (tmp == NULL) {
      return E_REALLOC_FAIL;
    }
    ur_list->reads = tmp;
  }
  ur_list->reads[ur_list->n] = read;
  ur_list->n++;
  return E_SUCCESS;
}
int add_read_to_unique_read_list(struct unique_read_list *ur_list, u64 start,
                                 const char *seq) {
  struct unique_read *read = NULL;
  for (size_t i = 0; i < ur_list->n; i++) {
    read = ur_list->reads[i];
    if (start == read->start && strcmp(seq, read->seq) == 0) {
      read->count++;
      return E_SUCCESS;
    }
  }
  int err;
  err = create_unique_read(&read, start, seq);
  if (err != E_SUCCESS) {
    return err;
  }
  err = append_unique_read_list(ur_list, read);
  if (err != E_SUCCESS) {
    return err;
  }
  return E_SUCCESS;
}

int free_unique_read_list(struct unique_read_list *ur_list) {
  if (ur_list == NULL) {
    return E_SUCCESS;
  }
  for (size_t i = 0; i < ur_list->n; i++) {
    free_unique_read(ur_list->reads[i]);
  }
  free(ur_list->reads);
  free(ur_list);
  return E_SUCCESS;
}
int free_unique_read(struct unique_read *read) {
  free(read->seq);
  free(read);
  return E_SUCCESS;
}