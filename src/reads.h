

#ifndef READS_H
#define READS_H

#include "coverage.h"
#include "util.h"
#include "defs.h"
#include "parse_sam.h"
#include "candidates.h"

/* forward declaration of struct in coverage.h */
struct extended_candidate;
struct candidate_subsequence;

struct unique_read {
  u64 start;
  u64 end;
  char *seq;
  u32 count;
};
struct unique_read_list {
  struct unique_read **reads;
  size_t n;
  size_t capacity;
};

int count_unique_reads(struct extended_candidate *ecand, struct sam_file *sam);
int check_subsequence_match(struct sam_entry *entry,
                            struct micro_rna_candidate *cand,
                            struct candidate_subsequence *sseq);
int create_unique_read(struct unique_read **read, u64 start, const char *seq);
int create_unique_read_list(struct unique_read_list **ur_list);
int append_unique_read_list(struct unique_read_list *ur_list,
                            struct unique_read *read);
int add_read_to_unique_read_list(struct unique_read_list *ur_list, u64 start,
                                 const char *seq);

int free_unique_read_list(struct unique_read_list *ur_list);
int free_unique_read(struct unique_read *read);

#endif