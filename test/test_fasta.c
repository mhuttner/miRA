#include "testerino.h"
#include "../src/fasta.h"
#include "../src/errors.h"
#include "../src/uthash.h"
#include <string.h>

void test_strip_newlines(struct test *t) {
  t_set_msg(t, "Testing if newline chars get stripped correctly...");
  char sample_line[100] = "this\nis\n\na\ntest";
  int size = strlen(sample_line);
  char *result = NULL;
  size_t result_size;

  strip_newlines(&result, &result_size, sample_line, size);

  t_assert_msg(t, result_size == 11, "Wrong size of new String");
  t_assert_msg(t, strncmp(result, "thisisatest", 11) == 0, "Stripping failed");
}

void test_read_fasta_file(struct test *t) {
  t_set_msg(t, "Testing reading a fasta File...");
  struct genome_sequence *sequence_table = NULL;
  char filename[30] = "test/test.fasta";
  int err = read_fasta_file(&sequence_table, filename);
  t_assert_msg(t, err == E_SUCCESS, "Parsing failed");
  t_log(t, "Error: %d\n", err);

  struct genome_sequence *s = NULL;
  struct genome_sequence *tmp = NULL;
  HASH_ITER(hh, sequence_table, s, tmp) {
    s->data[s->n] = 0;
    t_log(t, "%s %s \n\n", s->chrom, s->data);
  }

  HASH_FIND_STR(sequence_table, "HSBGPG", s);
  if (s == NULL) {
    t_fail(t, "A inserted entry was not found");
    return;
  }
  t_assert_msg(t, strncmp(s->data, "GGCAGATTCCCCCTAGACCCGCCC", 20) == 0,
               "Entry was read wrong");

  free_sequence_table(sequence_table);
}