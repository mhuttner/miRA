#include "testerino.h"
#include "../src/vfold.h"

void test_reverse_complement(struct test *t) {
  t_set_msg(t, "Testing reverse complement function...");
  struct foldable_sequence s;
  char testseq[] = "GGCAGATTCCCCCTAGACCCGCCCGCACCATGGTCAGGCATGCCCCTCCTCATCGCTGG"
                   "GCACAGCCCAGAGGGT";
  char expected[] = "ACCCTCTGGGCTGTGCCCAGCGATGAGGAGGGGCATGCCTGACCATGGTGCGGGCGGG"
                    "TCTAGGGGGAATCTGCC";
  s.n = strlen(testseq);
  s.seq = (char *)malloc((s.n + 1) * sizeof(char));
  memcpy(s.seq, testseq, s.n + 1);
  s.seq[s.n] = 0;

  reverse_complement(&s);
  t_log(t, "Got:      %s\n", s.seq);
  t_log(t, "Expected: %s\n", expected);
  t_assert_msg(t, strncmp(s.seq, expected, s.n) == 0,
               "Reversing the Sequence failed");
  // t_fail(t, "always");
  free(s.seq);
}