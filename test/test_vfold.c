#include "testerino.h"
#include "../src/vfold.h"

void test_reverse_complement(struct test *t) {
  t_set_msg(t, "Testing reverse complement function...");
  struct foldable_sequence s;
  char testseq[] = "GGCAGATTCCCCCTAGACCCGCCCGCACCATGGTCAGGCATGCCCCTCCTCATCGCTGG"
                   "GCACAGCCCAGAGGGT";
  char expected[] = "ACCCTCTGGGCTGTGCCCAGCGATGAGGAGGGGCATGCCTGACCATGGTGCGGGCGGG"
                    "TCTAGGGGGAATCTGCC";
  s.n = strlen(testseq) + 1;
  s.seq = (char *)malloc((s.n) * sizeof(char));
  memcpy(s.seq, testseq, s.n);
  s.seq[s.n - 1] = 0;

  reverse_complement(&s);
  t_log(t, "Got:      %s\n", s.seq);
  t_log(t, "Expected: %s\n", expected);
  t_assert_msg(t, strncmp(s.seq, expected, s.n) == 0,
               "Reversing the Sequence failed");
  // t_fail(t, "always");
  free(s.seq);
}

void test_folding(struct test *t) {
  t_set_msg(t, "Testing libRNA folding...");
  char *fake_argv[] = {"fold", "test/data/contigs.bed", "test/data/Chlre3.fa"};
  int fake_argc = 3;
  vfold(fake_argc, fake_argv);
  t_fail(t, "always");
}