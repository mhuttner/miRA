#include "testerino.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../src/parse_sam.h"
#include "../src/errors.h"

void test_valid_line(struct test *t) {
  t_set_msg(t, "Testing reading a valid Sam line...");
  char sample_line[1000] =
      "Seq23599_x1\t272\tscaffold_1\t67\t0\t"
      "21M\t*\t0\t0\tGCCACCCATGCCGCATCCACA\t"
      "IIIIIIIIIIIIIIIIIIIII\tAS:i:-6\tXN:i:0\tXM:i:1\tXO:i:0\t"
      "XG:i:0\tNM:i:1\tMD:Z:12A8\tYT:Z:UU\tNH:i:16\tCC:Z:=\t"
      "CP:i:250\tHI:i:0\n";
  struct sam_entry test_entry;
  int result = parse_line(&test_entry, sample_line);
  t_assert_msg(t, result == E_SUCCESS, "Parsing the line failed");
  t_assert_msg(t, test_entry.flag == 272, "Flag parsed wrong");
  t_assert_msg(t, test_entry.pos == 67, "Position parsed wrong");
  t_log(t, "%s", test_entry.seq);
  t_assert_msg(t, strcmp(test_entry.seq, "GCCACCCATGCCGCATCCACA") == 0,
               "Sequence parsed wrong");
  free_sam_entry(&test_entry);
}

void test_header_line(struct test *t) {
  t_set_msg(t, "Testing reading a header Sam line...");
  char sample_line[100] = "@SQ\tSN:scaffold_1008\tLN:3885";
  struct sam_entry test_entry;
  int result = parse_line(&test_entry, sample_line);
  t_assert_msg(t, result == E_SAM_HEADER_LINE, "Header line not regognized");
}

void test_invalid_line(struct test *t) {
  t_set_msg(t, "Testing reading an invalid Sam line...");
  char sample_line[1000] =
      "Seq23599_x1\t272\tscaffold_1\t6a7\t0\t"
      "21M\t*\t0\t0\tGCCACCCATGCCGCATCCACA\t"
      "IIIIIIIIIIIIIIIIIIIII\tAS:i:-6\tXN:i:0\tXM:i:1\tXO:i:0\t"
      "XG:i:0\tNM:i:1\tMD:Z:12A8\tYT:Z:UU\tNH:i:16\tCC:Z:=\t"
      "CP:i:250\tHI:i:0";
  struct sam_entry test_entry;
  int result = parse_line(&test_entry, sample_line);
  t_assert_msg(t, result == E_INVALID_SAM_LINE,
               "Did not fail on invalid entry");
}

void test_multiple_consecutive_tabs(struct test *t) {
  t_set_msg(t, "Testing reading a valid Sam line with consecutive tabs...");
  char sample_line[1000] =
      "Seq23599_x1\t\t272\t\tscaffold_1\t67\t0\t"
      "21M\t*\t0\t0\t\t\t\tGCCACCCATGCCGCATCCACA\t"
      "IIIIIIIIIIIIIIIIIIIII\tAS:i:-6\tXN:i:0\tXM:i:1\tXO:i:0\t"
      "XG:i:0\tNM:i:1\tMD:Z:12A8\tYT:Z:UU\tNH:i:16\tCC:Z:=\t"
      "CP:i:250\tHI:i:0\n";
  struct sam_entry test_entry;
  int result = parse_line(&test_entry, sample_line);
  t_assert_msg(t, result == E_SUCCESS, "Parsing the line failed");
  t_assert_msg(t, test_entry.flag == 272, "Flag parsed wrong");
  t_assert_msg(t, test_entry.pos == 67, "Position parsed wrong");
  t_assert_msg(t, strcmp(test_entry.seq, "GCCACCCATGCCGCATCCACA") == 0,
               "Sequence parsed wrong");
  free_sam_entry(&test_entry);
}