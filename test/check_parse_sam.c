#include <stdlib.h>
#include <check.h>
#include <stdio.h>
#include "../src/parse_sam.h"
#include "../src/errors.h"

START_TEST(test_valid_line) {
  char sample_line[1000] =
      "Seq23599_x1\t272\tscaffold_1\t67\t0\t"
      "21M\t*\t0\t0\tGCCACCCATGCCGCATCCACA\t"
      "IIIIIIIIIIIIIIIIIIIII\tAS:i:-6\tXN:i:0\tXM:i:1\tXO:i:0\t"
      "XG:i:0\tNM:i:1\tMD:Z:12A8\tYT:Z:UU\tNH:i:16\tCC:Z:=\t"
      "CP:i:250\tHI:i:0";
  struct sam_entry test_entry;
  int result = parse_line(&test_entry, sample_line);

  ck_assert_int_eq(result, E_SUCCESS);
  ck_assert_int_eq(test_entry.flag, 272);
  ck_assert(test_entry.pos == 67);
  ck_assert_str_eq(test_entry.seq, "GCCACCCATGCCGCATCCACA");
}
END_TEST

START_TEST(test_header_line) {
  char sample_line[100] = "@SQ\tSN:scaffold_1008\tLN:3885";
  struct sam_entry test_entry;
  int result = parse_line(&test_entry, sample_line);

  ck_assert_int_eq(result, E_SAM_HEADER_LINE);
}
END_TEST

START_TEST(test_invalid_line) {
  char sample_line[1000] =
      "Seq23599_x1\t272\tscaffold_1\t6a7\t0\t"
      "21M\t*\t0\t0\tGCCACCCATGCCGCATCCACA\t"
      "IIIIIIIIIIIIIIIIIIIII\tAS:i:-6\tXN:i:0\tXM:i:1\tXO:i:0\t"
      "XG:i:0\tNM:i:1\tMD:Z:12A8\tYT:Z:UU\tNH:i:16\tCC:Z:=\t"
      "CP:i:250\tHI:i:0";
  struct sam_entry test_entry;
  int result = parse_line(&test_entry, sample_line);

  ck_assert_int_eq(result, E_INVALID_SAM_LINE);
}
END_TEST

START_TEST(test_line_without_tabs) {
  char sample_line[100] = "ajksndljka\taksdjnlandkja\tsdlandjkas";
  struct sam_entry test_entry;
  int result = parse_line(&test_entry, sample_line);

  ck_assert_int_eq(result, E_INVALID_SAM_LINE);
}
END_TEST

Suite *hello_suite(void) {
  Suite *s = suite_create("Parse SAM");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_invalid_line);
  tcase_add_test(tc_core, test_valid_line);
  tcase_add_test(tc_core, test_header_line);
  tcase_add_test(tc_core, test_line_without_tabs);
  suite_add_tcase(s, tc_core);

  return s;
}

int main(void) {
  int number_failed;
  Suite *s = hello_suite();
  SRunner *sr = srunner_create(s);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
