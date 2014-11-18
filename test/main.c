#include <stdio.h>
#include "testerino.h"
#include "test_parse_sam.h"
#include "test_cluster.h"

int main(int argc, char const *argv[]) {
  struct test_suite *s = NULL;
  create_test_suite(&s);
  suite_add_test(s, test_valid_line);
  suite_add_test(s, test_header_line);
  suite_add_test(s, test_invalid_line);
  suite_add_test(s, test_multiple_consecutive_tabs);
  suite_add_test(s, test_sort_clusters);
  suite_add_test(s, test_merge_clusters);
  suite_add_test(s, test_filter_clusters);
  suite_add_test(s, test_merge_extended_clusters);
  suite_run_all_tests(s);
  free_suite(s);
  return 0;
}