#include <stdio.h>
#include "testerino.h"
#include "sample_tests.h"

int main(int argc, char const *argv[]) {
  struct test_suite *s = NULL;
  create_test_suite(&s);
  suite_add_test(s, test_nothing);
  suite_add_test(s, test_nothing_and_fail);
  suite_add_test(s, test_failing_with_long_log);
  suite_add_test(s, test_long_test_msg);
  suite_run_all_tests(s);
  free_suite(s);
  return 0;
}