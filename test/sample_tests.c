#include "testerino.h"

void test_nothing(struct test *t) {
  t_set_msg(t, "Testing nothing...");
  t_log(t, "This test should always succeed");
  t_assert(t, 1 == 1);
}

void test_nothing_and_fail(struct test *t) {
  t_set_msg(t, "Testing failing at nothing...");
  t_log(t, "This test should always fail");
  t_assert_msg(t, 1 == 2, "1 does not equal 2");
}

void test_failing_with_long_log(struct test *t) {
  t_set_msg(t, "Testing failing with a long log...");
  t_log(t, "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do "
           "eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut "
           "enim ad minim veniam, quis nostrud exercitation ullamco laboris "
           "nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor "
           "in reprehenderit in voluptate velit esse cillum dolore eu fugiat "
           "nulla pariatur. Excepteur sint occaecat cupidatat non proident, "
           "sunt in culpa qui officia deserunt mollit anim id est laborum.");
  t_assert_msg(t, 1 == 2, "this is expected");
}

void test_long_test_msg(struct test *t) {
  t_set_msg(t, "Testing a very long test message to see if spacing works "
               "correcly... bla bal bla");
  t_log(t, "This test should always succeed");
  t_assert(t, 1 == 1);
}