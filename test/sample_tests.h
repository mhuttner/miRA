#ifndef SAMPLE_TESTS_H
#define SAMPLE_TESTS_H

void test_nothing(struct test *t);
void test_nothing_and_fail(struct test *t);
void test_failing_with_long_log(struct test *t);
void test_long_test_msg(struct test *t);

#endif