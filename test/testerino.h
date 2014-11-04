#include <stddef.h>

#ifndef TESTERINO_H
#define TESTERINO_H

enum _testcode { TEST_PASSED = 0, TEST_FAILED = 1, UNDEFINED = -1 };

struct log_buffer {
  char *buffer;
  char *start;
  size_t capacity;
  size_t free_space;
};

struct test {
  struct log_buffer *log;
  enum _testcode status;
  char *error_msg;
  char *test_msg;
  size_t msg_size;
};

struct test_suite {
  void (**tests)(struct test *t);
  size_t capacity;
  size_t n;
  size_t failed;
};

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_RESET "\x1b[0m"

enum _errorcode { E_SUCCESS = 0, E_FAIL = 1 };
void t_set_msg(struct test *t, const char *msg);
void t_log(struct test *t, const char *msg, ...);
void t_assert(struct test *t, int condition);
void t_assert_msg(struct test *t, int condition, const char *msg);
int suite_add_test(struct test_suite *suite, void (*test)(struct test *));
int suite_run_all_tests(struct test_suite *suite);
int print_test_result(struct test *t);

int create_test_suite(struct test_suite **suite);
int create_test(struct test **t);
int create_buffer(struct log_buffer **log, size_t cap);
int reset_buffer(struct log_buffer *log);
int free_buffer(struct log_buffer *log);
int free_test(struct test *t);
int free_suite(struct test_suite *suite);

#endif