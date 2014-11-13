/**
 *  testerino, a simple unit testing framework for C.
 *
 *   The MIT License (MIT)
 *
 *  Copyright (c) 2014 Michael Huttner
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to
 *  deal in the Software without restriction, including without limitation the
 *  rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 *  sell copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 *  IN THE SOFTWARE.
 */

#include <stddef.h>

#ifndef TESTERINO_H
#define TESTERINO_H

enum _testcode { TEST_PASSED = 0, TEST_FAILED = 1, UNDEFINED = -1 };

struct log_buffer {
  char *buffer; /* Buffer holding log messages, points to curent write point */
  char *start;  /* Always points to the start of the buffer */
  size_t capacity;
  size_t free_space;
};

struct test {
  struct log_buffer *log;
  enum _testcode status;
  char *error_msg; /* message shown if the test fails */
  char *test_msg;  /* message always shown when test is run */
  size_t msg_size; /* maximum size of error and test messages */
};

struct test_suite {
  void (**tests)(struct test *t); /*list of function pointers to tests */
  size_t capacity;
  size_t n;
};

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_RESET "\x1b[0m"

void t_set_msg(struct test *t, const char *msg);
void t_log(struct test *t, const char *msg, ...);
void t_assert(struct test *t, int condition);
void t_assert_msg(struct test *t, int condition, const char *msg);
void t_fail(struct test *t, const char *msg);
int suite_add_test(struct test_suite *suite, void (*test)(struct test *));
int suite_run_all_tests(struct test_suite *suite);
int print_test_result(struct test *t);

int create_test_suite(struct test_suite **suite);
int create_test(struct test **t);
int reset_test(struct test *t);
int create_buffer(struct log_buffer **log, size_t cap);
int reset_buffer(struct log_buffer *log);
int free_buffer(struct log_buffer *log);
int free_test(struct test *t);
int free_suite(struct test_suite *suite);

#endif