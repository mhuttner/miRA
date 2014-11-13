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
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "testerino.h"

static const int E_SUCCESS = 0;
static const int E_FAIL = 1;

void t_set_msg(struct test *t, const char *msg) {
  snprintf(t->test_msg, t->msg_size, "%s", msg);
}
/**
 * \brief Conditional logging in a test case.
 *
 * Works like printf in a test case, but the output is only written to console
 * if the test fails.
 *
 * Typical usage:
 * \code
 *   t_log(t,"Info i need when the test fails %d",42);
 * \endcode
 *
 * \param t test data passed to test function
 * \param msg format string
 * \param ... Arguments of format string
 *
 */
void t_log(struct test *t, const char *msg, ...) {
  struct log_buffer *log = t->log;
  if (log->free_space <= 0)
    return;
  int n;
  va_list fmtargs;
  va_start(fmtargs, msg);
  n = vsnprintf(log->buffer, log->free_space, msg, fmtargs);
  if (n > 0) {
    if (log->free_space < n) {
      log->capacity *= 2;
      char *tmp_buffer =
          (char *)realloc(log->start, log->capacity * sizeof(char));
      if (tmp_buffer == NULL) {
        log->free_space = 0;
      } else {
        log->free_space += log->capacity;
        long offset = log->buffer - log->start;
        log->start = tmp_buffer;
        log->buffer = log->start + offset;
      }
    } else {
      log->buffer += n;
      log->free_space -= n;
    }
  }
  va_end(fmtargs);
}

/**
 * \brief Asserting conditions in test.
 *
 * Assert boolean conditions in a test case, if one assert fails the test fails.
 *
 * Typical usage:
 * \code
 *   t_assert(t,1==2)
 * \endcode
 *
 * \param t test data passed to test function
 * \param condition boolean condition, if false assert has failed
 *
 */
void t_assert(struct test *t, int condition) {
  t_assert_msg(t, condition, "Assert failed");
}
void t_assert_msg(struct test *t, int condition, const char *msg) {
  if (t->status == UNDEFINED) {
    t->status = TEST_PASSED;
  }
  if (!condition) {
    t->status = TEST_FAILED;
    snprintf(t->error_msg, t->msg_size, "%s", msg);
  }
}

/**
 * \brief fails the test.
 *
 * Fails the Test case without condition, use e.g if the test setup has failed.
 *
 * Typical usage:
 * \code
 *   t_fail(t,"This test failed")
 * \endcode
 *
 * \param t test data passed to test function
 * \param msg error mesage to be given for fail
 *
 */
void t_fail(struct test *t, const char *msg) { t_assert_msg(t, 0, msg); }

/**
 * \brief add a test to the testsuite.
 *
 * Adds the function pointer of the test to the list of test cases in this
 * suite.
 * One test can be added multiple times and to different test suites.
 *
 * Typical usage:
 * \code
 *   suite_add_test(suite,test_function);
 * \endcode
 *
 * \param suite test suite
 * \param test function pointer to test function
 *
 * \returns E_SUCCESS on success, E_FAIL if memory allocation fails
 */
int suite_add_test(struct test_suite *suite, void (*test)(struct test *)) {
  if (suite->n >= suite->capacity) {
    suite->capacity *= 2;
    void (**tmpf)(struct test *t) = (void (**)(struct test *))realloc(
        suite->tests, suite->capacity * sizeof(void (*)(struct test *)));
    if (tmpf == NULL) {
      return E_FAIL;
    }
    suite->tests = tmpf;
  }
  suite->tests[suite->n] = test;
  suite->n++;
  return E_SUCCESS;
}

/**
 * \brief run all test cases in this suite.
 *
 * Sets up testing by creating an empty test struct and passing it to every test
 * in the order they were added.
 * All tests are executed even if some fail.
 *
 * Typical usage:
 * \code
 *   suite_run_all_tests(suite);
 * \endcode
 *
 * \param suite test suite
 *
 * \returns E_SUCCESS on success, E_FAIL if create_test fails.
 */
int suite_run_all_tests(struct test_suite *suite) {
  struct test *t = NULL;
  int err;
  err = create_test(&t);
  if (err) {
    return E_FAIL;
  }
  for (int i = 0; i < suite->n; i++) {
    void (*test)(struct test *) = suite->tests[i];
    test(t);
    print_test_result(t);
    reset_test(t);
  }
  free_test(t);
  return E_SUCCESS;
}

/**
 * \brief pretty prints the results of a test.
 *
 * If a test passed it's message is printed with an [OK] netxt to it.
 * If a test failed it's its message, error message and log are printed.
 * Used internally by suite_run_all_tests(...)
 *
 * \param t test data
 *
 * \returns E_SUCCESS currently always succedes.
 */
int print_test_result(struct test *t) {
  int free_space = 70 - strnlen(t->test_msg, t->msg_size);
  while (free_space < 0) {
    free_space += 80;
  }
  if (t->status == TEST_PASSED) {
    printf("%s %*s" ANSI_COLOR_GREEN "[OK]" ANSI_COLOR_RESET "\n", t->test_msg,
           free_space, "");
  } else {
    printf("%s %*s" ANSI_COLOR_RED "[FAIL]" ANSI_COLOR_RESET "\n", t->test_msg,
           free_space, "");
    printf("%s\n\n", t->error_msg);
    printf("%s\n\n", t->log->start);
  }
  return E_SUCCESS;
}

/**
 * \brief Create a new test suite.
 *
 * Sets up all nessesary data for a new test suite.
 *
 * Typical usage:
 * \code
 *   struct test_suite *suite = NULL;
 *   create_test_suite(&suite);
 * \endcode
 *
 * \param suite pointer to where the suite should be created.
 *
 * \returns E_SUCCESS on success, E_FAIL if memory allocation fails.
 */
int create_test_suite(struct test_suite **suite) {
  const int initial_size = 128;
  struct test_suite *tmp =
      (struct test_suite *)malloc(sizeof(struct test_suite));
  if (tmp == NULL) {
    return E_FAIL;
  }
  tmp->capacity = initial_size;
  tmp->n = 0;
  void (**tmpf)(struct test *t) = (void (**)(struct test *))malloc(
      tmp->capacity * sizeof(void (*)(struct test *)));
  if (tmpf == NULL) {
    free(tmp);
    return E_FAIL;
  }
  tmp->tests = tmpf;
  *suite = tmp;

  return E_SUCCESS;
}

/**
 * \brief Setup test data.
 *
 * Creates all field to be written by the user test case functions.
 * The test object will hold all state regarding the current test case.
 *
 * Used internally by suite_run_all_tests(...)
 *
 * \param t pointer to where test data should be created.
 *
 * \returns E_SUCCESS on success, E_FAIL memory allocation fails.
 */
int create_test(struct test **t) {
  const size_t buffer_size = 4096;
  const size_t msg_size = 1024;
  struct test *tmp_test = (struct test *)malloc(sizeof(struct test));
  if (tmp_test == NULL) {
    return E_FAIL;
  }
  int err;
  tmp_test->status = UNDEFINED;
  tmp_test->msg_size = msg_size;

  tmp_test->error_msg = (char *)malloc(msg_size * sizeof(char));
  if (tmp_test->error_msg == NULL) {
    goto emsg_error;
  }
  tmp_test->test_msg = (char *)malloc(msg_size * sizeof(char));
  if (tmp_test->test_msg == NULL) {
    goto tmsg_error;
  }

  err = create_buffer(&(tmp_test->log), buffer_size);
  if (err != E_SUCCESS) {
    goto buffer_error;
  }
  *t = tmp_test;
  return E_SUCCESS;
buffer_error:
  free(tmp_test->test_msg);
tmsg_error:
  free(tmp_test->error_msg);
emsg_error:
  free(tmp_test);
  return E_FAIL;
}

/**
 * \brief Reset test data.
 *
 * Clears all messages, buffers and test state from test data.
 *
 * Used internally by suite_run_all_tests(...)
 *
 * \param t test data to clear
 *
 * \returns E_SUCCESS
 */
int reset_test(struct test *t) {
  reset_buffer(t->log);
  t->status = UNDEFINED;
  memset(t->error_msg, 0, t->msg_size);
  memset(t->test_msg, 0, t->msg_size);
  return E_SUCCESS;
}

/**
 * \brief Creates expandable buffer for storing log statements.
 *
 *
 *
 * Used internally by create_test(...)
 *
 * \param log pointer to where the buffer should be created
 * \param cap initial capacity of the buffer
 *
 * \returns E_SUCCESS, E_FAIL if memory allocation fails.
 */
int create_buffer(struct log_buffer **log, size_t cap) {
  struct log_buffer *tmp =
      (struct log_buffer *)malloc(sizeof(struct log_buffer));
  if (tmp == NULL) {
    return E_FAIL;
  }
  tmp->buffer = (char *)malloc(cap * sizeof(char));
  if (tmp->buffer == NULL) {
    free(tmp);
    return E_FAIL;
  }
  tmp->capacity = cap;
  tmp->free_space = cap;
  tmp->start = tmp->buffer;

  *log = tmp;
  return E_SUCCESS;
}

/**
 * \brief Clears the buffer of all data.
 *
 *
 *
 * Used internally by reset_test(...)
 *
 * \param log the buffer to clear
 *
 * \returns E_SUCCESS
 */
int reset_buffer(struct log_buffer *log) {
  log->buffer = log->start;
  log->free_space = log->capacity;
  memset(log->buffer, 0, log->capacity);
  return E_SUCCESS;
}

/**
 * \brief Cleanup for the buffer.
 *
 * Frees all allocated memory for the buffer.
 *
 * Used internally by free_test(...)
 *
 * \param log buffer to free
 *
 * \returns E_SUCCESS
 */
int free_buffer(struct log_buffer *log) {
  free(log->start);
  log->buffer = NULL;
  free(log);
  return E_SUCCESS;
}

/**
 * \brief Cleanup for test data.
 *
 * Frees all allocated memory for the test data.
 *
 * Used internally by suite_run_all_tests(...)
 *
 * \param t test data to free
 *
 * \returns E_SUCCESS
 */
int free_test(struct test *t) {
  free_buffer(t->log);
  free(t->error_msg);
  free(t->test_msg);
  free(t);
  return E_SUCCESS;
}

/**
 * \brief Cleanup for test suite.
 *
 * Frees all memory for the given test suite.
 *
 * Typical usage:
 * \code
 *   struct test_suite *suite = NULL;
 *   create_test_suite(&suite);
 *     ... testing ...
 *   free_suite(suite);
 * \endcode
 *
 * \param suite suite to free
 *
 * \returns E_SUCCESS
 */
int free_suite(struct test_suite *suite) {
  free(suite->tests);
  free(suite);
  return E_SUCCESS;
}
