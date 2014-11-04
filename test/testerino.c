#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "testerino.h"

void t_set_msg(struct test *t, const char *msg) {
  snprintf(t->test_msg, t->msg_size, "%s", msg);
}

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
      log->free_space = 0;
    } else {
      log->buffer += n;
      log->free_space -= n;
    }
  }
  va_end(fmtargs);
}
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

int suite_add_test(struct test_suite *suite, void (*test)(struct test *)) {
  if (suite->n >= suite->capacity) {
    // realloc
    return E_FAIL;
  }
  suite->tests[suite->n] = test;
  suite->n++;
  return E_SUCCESS;
}

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
    t->status = UNDEFINED;
    reset_buffer(t->log);
  }
  free_test(t);
  return E_SUCCESS;
}
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

int create_test_suite(struct test_suite **suite) {
  const int initial_size = 128;
  struct test_suite *tmp =
      (struct test_suite *)malloc(sizeof(struct test_suite));
  if (tmp == NULL) {
    return E_FAIL;
  }
  tmp->capacity = initial_size;
  tmp->n = 0;
  tmp->failed = 0;
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

int reset_buffer(struct log_buffer *log) {
  log->buffer = log->start;
  log->free_space = log->capacity;
  return E_SUCCESS;
}

int free_buffer(struct log_buffer *log) {
  free(log->start);
  log->buffer = NULL;
  free(log);
  return E_SUCCESS;
}

int free_test(struct test *t) {
  free_buffer(t->log);
  free(t->error_msg);
  free(t->test_msg);
  free(t);
  return E_SUCCESS;
}

int free_suite(struct test_suite *suite) {
  free(suite->tests);
  free(suite);
  return E_SUCCESS;
}
