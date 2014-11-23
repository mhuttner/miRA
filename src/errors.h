#ifndef ERRORS_H
#define ERRORS_H

enum _error {
  E_SUCCESS = 0,
  E_FILE_NOT_FOUND = -1,
  E_NO_FILE_SPECIFIED = -2,
  E_INVALID_SAM_LINE = -3,
  E_SAM_HEADER_LINE = -4,
  E_INVALID_ARGUMENT = -5,
  E_SAM_NON_SQ_HEADER = -6,
  E_CHROMOSOME_NOT_FOUND = -7,
  E_FILE_WRITING_FAILED = -8,
  E_INVALID_BED_LINE = -9,
  E_UNKNOWN_FILE_IO_ERROR = -10,
  E_INVALID_FASTA_FILE = -11,

  E_MALLOC_FAIL = -30,
  E_REALLOC_FAIL = -31,
  E_UNKNOWN = -99
};

struct _errordesc {
  int code;
  const char *message;
};

int print_error(int err);

#endif