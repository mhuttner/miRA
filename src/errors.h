#ifndef ERRORS_H
#define ERRORS_H 

enum _error {
    E_SUCCESS = 0,
    E_FILE_NOT_FOUND = -1,
    E_INVALID_FILE_FORMAT = -2,
    E_INVALID_SAM_LINE = -3,
    E_SAM_HEADER_LINE = -4
};

struct _errordesc {
    int code;
    const char *message;
};

#endif