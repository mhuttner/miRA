#include "errors.h"

struct _errordesc errordesc[] = {
    { E_SUCCESS, "No Error"},
    { E_FILE_NOT_FOUND, "The File was not found"},
    { E_INVALID_FILE_FORMAT, "The specified File was invalid"},
    { E_INVALID_SAM_LINE, "The Line of the SAM File is invalid"},
    { E_SAM_HEADER_LINE, "The Line is a Header Line and was not parsed"}
};