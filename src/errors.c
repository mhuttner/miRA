#include "errors.h"

struct _errordesc errordesc[] = {
    { E_SUCCESS, "No Error"},
    { E_FILE_NOT_FOUND, "The File was not found"},
    { E_NO_FILE_SPECIFIED, "No file was given"},
    { E_INVALID_SAM_LINE, "The Line of the SAM File is invalid"},
    { E_SAM_HEADER_LINE, "The Line is a Header Line and was not parsed"}
};