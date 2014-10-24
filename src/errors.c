#include <stdio.h>
#include "errors.h"

struct _errordesc errordesc[] = {
    { E_SUCCESS, "No Error"},
    { E_FILE_NOT_FOUND, "The File was not found"},
    { E_NO_FILE_SPECIFIED, "No file was given"},
    { E_INVALID_SAM_LINE, "The Line of the SAM File is invalid"},
    { E_SAM_HEADER_LINE, "The Line is a Header Line and was not parsed"},
    { E_MALLOC_FAIL , "malloc failed, check available memory"},
    { E_REALLOC_FAIL, "realloc failed, check available memory"},
    { E_UNKNOWN, "An unknown error occured"}
};

int print_error(int err){
    for(int i=0;errordesc[i].code != E_UNKNOWN;i++){
        if(err == errordesc[i].code){
            fprintf(stderr, "%s\n",errordesc[i].message);
            break;
        }
    }
    return E_SUCCESS;
}