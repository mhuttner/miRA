#include <stdio.h>
#include <stdarg.h>
#include "errors.h"

struct _errordesc errordesc[] = {
    {E_SUCCESS, "No Error"},
    {E_FILE_NOT_FOUND, "The File was not found"},
    {E_NO_FILE_SPECIFIED, "No file was given"},
    {E_INVALID_SAM_LINE, "The Line of the SAM File is invalid"},
    {E_SAM_HEADER_LINE, "The Line is a Header Line"},
    {E_INVALID_ARGUMENT, "Invalid commad line argument"},
    {E_SAM_NON_SQ_HEADER, "The Line is a Header Line that was ignored"},
    {E_CHROMOSOME_NOT_FOUND, "The Chromosome was not found as a @SQ entry"},
    {E_FILE_WRITING_FAILED, "Failed to write the output File"},
    {E_INVALID_BED_LINE, "The line of the BED file is invalid"},
    {E_UNKNOWN_FILE_IO_ERROR, "An unknown error occurred reading a file"},
    {E_INVALID_FASTA_FILE, "Given fasta file is invalid"},
    {E_NEEDED_SEQUENCE_NOT_FOUND,
     "A fasta sequence needed for folding was not found"},
    {E_CREATING_PIPE_FAILED, "A file descriptor pipe could not be created"},
    {E_INVALID_FASTA_SEQUENCE_LENGTH, "Cluster out of genome range"},
    {E_NO_OPTIMAL_STRUCTURE_FOUND,
     "No secondary Structure containing the core region found"},
    {E_MALLOC_FAIL, "malloc failed, check available memory"},
    {E_REALLOC_FAIL, "realloc failed, check available memory"},
    {E_NO_CLUSTERS_LEFT, "No clusters left to work with."},
    {E_UNKNOWN, "An unknown error occured"}};

int print_error(int err) {
  for (int i = 0; errordesc[i].code != E_UNKNOWN; i++) {
    if (err == errordesc[i].code) {
      fprintf(stderr, "%s\n", errordesc[i].message);
      break;
    }
  }
  return E_SUCCESS;
}
