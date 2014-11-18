#include <stdio.h>
#include "vfold.h"
#include "errors.h"

static int print_help();

int vfold(int argc, char *argv[]) {
  print_help();
  return E_SUCCESS;
}

static int print_help() {
  printf("Description:\n"
         "    fold tries to fold rna sequences and calculates secondary\n"
         "    structure information \n"
         "Usage: miRA fold <input BED file> <input FASTA file>\n");
  return E_SUCCESS;
}