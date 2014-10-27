#include "../config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "help.h"
#include "cluster.h"
#include "errors.h"

int main(int argc, char **argv) {
  /* List of all available operations */
  const char *operations[] = {"help", "cluster"};
  const int num_operations = 2;

  int operation_type = 0;
  if (argc >= 2) {
    char *operation = argv[1];
    for (int i = 0; i < num_operations; i++) {
      if (strcasecmp(operation, operations[i]) == 0) {
        operation_type = i;
      }
    }
  }

  switch (operation_type) {
  case 0: /* help */
    return help(argc, argv);
  case 1: /* cluster */
    return cluster(argc - 1, argv + 1);
  default:
    break;
  }

  return E_SUCCESS;
}
