#include "../config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "help.h"
#include "cluster.h"
#include "vfold.h"
#include "errors.h"
#include "coverage.h"
#include "full.h"
#include "batch.h"

int main(int argc, char **argv) {
  /* List of all available operations */
  const char *operations[] = {"help",     "cluster", "fold",
                              "coverage", "full",    "batch"};
  const int num_operations = 6;

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
  case 2: /* fold */
    return vfold(argc - 1, argv + 1);
  case 3: /* coverage */
    return coverage(argc - 1, argv + 1);
  case 4: /* full */
    return full(argc - 1, argv + 1);
  case 5: /*batch */
    return batch(argc - 1, argv + 1);
  default:
    break;
  }

  return E_SUCCESS;
}
