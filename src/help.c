#include <stdio.h>
#include "../config.h"
#include "help.h"
#include "errors.h"

int help(int argc, char **argv) {
  printf("Package: %s\n"
         "Version: %s\n"
         "Description: \n"
         "    Identification of novel miRNAs software implementation"
         "Authors:\n"
         "    Maurits Evers (maurits.evers@ur.de)\n"
         "    Michael Huttner (michael@mhuttner.com) \n"
         "Maintainers:\n"
         "    Michael Huttner (michael@mhuttner.com) \n"
         "\n"
         "Usage: miRA <command> [<args>]\n"
         "\n"
         "The most commonly used commands are:\n"
         "    cluster    generate list of main expression contigs based on \n"
         "               alignment data.\n"
         "    fold       calculate secondary rna structure\n"
         "    help       show this help message\n"
         "\n",
         PACKAGE, PACKAGE_VERSION);
  return E_SUCCESS;
}