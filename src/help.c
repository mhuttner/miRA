#include <stdio.h>
#include "../config.h"
#include "help.h"
#include "errors.h"

int help(int argc, char **argv) {
  printf(
      "Package: %s\n"
      "Version: %s\n"
      "Description: \n"
      "    Identification of novel miRNAs software implementation \n"
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
      "    coverage   coverage based verification and reporting\n"
      "    full       run the full miRA algorithm cluster, fold and \n"
      "               coverage in sequence\n"
      "    help       show this help message\n"
      "\n"
      "Example Usage:\n"
      "    You can test running miRA on provided example date in the folder \n"
      "    \"example\" with the command:\n\n"
      "    >./miRA full -c example/sample_configuration.config  \n"
      "    example/sample_reads.sam example/sample_sequence.fasta \n"
      "    example/sample_output/ \n"
      "\n",
      PACKAGE, PACKAGE_VERSION);
  return E_SUCCESS;
}