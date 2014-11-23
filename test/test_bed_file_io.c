#include "testerino.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../src/bed.h"
#include "../src/cluster.h"
#include "../src/errors.h"

void test_valid_bed_line(struct test *t) {
  t_set_msg(t, "Testing reading a valid Bed line...");
  char sample_line[1000] =
      "scaffold_1\t186203\t186647\tCluster_0\t0\t+\t186403\t186447\t0\t"
      "90\n";
  struct cluster *c = NULL;
  int result = parse_bed_line(&c, sample_line);
  t_assert_msg(t, result == E_SUCCESS, "parsing failed");
  if (result != E_SUCCESS)
    return;
  t_assert_msg(t, strcmp(c->chrom, "scaffold_1") == 0, "Chromsome read wrong");
  t_assert_msg(t, c->end == 186647, "Cluster end wrong");
  t_assert_msg(t, c->id == 0, "Cluster Id wrong");
  t_assert_msg(t, c->strand == '+', "Strand wrong");
  free_cluster(c);
}

void test_invalid_start_bed_line(struct test *t) {
  t_set_msg(t, "Testing reading a BED line with invalid start...");
  char sample_line[1000] =
      "scaffold_1\t186v203\t186647\tCluster_0\t0\t+\t186403\t186447\t0\t"
      "90\n";
  struct cluster *c = NULL;
  int result = parse_bed_line(&c, sample_line);
  t_assert_msg(t, result == E_INVALID_BED_LINE, "Invalid line got parsed");
}

void test_invalid_id_bed_line(struct test *t) {
  t_set_msg(t, "Testing reading a BED line with invalid id...");
  char sample_line[1000] =
      "scaffold_1\t186203\t186647\tClusterino_0\t0\td\t186403\t186447\t0\t"
      "90\n";
  struct cluster *c = NULL;
  int result = parse_bed_line(&c, sample_line);
  t_assert_msg(t, result == E_INVALID_BED_LINE, "Invalid line got parsed");
}
