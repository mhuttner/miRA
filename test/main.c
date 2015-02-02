#include <stdio.h>
#include "testerino.h"
#include "test_parse_sam.h"
#include "test_cluster.h"
#include "test_bed_file_io.h"
#include "test_fasta.h"
#include "test_util.h"
#include "test_vfold.h"

int main(int argc, char const *argv[]) {
  struct test_suite *s = NULL;
  create_test_suite(&s);
  suite_add_test(s, test_valid_line);
  suite_add_test(s, test_header_line);
  suite_add_test(s, test_invalid_line);
  suite_add_test(s, test_multiple_consecutive_tabs);
  suite_add_test(s, test_sort_clusters);
  suite_add_test(s, test_merge_clusters);
  suite_add_test(s, test_filter_clusters);
  suite_add_test(s, test_merge_extended_clusters);
  suite_add_test(s, test_valid_bed_line);
  suite_add_test(s, test_invalid_start_bed_line);
  suite_add_test(s, test_invalid_id_bed_line);
  suite_add_test(s, test_strip_newlines);
  suite_add_test(s, test_read_fasta_file);
  suite_add_test(s, test_reverse_complement);
  suite_add_test(s, test_mean);
  suite_add_test(s, test_sd);
  suite_add_test(s, test_pvalue);
  suite_add_test(s, test_config_parsing);
  // suite_add_test(s, test_folding);
  suite_run_all_tests(s);
  free_suite(s);
  return 0;
}