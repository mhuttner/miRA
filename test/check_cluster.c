#include <stdlib.h>
#include <check.h>
#include "../src/parse_sam.h"
#include "../src/cluster.h"
#include "../src/errors.h"

START_TEST(test_create_clusters) {
  struct sam_file sam;
  sam.n = 4;
  sam.capacity = 4;
  struct sam_entry test_entries[] = {
      {"abc", 0, "chrom_2", 200, 0, "", "", 0, 0, "ABCDEFGHIJKLMN", ""},
      {"abc", 0, "chrom_2", 100, 0, "", "", 0, 0, "ABCDEFGHIJKLMN", ""},
      {"abc", 0x10, "chrom_2", 150, 0, "", "", 0, 0, "ABCDEFGHIJKLMN", ""},
      {"abc", 0, "chrom_1", 10, 0, "", "", 0, 0, "ABCDEFGHIJKLMN", ""}};
  sam.entries = test_entries;

  struct cluster_list *list = NULL;
  int err = create_clusters(&list, &sam);
  ck_assert_int_eq(err, E_SUCCESS);
  free_clusters(list);
}
END_TEST

START_TEST(test_sort_clusters) {
  struct cluster_list list;
  list.n = 4;
  list.capacity = 4;
  struct cluster test_clusters[] = {{0, '-', "chrom_2", 10l, 20l, 1},
                                    {1, '+', "chrom_2", 30l, 40l, 1},
                                    {2, '+', "chrom_2", 15l, 25l, 1},
                                    {3, '+', "chrom_1", 10l, 20l, 1}};
  list.clusters = test_clusters;
  sort_clusters(&list, compare_strand_chrom_start);
  ck_assert_int_eq(list.clusters[0].id, 3);
  ck_assert_int_eq(list.clusters[1].id, 2);
  ck_assert_int_eq(list.clusters[2].id, 1);
  ck_assert_int_eq(list.clusters[3].id, 0);
}
END_TEST

Suite *hello_suite(void) {
  Suite *s = suite_create("Cluster");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_create_clusters);
  tcase_add_test(tc_core, test_sort_clusters);
  suite_add_tcase(s, tc_core);

  return s;
}

int main(void) {
  int number_failed;
  Suite *s = hello_suite();
  SRunner *sr = srunner_create(s);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}