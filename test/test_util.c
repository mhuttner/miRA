#include "testerino.h"
#include "../src/util.h"
#include <math.h>

void test_mean(struct test *t) {
  t_set_msg(t, "Testing the mean function...");
  double test_list[] = {3.0, 5.0, 3.0, 6.0, 3.0, 2.0, 4.0, 6.0, 7.0, 4.0, 3.0};
  double result = mean(test_list, 11);
  t_assert_msg(t, fabs(result - 4.181818) < 10e-4, "mean function incorrect");
}

void test_sd(struct test *t) {
  t_set_msg(t, "Testing the standard deviation function...");
  double test_list[] = {3.0, 5.0, 3.0, 6.0, 3.0, 2.0, 4.0, 6.0, 7.0, 4.0, 3.0};
  double test_mean = mean(test_list, 11);
  double result = sd(test_list, 11, test_mean);
  t_assert_msg(t, fabs(result - 1.601136) < 10e-4, "sd function incorrect");
}

void test_pvalue(struct test *t) {
  t_set_msg(t, "Testing the pvalue function...");
  double test_list[] = {3.0, 5.0, 3.0, 6.0, 3.0, 2.0, 4.0, 6.0, 7.0, 4.0, 3.0};
  double test_mean = mean(test_list, 11);
  double test_sd = sd(test_list, 11, test_mean);
  double result = pvalue(test_mean, test_sd, 1.0);
  t_log(t, "pvalue: %lf", result);
  t_assert_msg(t, fabs(result - 0.023448) < 10e-4, "pvalue function incorrect");
}

void test_config_parsing(struct test *t) {
  t_set_msg(t, "Testing parsing the config file...");
  const char *file = "test/data/test.config";
  struct configuration_params *config = NULL;
  initialize_configuration(&config, (char *)file);

  t_log(t, "log_level %d\n", config->log_level);
  t_log(t, "max_precursor_length %d\n", config->max_precursor_length);
  t_log(t, "min_precursor_length %d\n", config->min_precursor_length);
  t_log(t, "max_mfe_per_nt %lf\n", config->max_mfe_per_nt);
  t_log(t, "max_hairpin_count %d\n", config->max_hairpin_count);
  t_log(t, "min_double_strand_length %d\n", config->min_double_strand_length);
  t_log(t, "permutation_count %d\n", config->permutation_count);
  t_log(t, "max_pvalue %lf\n", config->max_pvalue);
  t_log(t, "min_coverage %lf\n", config->min_coverage);
  t_log(t, "min_paired_fraction %d\n", config->min_paired_fraction);
  t_log(t, "min_duplex_length %d\n", config->min_duplex_length);
  t_log(t, "max_duplex_length %d\n", config->max_duplex_length);
  t_log(t, "allow_three_mismatches %d\n", config->allow_three_mismatches);
  t_log(t, "allow_two_terminal_mismatches %d\n",
        config->allow_two_terminal_mismatches);

  t_assert_msg(t, config->log_level == 1, "log_level wrong");
  t_assert_msg(t, config->cluster_gap_size == -37, "cluster_gap_size wrong");
  t_assert_msg(t, config->cluster_min_reads == -38, "cluster_min_reads wrong");
  t_assert_msg(t, config->cluster_flank_size == -39,
               "cluster_flank_size wrong");
  t_assert_msg(t, config->cluster_max_length == -40,
               "cluster_max_length wrong");

  t_assert_msg(t, config->max_precursor_length == -41,
               "max_precursor_length wrong");
  t_assert_msg(t, config->min_precursor_length == 50,
               "min_precursor_length wrong");
  t_assert_msg(t, config->max_mfe_per_nt == -0.24, "max_mfe_per_nt wrong");
  t_assert_msg(t, config->max_hairpin_count == -42, "max_hairpin_count wrong");
  t_assert_msg(t, config->min_double_strand_length == -43,
               "min_double_strand_length wrong");
  t_assert_msg(t, config->permutation_count == -44, "permutation_count wrong");
  t_assert_msg(t, config->max_pvalue == -0.42, "max_pvalue wrong");
  t_assert_msg(t, config->min_coverage == -0.43, "min_coverage wrong");
  t_assert_msg(t, config->min_paired_fraction == -0.44,
               "min_paired_fraction wrong");
  t_assert_msg(t, config->min_duplex_length == -45,
               "min_duplex_length wrong");
  t_assert_msg(t, config->max_duplex_length == -46,
               "max_duplex_length wrong");
  t_assert_msg(t, config->allow_three_mismatches == -47,
               "allow_three_mismatches wrong");
  t_assert_msg(t, config->allow_two_terminal_mismatches == -48,
               "allow_two_terminal_mismatches wrongd\n");
}
