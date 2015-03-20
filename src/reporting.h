#include "candidates.h"
#include "errors.h"
#include "coverage.h"

int report_valid_candiates(struct extended_candidate_list *ec_list,
                           struct chrom_coverage **coverage_table,
                           const char *executable_path, const char *output_path,
                           struct configuration_params *config);
int create_output_directory_structure(char **cov_plot_ouput_path,
                                      char **structure_output_path,
                                      char **coverage_output_path,
                                      char **report_output_path,
                                      const char *output_path);
int create_directory_if_ne(const char *path);
int create_candidate_report(
    struct extended_candidate *ecand, struct chrom_coverage *chrom_cov,
    const char *executable_path, const char *cov_plot_ouput_path,
    const char *structure_output_path, const char *coverage_output_path,
    const char *report_output_path, struct configuration_params *config);
int create_coverage_plot(char **result_file, struct extended_candidate *ecand,
                         struct chrom_coverage *chrom_cov,
                         const char *output_path);
int create_structure_image(char **result_file, struct extended_candidate *ecand,
                           const char *executable_path,
                           const char *output_path);
int create_coverage_image(char **result_file, struct extended_candidate *ecand,
                          const char *executable_path,
                          struct chrom_coverage *chrom_cov,
                          const char *output_path);
int create_latex_template(char **tex_file, struct extended_candidate *ecand,
                          struct chrom_coverage *chrom_cov,
                          const char *cov_plot_file, const char *structure_file,
                          const char *coverage_file, const char *output_path);

int write_unique_reads_to_tex(FILE *fp, struct micro_rna_candidate *cand,
                              struct candidate_subsequence *subseq,
                              struct unique_read_list *reads);
int compile_tex_file(const char *tex_file_path, const char *output_path);
int map_coverage_to_color_index(u32 *result, u32 coverage);
int write_bed_lines(FILE *fp, struct extended_candidate *ecand);
int write_gtf_line(FILE *fp, struct extended_candidate *ecand);
int write_json_entry(FILE *fp, struct extended_candidate *ecand);
int inititalize_html_report(FILE *fp);
int write_html_table_row(FILE *fp, struct extended_candidate *ecand);
int finalize_html_report(FILE *fp);
int cleanup_auxiliary_files(char *cov_plot_file, char *structure_file,
                            char *coverage_file, char *tex_file,
                            struct configuration_params *config);
