#include "candidates.h"
#include "errors.h"
#include "coverage.h"

int report_valid_candiates(struct extended_candidate_list *ec_list,
                           struct chrom_coverage **coverage_table);
int create_candidate_report(struct extended_candidate *ecand,
                            struct chrom_coverage *chrom_cov,
                            const char *output_path);
int create_coverage_plot(char **result_file, struct extended_candidate *ecand,
                         struct chrom_coverage *chrom_cov,
                         const char *output_path);
int create_structure_image(char **result_file, struct extended_candidate *ecand,
                           const char *output_path);
int create_coverage_image(char **result_file, struct extended_candidate *ecand,
                          struct chrom_coverage *chrom_cov,
                          const char *output_path);
int create_latex_template(struct extended_candidate *ecand,
                          struct chrom_coverage *chrom_cov,
                          const char *cov_plot_file, const char *structure_file,
                          const char *coverage_file, const char *output_path);
int write_unique_reads_to_tex(FILE *fp, struct micro_rna_candidate *cand,
                              struct candidate_subsequence *subseq,
                              struct unique_read_list *reads);
int map_coverage_to_color_index(u32 *result, u32 coverage);
int create_file_path(char **file_path, const char *path, const char *filename);
const char *path_to_filename(const char *path);