#include "candidates.h"
#include "errors.h"
#include "coverage.h"

int report_valid_candiates(struct extended_candidate_list *ec_list);
int create_candidate_report(struct extended_candidate *ecand,
                            const char *output_path);
int create_structure_image(struct extended_candidate *ecand,
                           const char *output_path);
int create_file_path(char **file_path, const char *path, const char *filename);