#include "vfold.h"

#ifndef STRUCTURE_EVALUATION_H
#define STRUCTURE_EVALUATION_H

int evaluate_structure(struct foldable_sequence *fs);
int determine_external_loop_count(struct structure_info *si);
int determine_paired_fraction(struct structure_info *si);
int determine_max_stem_length(int *stem_start, int *stem_end,
                              struct structure_info *si,
                              int max_mismatch_count);

#endif