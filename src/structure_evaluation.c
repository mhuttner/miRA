#include "vfold.h"
#include "errors.h"
#include "structure_evaluation.h"

int evaluate_structure(struct structure_info *si) {
  if (si == NULL) {
    return E_NO_STRUCTURE;
  }
  determine_external_loop_count(si);
  determine_paired_fraction(si);
  determine_max_stem_length(&si->stem_start, &si->stem_end, si, 0);
  determine_max_stem_length(&si->stem_start_with_mismatch,
                            &si->stem_end_with_mismatch, si, 1);

  return E_SUCCESS;
}

int determine_external_loop_count(struct structure_info *si) {
  return get_loop_count_for_substructure(&(si->external_loop_count),
                                         si->structure_string, 0, si->n);
}
int get_loop_count_for_substructure(int *result, char *structure_string,
                                    size_t start, size_t end) {
  int loop_count = 0;
  int status = 0;
  for (int i = start; i < end; i++) {
    if (structure_string[i] == '(') {
      status = 1;
    }
    if (status == 1 && structure_string[i] == ')') {
      status = 0;
      loop_count++;
    }
  }
  *result = loop_count;
  return E_SUCCESS;
}

int determine_paired_fraction(struct structure_info *si) {
  return get_paired_fraction_for_substructure(&(si->paired_fraction),
                                              si->structure_string, 0, si->n);
}
int get_paired_fraction_for_substructure(double *result, char *structure_string,
                                         size_t start, size_t end) {
  if (end <= start) {
    *result = 0.0;
    return E_INVALID_RANGE;
  }
  size_t dot_count = 0;
  for (size_t i = start; i < end; i++) {
    if (structure_string[i] == '.') {
      dot_count++;
    }
  }
  *result = 1.0 - (double)dot_count / (end - start);
  return E_SUCCESS;
}

int determine_max_stem_length(int *stem_start, int *stem_end,
                              struct structure_info *si,
                              int max_mismatch_count) {
  int status = 0;
  int stem_length = 0;
  int mismatches = 0;
  int start = 0;

  int longest_stem = 0;
  int longest_start = 0;
  int longest_end = 0;

  for (int i = 0; i < si->n; i++) {
    if (status == '(' || status == ')') {
      if (si->structure_string[i] == status) {
        stem_length++;
        continue;
      }
      if (si->structure_string[i + 1] == status &&
          (mismatches < max_mismatch_count)) {
        mismatches++;
        stem_length++;
        continue;
      }
      if (stem_length > longest_stem) {
        longest_stem = stem_length;
        longest_end = i - 1;
        longest_start = start;
      }
      stem_length = 0;
      mismatches = 0;
    }
    status = si->structure_string[i];
    start = i;
  }
  *stem_start = longest_start;
  *stem_end = longest_end;
  return E_SUCCESS;
}
