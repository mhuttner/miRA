#include "reporting.h"
#include "coverage.h"
#include "candidates.h"
#include "errors.h"
#include "defs.h"
#include <stdio.h>
#include "../config.h"

int report_valid_candiates(struct extended_candidate_list *ec_list) {
  char output_path[] = "./data/output";
  struct extended_candidate *ecand = NULL;
  for (size_t i = 0; i < ec_list->n; i++) {
    ecand = ec_list->candidates[i];
    if (ecand->is_valid != 1) {
      continue;
    }
    create_candidate_report(ecand, output_path);
  }

  return E_SUCCESS;
}

int create_candidate_report(struct extended_candidate *ecand,
                            const char *output_path) {

  int err = 1;
#ifdef HAVE_JAVA
  err = create_structure_image(ecand, output_path);
#endif /* HAVE_JAVA */
  if (err != E_SUCCESS) {
    return E_CREATING_REPORT_FAILED;
  }

  return E_SUCCESS;
}

int create_structure_image(struct extended_candidate *ecand,
                           const char *output_path) {
  const int SYS_CALL_MAX_LENGHT = 4096;
  const char varna_path[] = "./VARNAv3-91.jar";
  const char varna_class_name[] = "fr.orsay.lri.varna.applications.VARNAcmd";
  const char varna_algorith[] = "naview";
  const char varna_auto_interior_loops[] = "True";
  const char varna_auto_terminal_loops[] = "True";
  const int varna_resolution = 10;
  const char hex_color_red[] = "#FF0000";
  const char hex_color_blue[] = "#0000FF";

  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence *mature_mirna = ecand->mature_micro_rna;
  struct candidate_subsequence *star_mirna = ecand->star_micro_rna;

  char file_name[265];
  printf("1\n");
  sprintf(file_name, "Cluster_%lld_structure_image.eps", cand->id);
  char *file_path = NULL;
  printf("2\n");
  create_file_path(&file_path, output_path, file_name);

  u32 m_start = mature_mirna->start + 1;
  u32 m_end = mature_mirna->end;
  u32 s_start = star_mirna->start + 1;
  u32 s_end = star_mirna->end;
  char *seq = cand->sequence;
  char *structure = cand->structure;

  char java_cmd[20] = "java";
#ifdef JAVA_VM_COMMAND
  sprintf(java_cmd, "%s", JAVA_VM_COMMAND);
#endif

  char java_system_call[4096];
  snprintf(java_system_call, SYS_CALL_MAX_LENGHT,
           "java -cp %s %s -algorithm %s -autoInteriorLoops %s "
           "-autoTerminalLoops %s -resolution %d -highlightRegion "
           "\"%d-%d:fill=%s;%d-%d:fill=%s\" -title \"Cluster %lld\" "
           "-sequenceDBN \"%s\" "
           "-structureDBN \"%s\" -o %s",
           varna_path, varna_class_name, varna_algorith,
           varna_auto_interior_loops, varna_auto_terminal_loops,
           varna_resolution, m_start, m_end, hex_color_red, s_start, s_end,
           hex_color_blue, cand->id, seq, structure, file_path);
  int err = system(java_system_call);
  if (err != 0) {
    return E_JAVA_SYSTEM_CALL_FAILED;
  }
  return E_SUCCESS;
}

int create_file_path(char **file_path, const char *path, const char *filename) {
  size_t path_n = strnlen(path, 1024);
  size_t file_n = strnlen(filename, 1024);
  if (path_n >= 1024 || file_n >= 1024) {
    return E_FILE_DESCRIPTOR_TOO_LONG;
  }
  char *file_path_tmp = (char *)malloc((path_n + file_n + 2) * sizeof(char));
  if (file_path_tmp == NULL) {
    return E_MALLOC_FAIL;
  }
  if (path[path_n - 1] == '/') {
    sprintf(file_path_tmp, "%s%s", path, filename);
  } else {
    sprintf(file_path_tmp, "%s/%s", path, filename);
  }
  *file_path = file_path_tmp;
  return E_SUCCESS;
}
