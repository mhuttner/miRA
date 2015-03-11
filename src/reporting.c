#include "reporting.h"
#include "coverage.h"
#include "candidates.h"
#include "errors.h"
#include "defs.h"
#include <stdio.h>
#include "../config.h"
#include "math.h"
#include "util.h"

int report_valid_candiates(struct extended_candidate_list *ec_list,
                           struct chrom_coverage **coverage_table,
                           char *output_path,
                           struct configuration_params *config) {
  char bed_result_filename[] = "final_candidates.bed";
  char *bed_file = NULL;
  create_file_path(&bed_file, output_path, bed_result_filename);
  FILE *bed_fp = fopen(bed_file, "w");
  free(bed_file);
  char json_result_filename[] = "final_candidates.json";
  char *json_file = NULL;
  create_file_path(&json_file, output_path, json_result_filename);
  FILE *json_fp = fopen(json_file, "w");
  free(json_file);
  fprintf(json_fp, "{\n");

  struct extended_candidate *ecand = NULL;
  struct chrom_coverage *chrom_cov = NULL;
  for (size_t i = 0; i < ec_list->n; i++) {
    ecand = ec_list->candidates[i];
    if (ecand->is_valid != 1) {
      continue;
    }
    HASH_FIND_STR(*coverage_table, ecand->cand->chrom, chrom_cov);
    if (chrom_cov == NULL) {
      return E_CHROMOSOME_NOT_FOUND;
    }
    create_candidate_report(ecand, chrom_cov, output_path, config);
    if (bed_fp != NULL) {
      write_bed_lines(bed_fp, ecand);
    }
    if (json_fp != NULL) {
      write_json_entry(json_fp, ecand);
      if (i < ec_list->n - 1) {
        fprintf(json_fp, ",\n");
      }
    }
  }
  if (bed_fp != NULL) {
    fclose(bed_fp);
  }
  if (json_fp != NULL) {
    fprintf(json_fp, "}");
    fclose(json_fp);
  }

  return E_SUCCESS;
}

int create_candidate_report(struct extended_candidate *ecand,
                            struct chrom_coverage *chrom_cov,
                            const char *output_path,
                            struct configuration_params *config) {
  char *cov_plot_file = NULL;
  char *structure_file = NULL;
  char *coverage_file = NULL;
  char *tex_file = NULL;
  int err;

#ifdef HAVE_GNUPLOT
  if (config->create_coverage_plots) {
    err = create_coverage_plot(&cov_plot_file, ecand, chrom_cov, output_path);
    if (err != E_SUCCESS) {
      cov_plot_file = NULL;
    }
  }
#endif /* HAVE_GNUPLOT */
#ifdef HAVE_JAVA
  if (config->create_structure_images) {
    err = create_structure_image(&structure_file, ecand, output_path);
    if (err != E_SUCCESS) {
      structure_file = NULL;
    }
  }
#endif /* HAVE_JAVA */

#ifdef HAVE_JAVA
  if (config->create_coverage_images) {
    err = create_coverage_image(&coverage_file, ecand, chrom_cov, output_path);
    if (err != E_SUCCESS) {
      coverage_file = NULL;
    }
  }
#endif /* HAVE_JAVA */
  err = create_latex_template(&tex_file, ecand, chrom_cov, cov_plot_file,
                              structure_file, coverage_file, output_path);
  if (err != E_SUCCESS) {
    cleanup_auxiliary_files(cov_plot_file, structure_file, coverage_file,
                            tex_file, config);
    return E_CREATING_REPORT_FAILED;
  }

#ifdef HAVE_LATEX
  err = compile_tex_file(tex_file, output_path);
  if (err != E_SUCCESS) {
    cleanup_auxiliary_files(cov_plot_file, structure_file, coverage_file,
                            tex_file, config);
    return E_CREATING_REPORT_FAILED;
  }
#endif /* HAVE_LATEX */
  cleanup_auxiliary_files(cov_plot_file, structure_file, coverage_file,
                          tex_file, config);

  return E_SUCCESS;
}
int create_coverage_plot(char **result_file, struct extended_candidate *ecand,
                         struct chrom_coverage *chrom_cov,
                         const char *output_path) {
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence *mature_mirna = ecand->mature_micro_rna;
  struct candidate_subsequence *star_mirna = ecand->star_micro_rna;

  char data_file_name[265];
  sprintf(data_file_name, "Cluster_%lld_coverage_data.dat", cand->id);
  char *data_file_path = NULL;
  create_file_path(&data_file_path, output_path, data_file_name);

  char file_name[265];
  sprintf(file_name, "Cluster_%lld_coverage_plot.png", cand->id);
  char *file_path = NULL;
  create_file_path(&file_path, output_path, file_name);

  size_t l = cand->end - cand->start;
  u64 x_start = cand->start;
  u64 x_end = cand->end;
  u32 x2_start = 0;
  u32 x2_end = (u32)l;

  u64 box1_start = mature_mirna->start + cand->start + 1;
  u64 box1_end = mature_mirna->end + cand->start;

  u64 box2_start = star_mirna->start + cand->start + 1;
  u64 box2_end = star_mirna->end + cand->start;

  u32 *plus_cov = chrom_cov->coverage_plus + cand->start;
  u32 *minus_cov = chrom_cov->coverage_minus + cand->start;

  FILE *fp = fopen(data_file_path, "w");
  if (fp == NULL) {
    return E_FILE_WRITING_FAILED;
  }
  for (size_t i = 0; i < l; i++) {
    size_t global_index = cand->start + i;
    fprintf(fp, "%ld\t%d\t%d\t%ld\t%c\n", global_index + 1, plus_cov[i],
            minus_cov[i], i + 1, cand->sequence[i]);
  }
  fclose(fp);

  u32 y_max = 0;
  for (size_t i = 0; i < l; i++) {
    if (plus_cov[i] > y_max) {
      y_max = plus_cov[i];
    }
    if (minus_cov[i] > y_max) {
      y_max = minus_cov[i];
    }
  }
  y_max = (u32)(y_max * 1.5 * 1.5);

  char gnuplot_system_call[4096];
  sprintf(gnuplot_system_call,
          "gnuplot -p -e \""
          "set terminal png font 'Verdana,10';\n"
          "set output '%s';\n"
          " set xtics nomirror;\n"
          " set x2tics;\n"
          " set format x '%%s';\n"
          " set autoscale xfix;\n"
          " set autoscale x2fix;\n"
          " set style rect fc rgb '#084081' fs transparent solid 0.1 "
          "noborder;\n"
          " set logscale y;\n"
          " set obj rect from %lld, graph 0 to %lld, graph 1;\n"
          " set obj rect from %lld, graph 0 to %lld, graph 1;\n"
          " set yrange [0.1:%d];\n"
          " set xrange [%lld:%lld];\n"
          " set x2range [%d:%d];\n"
          " set xlabel 'Genome position (1-based)';\n"
          " set x2label 'Sequence position (1-based)';\n"
          " set ylabel 'Coverage (per nucleotide)';\n"
          " set key left top;\n"
          " plot '%s' using 1:2 with histeps lw 2 linecolor rgb 'red' title "
          "'plus','' using 4:2 with histeps lw 2 linecolor rgb 'red' axes x2y1 "
          "notitle, '' using 1:3 with histeps lw 2 linecolor rgb 'blue' title "
          "'minus', '' using 4:3 with histeps lw 2 linecolor rgb 'blue' "
          "notitle\"",
          file_path, box1_start, box1_end, box2_start, box2_end, y_max, x_start,
          x_end, x2_start, x2_end, data_file_path);
  int err = system(gnuplot_system_call);
  free(data_file_path);
  if (err != E_SUCCESS) {
    free(file_path);
    return E_GNUPLOT_SYSTEM_CALL_FAILED;
  }
  *result_file = file_path;
  return E_SUCCESS;
}

int create_structure_image(char **result_file, struct extended_candidate *ecand,
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
  sprintf(file_name, "Cluster_%lld_structure_image.eps", cand->id);
  char *file_path = NULL;
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
    free(file_path);
    return E_JAVA_SYSTEM_CALL_FAILED;
  }
  *result_file = file_path;
  return E_SUCCESS;
}

int create_coverage_image(char **result_file, struct extended_candidate *ecand,
                          struct chrom_coverage *chrom_cov,
                          const char *output_path) {
  const char varna_path[] = "./VARNAv3-91.jar";
  const char varna_class_name[] = "fr.orsay.lri.varna.applications.VARNAcmd";
  const char varna_algorith[] = "naview";
  const char varna_auto_interior_loops[] = "True";
  const char varna_auto_terminal_loops[] = "True";
  const int varna_resolution = 10;
  const char *coverage_colors[] = {"#EVEVFF", "#C2C2FF", "#ADADFF",
                                   "#7070FF", "#3333FF", "#0000F5",
                                   "#0000B8", "#00007A", "#00003D"};
  const char hex_color_white[] = "#FFFFFF";
  struct micro_rna_candidate *cand = ecand->cand;

  char *seq = cand->sequence;
  char *structure = cand->structure;

  u32 *cov_list = chrom_cov->coverage_plus;
  if (cand->strand == '-') {
    cov_list = chrom_cov->coverage_minus;
  }

  size_t segment_n = 50;
  size_t cand_n = cand->end - cand->start;
  u32 color_index = 0;

  char *highlight_string = (char *)malloc(segment_n * cand_n * sizeof(char));
  char *write_point = highlight_string;
  for (size_t i = 0; i < cand_n; i++) {
    size_t global_index = cand->start + i;
    map_coverage_to_color_index(&color_index, cov_list[global_index]);
    write_point +=
        sprintf(write_point, "%ld-%ld:fill=%s,outline=%s;", i + 1, i + 1,
                coverage_colors[color_index], hex_color_white);
  }
  char file_name[265];
  sprintf(file_name, "Cluster_%lld_coverage_image.eps", cand->id);
  char *file_path = NULL;
  create_file_path(&file_path, output_path, file_name);

  char java_cmd[20] = "java";
#ifdef JAVA_VM_COMMAND
  sprintf(java_cmd, "%s", JAVA_VM_COMMAND);
#endif
  size_t sys_call_n = 4096 + segment_n * cand_n;
  char *java_system_call = (char *)malloc(sys_call_n * sizeof(char));
  snprintf(java_system_call, sys_call_n,
           "java -cp %s %s -algorithm %s -autoInteriorLoops %s "
           "-autoTerminalLoops %s -resolution %d -highlightRegion "
           "\"%s\" -title \"Cluster %lld\" "
           "-sequenceDBN \"%s\" "
           "-structureDBN \"%s\" -o %s",
           varna_path, varna_class_name, varna_algorith,
           varna_auto_interior_loops, varna_auto_terminal_loops,
           varna_resolution, highlight_string, cand->id, seq, structure,
           file_path);
  int err = system(java_system_call);

  if (err != 0) {
    free(file_path);
    return E_JAVA_SYSTEM_CALL_FAILED;
  }
  *result_file = file_path;
  return E_SUCCESS;
}

int create_latex_template(char **tex_file, struct extended_candidate *ecand,
                          struct chrom_coverage *chrom_cov,
                          const char *cov_plot_file, const char *structure_file,
                          const char *coverage_file, const char *output_path) {

  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence *mature_mirna = ecand->mature_micro_rna;
  struct candidate_subsequence *star_mirna = ecand->star_micro_rna;
  struct unique_read_list *mature_reads = ecand->mature_reads;
  struct unique_read_list *star_reads = ecand->star_reads;

  char file_name[265];
  sprintf(file_name, "Cluster_%lld_report.tex", cand->id);
  char *file_path = NULL;
  create_file_path(&file_path, output_path, file_name);
  FILE *fp = fopen(file_path, "w");
  if (fp == NULL) {
    free(file_path);
    return E_FILE_WRITING_FAILED;
  }
  fprintf(fp, "\\documentclass[a4paper]{article}\n"
              "\\usepackage[margin=2cm,nohead,nofoot]{geometry}\n"
              "\\usepackage{graphicx}\n"
              "\\usepackage{color}\n"
              "\\parindent 0pt\n"
              "\\begin{document}\n"
              "\\section*{Summary}\n"
              "\\verb$Cluster %lld$\\\\\n"
              "\\verb$%s$\\\\\n"
              "Precursor miRNA : start = \\verb$%lld$, stop = \\verb$%lld$, "
              "length = \\verb$%lld$\\\\\n"
              "Mature miRNA : start = \\verb$%lld$ (\\verb$%d$), stop = "
              "\\verb$%lld$ (\\verb$%d$), length = \\verb$%d$, arm = "
              "\\verb$%d'$\\\\\n"
              "Star miRNA : start = \\verb$%lld$ (\\verb$%d$), stop = "
              "\\verb$%lld$ (\\verb$%d$), length = \\verb$%d$\\\\\n"
              "{[}Note: Positions are 1-based and inclusive{]}\\\\\n \\\\\n",
          cand->id, cand->chrom, cand->start, cand->end,
          cand->start - cand->end, cand->start + mature_mirna->start,
          mature_mirna->start, cand->start + mature_mirna->end,
          mature_mirna->end, mature_mirna->end - mature_mirna->start,
          (int)mature_mirna->arm, cand->start + star_mirna->start,
          star_mirna->start, cand->start + star_mirna->end, star_mirna->end,
          star_mirna->end - star_mirna->start);
  fprintf(fp, "Mature miRNA sequence: \\verb$");
  for (u32 i = mature_mirna->start; i < mature_mirna->end; i++) {
    fprintf(fp, "%c", cand->sequence[i]);
  }
  fprintf(fp, "$\\\\\n");
  fprintf(fp, "Star miRNA sequence: \\verb$");
  for (u32 i = star_mirna->start; i < star_mirna->end; i++) {
    fprintf(fp, "%c", cand->sequence[i]);
  }
  fprintf(fp, "$\\\\\n");

  fprintf(fp, "\\\\\n"
              "MFE = \\verb$%9.7e$ kcal/mol/nt\\\\\n"
              "p-value = \\verb$%9.7e$\\\\\n"
              "Null distribution: Mean = \\verb$%9.7e$ kcal/mol/nt, "
              "StandDev = \\verb$%9.7e$ kcal/mol/nt\\\\\n"
              "\\\\\n"
              "Number of terminal loops = \\verb$%i$\\\\\n"
              "Longest ds segment (0 MM): start = \\verb$%i$, stop = "
              "\\verb$%i$, length = \\verb$%i$\\\\\n"
              "Longest ds segment (2 MM): start = \\verb$%i$, stop = "
              "\\verb$%i$, length = \\verb$%i$\\\\\n"
              "Total fraction of paired nucleotides = %5.4f\n",
          cand->mfe, cand->pvalue, cand->mean, cand->sd,
          cand->external_loop_count, cand->stem_start, cand->stem_end,
          cand->stem_end - cand->stem_start + 1, cand->stem_start_with_mismatch,
          cand->stem_end_with_mismatch,
          cand->stem_end_with_mismatch - cand->stem_start_with_mismatch + 1,
          cand->paired_fraction);

  u32 *cov_list = chrom_cov->coverage_plus;
  if (cand->strand == '-') {
    cov_list = chrom_cov->coverage_minus;
  }

  u64 total_coverage = 0;
  get_coverage_in_range(&total_coverage, cov_list, cand->start, cand->end);

  fprintf(fp,
          "\\section*{Coverage}\n"
          "Mean coverage precursor miRNA: \\verb$%5.4f$ per nucleotide\\\\\n"
          "Mean coverage mature miRNA: \\verb$%5.4f$ per nucleotide\\\\\n"
          "Mean coverage star miRNA: \\verb$%5.4f$ per nucleotide\\\\\n"
          "Total number of reads for precursor miRNA:\\\\\\n",
          (double)total_coverage / (cand->end - cand->start),
          (double)mature_mirna->coverage /
              (mature_mirna->end - mature_mirna->start),
          (double)star_mirna->coverage / (star_mirna->end - star_mirna->start));
  if (cov_plot_file != NULL) {
    fprintf(fp, "\\begin{center}\n"
                "\\includegraphics[width=\\textwidth]{%s}\n"
                "\\end{center}\\newpage\n",
            cov_plot_file);
  } else {
    fprintf(fp, "{[}Note: Coverage Plot was not included because gnuplot could "
                "not be found or the gnuplot call failed{]}\\\\\n \\\\\n"
                "\\newpage\n");
  }

  fprintf(fp, "\\subsection*{Read alignment: Mature sequence}\n");
  write_unique_reads_to_tex(fp, cand, mature_mirna, mature_reads);
  fprintf(fp, "\\newpage\n");
  fprintf(fp, "\\subsection*{Read alignment: Star sequence}\n");
  fprintf(fp, "\\section*{Structure}\n");
  write_unique_reads_to_tex(fp, cand, star_mirna, star_reads);
  if (structure_file != NULL) {
    fprintf(fp, "\\begin{center}\n"
                "\\includegraphics[height=\\textheight,width=\\textwidth,"
                "keepaspectratio]{%s}\n"
                "\\end{center}\n",
            structure_file);
  } else {
    fprintf(fp, "{[}Note: Structure was not included because java could "
                "not be found or the Varna call failed{]}\\\\\n \\\\\n");
  }
  if (coverage_file != NULL) {
    fprintf(fp, "\\begin{center}\n"
                "\\includegraphics[height=\\textheight,width=\\textwidth,"
                "keepaspectratio]{%s}\n"
                "\\end{center}\n",
            coverage_file);
  } else {
    fprintf(fp,
            "{[}Note: Coverage Structure was not included because java could "
            "not be found or the Varna call failed{]}\\\\\n \\\\\n");
  }
  fprintf(fp, "\\end{document}");

  fclose(fp);
  *tex_file = file_path;
  return E_SUCCESS;
}

int write_unique_reads_to_tex(FILE *fp, struct micro_rna_candidate *cand,
                              struct candidate_subsequence *subseq,
                              struct unique_read_list *reads) {
  size_t total_length = cand->end - cand->start;
  const int PRINT_FLANK = 30;
  u32 start;
  if (subseq->start < PRINT_FLANK) {
    start = 0;
  } else {
    start = subseq->start - PRINT_FLANK;
  }
  u32 end = subseq->end + PRINT_FLANK;
  if (end > total_length) {
    end = total_length - 1;
  }
  fprintf(fp, "\\verb$");
  for (u32 i = start; i < end; i++) {
    if (i >= subseq->start && i < subseq->end) {
      fprintf(fp, "%c", cand->sequence[i]);
    } else {
      fprintf(fp, ".");
    }
  }
  fprintf(fp, "$\\\\\n \\verb$");
  for (u32 i = start; i < end; i++) {
    fprintf(fp, "%c", cand->sequence[i]);
  }
  fprintf(fp, "$\\\\\n");
  struct unique_read *read = NULL;
  for (size_t i = 0; i < reads->n; i++) {
    read = reads->reads[i];
    fprintf(fp, "\\verb$");
    for (u32 i = start; i < end; i++) {
      if (i >= read->start - cand->start && i < read->end - cand->start) {
        u32 local_index = i - (read->start - cand->start);
        fprintf(fp, "%c", read->seq[local_index]);
      } else {
        fprintf(fp, ".");
      }
    }
    fprintf(fp, "$[$\\times$ %i]\\\\\n", read->count);
  }

  return E_SUCCESS;
}

int compile_tex_file(const char *tex_file_path, const char *output_path) {
  char latex_system_call[2048];
  sprintf(latex_system_call,
          "pdflatex -interaction=batchmode -output-directory=\"%s\" %s",
          output_path, tex_file_path);
  int err = system(latex_system_call);
  if (err != 0) {
    return E_LATEX_SYSTEM_CALL_FAILED;
  }
  return E_SUCCESS;
}

int map_coverage_to_color_index(u32 *result, u32 coverage) {
  if (coverage < 1) {
    coverage = 1;
  }
  double log_cov = log10(coverage);
  int index = (int)floor(log_cov * 5.0 + 0.5);
  if (index < 0) {
    index = 0;
  }
  if (index > 8) {
    index = 8;
  }
  *result = index;
  return E_SUCCESS;
}

int write_bed_lines(FILE *fp, struct extended_candidate *ecand) {
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence *mature_mirna = ecand->mature_micro_rna;
  struct candidate_subsequence *star_mirna = ecand->star_micro_rna;
  struct unique_read_list *mature_reads = ecand->mature_reads;
  struct unique_read_list *star_reads = ecand->star_reads;

  u32 mature_read_count = 0;
  for (size_t i = 0; i < mature_reads->n; i++) {
    mature_read_count += mature_reads->reads[i]->count;
  }
  u32 star_read_count = 0;
  for (size_t i = 0; i < star_reads->n; i++) {
    star_read_count += star_reads->reads[i]->count;
  }

  fprintf(fp, "%s\t%llu\t%llu\tCluster_%lld_precursor\t%d\t%c\t%llu\t%llu\t%"
              "d\t%d\t%d\t%7.5f\t%9.7e\t%7.5e\t%7.5e\n",
          cand->chrom, cand->start, cand->end, cand->id, 0, cand->strand,
          cand->start, cand->end, 0, 0, 0, cand->mfe, cand->pvalue, cand->mean,
          cand->sd);
  fprintf(fp, "%s\t%llu\t%llu\tCluster_%lld_%s_[mature]\t%d\t%c\t%llu\t%llu\t%"
              "d\t%d\t%d\n",
          cand->chrom, cand->start, cand->end, cand->id,
          (mature_mirna->arm == 5) ? "5p" : "3p", 0, cand->strand,
          cand->start + mature_mirna->start, cand->start + mature_mirna->end, 0,
          mature_read_count, 0);
  fprintf(fp, "%s\t%llu\t%llu\tCluster_%lld_%s_[star]\t%d\t%c\t%llu\t%llu\t%"
              "d\t%d\t%d\n",
          cand->chrom, cand->start, cand->end, cand->id,
          (mature_mirna->arm == 5) ? "3p" : "5p", 0, cand->strand,
          cand->start + star_mirna->start, cand->start + star_mirna->end, 0,
          star_read_count, 0);
  return E_SUCCESS;
}

int write_json_entry(FILE *fp, struct extended_candidate *ecand) {
  struct micro_rna_candidate *cand = ecand->cand;
  struct candidate_subsequence *mature_mirna = ecand->mature_micro_rna;
  struct candidate_subsequence *star_mirna = ecand->star_micro_rna;

  fprintf(fp, "\"Cluster_%lld_%s\":{\n", cand->id,
          (cand->strand == '-') ? "minus" : "plus");
  fprintf(fp, "\t\"chromosome\":\"%s\",\n", cand->chrom);
  fprintf(fp, "\t\"candidate_start\":%lld,\n", cand->start);
  fprintf(fp, "\t\"candidate_end\":%lld,\n", cand->end);
  fprintf(fp, "\t\"candidate_length\":%lld,\n", cand->end - cand->start);
  fprintf(fp, "\t\"strand\":\"%c\",\n", cand->strand);
  fprintf(fp, "\t\"mature_sequence\":\"");
  for (u32 i = mature_mirna->start; i < mature_mirna->end; i++) {
    fprintf(fp, "%c", cand->sequence[i]);
  }
  fprintf(fp, "\",\n");
  fprintf(fp, "\t\"mature_start\":%lld,\n", cand->start + mature_mirna->start);
  fprintf(fp, "\t\"mature_end\":%lld,\n", cand->start + mature_mirna->end);
  fprintf(fp, "\t\"arm\":%d,\n", mature_mirna->arm);

  fprintf(fp, "\t\"star_sequence\":\"");
  for (u32 i = star_mirna->start; i < star_mirna->end; i++) {
    fprintf(fp, "%c", cand->sequence[i]);
  }
  fprintf(fp, "\",\n");
  fprintf(fp, "\t\"star_start\":%lld,\n", cand->start + star_mirna->start);
  fprintf(fp, "\t\"star_end\":%lld,\n", cand->start + star_mirna->end);
  fprintf(fp, "\t\"mfe\":\"%7.5lf kcal/mol/bp\",\n", cand->mfe);
  fprintf(fp, "\t\"paired_fraction\":%7.5lf,\n", cand->paired_fraction);
  fprintf(fp, "\t\"mfe_mean\":%7.5lf,\n", cand->mean);
  fprintf(fp, "\t\"mfe_sd\":%7.5lf,\n", cand->sd);
  fprintf(fp, "\t\"mfe_precise\":%8lf,\n", cand->mfe);
  fprintf(fp, "\t\"z_value\":%7.5lf,\n", (cand->mfe - cand->mean) / cand->sd);
  fprintf(fp, "\t\"p_value\":%7.5le,\n", cand->pvalue);
  fprintf(fp, "}");
  return E_SUCCESS;
}

int cleanup_auxiliary_files(char *cov_plot_file, char *structure_file,
                            char *coverage_file, char *tex_file,
                            struct configuration_params *config) {
  if (config->cleanup_auxiliary_files != 0) {
    if (cov_plot_file != NULL) {
      remove(cov_plot_file);
    }
    if (structure_file != NULL) {
      remove(structure_file);
    }
    if (coverage_file != NULL) {
      remove(coverage_file);
    }
    if (tex_file != NULL) {
      remove(tex_file);
      char *dot = strchr(tex_file, '.');
      if (dot != NULL) {
        *dot = 0;
        char tmp_file[1024];
        sprintf(tmp_file, "%s.aux", tex_file);
        remove(tmp_file);
        sprintf(tmp_file, "%s.log", tex_file);
        remove(tmp_file);
      }
    }
  }

  if (cov_plot_file != NULL) {
    free(cov_plot_file);
  }
  if (structure_file != NULL) {
    free(structure_file);
  }
  if (coverage_file != NULL) {
    free(coverage_file);
  }
  if (tex_file != NULL) {
    free(tex_file);
  }
  return E_SUCCESS;
}
