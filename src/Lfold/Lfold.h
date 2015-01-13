#include <stddef.h>

#ifndef __VIENNA_RNA_PACKAGE_LFOLD_H__
#define __VIENNA_RNA_PACKAGE_LFOLD_H__

/* custom methods */

struct secondary_structure {
  double mfe;
  int start;
  char *structure_string;
};

struct structure_list {
  struct secondary_structure **structures;
  size_t capacity;
  size_t n;
};

int create_structure_list(struct structure_list **list);
int save_secondary_structure(struct structure_list *list,
                             char *structure_string, double mfe, int start);

void free_structure_list(struct structure_list *list);
void free_secondary_structure(struct secondary_structure *s);

/**
 */

/**
 *  \addtogroup local_fold
 *
 *  Local structures can be predicted by a modified version of the
 *  fold() algorithm that restricts the span of all base pairs.
 *  @{
 *    \file Lfold.h
 *    \brief Predicting local MFE structures of large sequences
 *
 *  @}
 */

/**
 *  \addtogroup local_mfe_fold
 *  @{
 *
 *  @}
 */

/**
 *  \brief The local analog to fold().
 *
 *  Computes the minimum free energy structure including only base pairs
 *  with a span smaller than 'maxdist'
 *
 *  \ingroup local_mfe_fold
 *
 *  \param string
 *  \param structure
 *  \param maxdist
 */
float Lfold(struct structure_list **result, const char *string, int maxdist);

/**
 *  \brief
 *
 *  \ingroup local_mfe_fold
 *
 *  \param string
 *  \param structure
 *  \param maxdist
 *  \param zsc
 *  \param min_z
 */
float Lfoldz(struct structure_list **result, const char *string, int maxdist,
             int zsc, double min_z);

/**
 *  \addtogroup local_consensus_fold
 *  @{
 *
 *  @}
 */

/**
 *  \brief
 *
 *  \ingroup local_consensus_fold
 *
 *  \param strings
 *  \param structure
 *  \param maxdist
 *  \return
 */
float aliLfold(const char **strings, char *structure, int maxdist);

#endif
