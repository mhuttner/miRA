#include <stddef.h>

#ifndef PARSE_SAM_H
#define PARSE_SAM_H 

struct sam_file
{
    size_t n;
    size_t capacity;
    struct sam_entry* entries;
};

struct sam_entry
{
    char* qname;
    int flag;
    char* rname;
    long pos;
    int mapq;
    char* cigar;
    char* rnext;
    long pnext;
    long tlen;
    char* seq;
    char* qual;
};

int parse_sam(struct sam_file** sam,char* file);
int parse_line(struct sam_entry* e, char* line, int maxlength);
int free_sam(struct sam_file* sam);
int free_sam_entry(struct sam_entry* e);

#endif