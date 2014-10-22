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

enum sam_flag {
    MULT_SEG = 0x1,
    PROP_ALIGN = 0x2,
    UNMAPPED = 0x4,
    NEXT_UNMAPPED = 0x8,
    REV_COMPLM = 0x10,
    NEXT_REV_COMPLM = 0x20,
    FIRST_SEG = 0x40,
    LAST_SEG = 0x80,
    SECONDARY_ALIGN = 0x100,
    BAD_QUALITY = 0x200,
    DUPLICATE = 0x400,
    SUPPLEMENTARY = 0x800
};

int parse_sam(struct sam_file** sam,char* file);
int parse_line(struct sam_entry* e, char* line, int maxlength);
int free_sam(struct sam_file* sam);
int free_sam_entry(struct sam_entry* e);

#endif