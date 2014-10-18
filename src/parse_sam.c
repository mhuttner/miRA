
#include <stdlib.h>
#include <stdio.h>
#include <regex.h>
#include <string.h>
#include "parse_sam.h"

int parse_sam(struct sam_file** sam,char* file)
{
    static const int MAXLINELENGHT = 2048;
    static const int STARTINGSIZE = 100;
    FILE* fp = fopen(file,"r");
    if(fp==NULL){
        return 1;
    }
    char line[MAXLINELENGHT];
    while(fgets(line,sizeof(line),fp) != NULL){
        parse_line(NULL,line,MAXLINELENGHT);
    }
    fclose(fp);

    return 0;
}
int parse_line(struct sam_entry* e, char* line, int maxlength){
    const char seperator = '\t';
    char* start = line;
    char* end = NULL;
    int line_done =0;
    while(!line_done){
        end = strchr(start, seperator);
        if(end == NULL){
            end = strchr(start,'\n');
            line_done=1;
        }
        printf("%ld\n",end-start);
        start = end+1;
    }
    printf("\n\n");
    return 0;
}


int free_sam(struct sam_file* sam)
{
    for(int i=0;i<sam->n;i++){
        free_sam_entry(sam->entries+i);
    }
    return 0;   
}

int free_sam_entry(struct sam_entry* e){
    free(e->qname);
    free(e->rname);
    free(e->cigar);
    free(e->rnext);
    free(e->seq);
    free(e->qual);
    free(e);
    return 0;
};


