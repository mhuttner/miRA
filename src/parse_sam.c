
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "parse_sam.h"
#include "errors.h"

int parse_sam(struct sam_file** sam,char* file)
{
    static const int MAXLINELENGHT = 2048;
    static const int STARTINGSIZE = 1024;
    struct sam_file* data = (struct sam_file*)malloc(sizeof(struct sam_file));

    data->capacity = STARTINGSIZE;
    data->entries = (struct sam_entry*)malloc(data->capacity * sizeof(struct sam_entry));
    data->n=0;

    FILE* fp = fopen(file,"r");
    if(fp==NULL){
        return E_FILE_NOT_FOUND;
    }
    char line[MAXLINELENGHT];
    while(fgets(line,sizeof(line),fp) != NULL){
        int result = parse_line(data->entries+data->n,line,MAXLINELENGHT);
        if(result != E_SUCCESS) continue;
        data->n++;
        if(data->n==data->capacity){
            data->capacity *=2;
            data->entries=(struct sam_entry*)realloc(data->entries,data->capacity*sizeof(struct sam_entry));
        }
    }
    fclose(fp);
    *sam = data;

    return E_SUCCESS;
}
int parse_line(struct sam_entry* e, char* line, int maxlength){
    const char seperator = '\t';
    const int num_entries = 11;
    const char header_line_marker = '@';

    char* start = line;
    char* end = NULL;
    if(line[0] == header_line_marker) return E_SAM_HEADER_LINE;
    char** tokens = (char**)malloc(num_entries*sizeof(char*));
    int line_done =0;
    int current_token =0;
    

    while(!line_done){
        if(current_token >= num_entries){
            break;
        }
        end = strchr(start, seperator);
        if(end == NULL){
            end = strchr(start,'\n');
            if(end == NULL) end = strchr(start,'\0');
            line_done=1;
        }
        int l = end - start;
        tokens[current_token] = (char*)malloc((l+1)*sizeof(char));
        memcpy(tokens[current_token],start,l);
        tokens[current_token][l]=0;
        current_token++;
        start = end+1;
    }
    int line_valid = 1;

    if(current_token<num_entries){
        for(int i=0;i<current_token;i++){
            free(tokens[i]);
        }
        free(tokens);
        return E_INVALID_SAM_LINE;
    }

    char* check = NULL;

    e->qname = tokens[0];
    e->rname = tokens[2];
    e->cigar = tokens[5];
    e->rnext = tokens[6];
    e->seq = tokens[9];
    e->qual = tokens[10];

    long flag = strtol(tokens[1],&check,10);
    if(check == tokens[1] || *check!=0){
        line_valid =0;
    }
    e->flag=(int)flag;
    long pos = strtol(tokens[3],&check,10);
    if(check == tokens[3] || *check!=0){
        line_valid =0;
    }
    e->pos=pos;
    long mapq = strtol(tokens[4],&check,10);
    if(check == tokens[4] || *check!=0){
        line_valid =0;
    }
    e->mapq=(int)mapq;
    long pnext = strtol(tokens[7],&check,10);
    if(check == tokens[7] || *check!=0){
        line_valid =0;
    }
    e->pnext=pnext;
    long tlen = strtol(tokens[8],&check,10);
    if(check == tokens[8] || *check!=0){
        line_valid =0;
    }
    e->tlen=tlen;


    if(line_valid){
        free(tokens[1]);
        free(tokens[3]);
        free(tokens[4]);
        free(tokens[7]);
        free(tokens[8]);
        free(tokens);
        return E_SUCCESS;
    } else {
        for(int i=0;i<num_entries;i++){
            free(tokens[i]);
        }
        free(tokens);
        return E_INVALID_SAM_LINE;
    }
}


int free_sam(struct sam_file* sam)
{
    for(int i=0;i<sam->n;i++){
        free_sam_entry(sam->entries+i);
    }
    free(sam->entries);
    return E_SUCCESS;   
}

int free_sam_entry(struct sam_entry* e){
    free(e->qname);
    free(e->rname);
    free(e->cigar);
    free(e->rnext);
    free(e->seq);
    free(e->qual);
    return E_SUCCESS;
};


