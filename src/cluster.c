#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "cluster.h"
#include "parse_sam.h"
#include "errors.h"
#include "string.h"

int cluster(int argc,char** argv)
{
    int c;
    char* test=NULL;
    while((c = getopt(argc, argv, "g:w:r:l:o:h")) != -1){
        switch(c){
        case 'g':
            printf("A set\n");
            break;
        case 'w':
            printf("B set\n");
            break;
        case 'r':
            test= optarg;
            printf("C set: %s\n", test);
            break;
        case 'l':
            break;
        case 'o':
            break;
        case 'h':
            break;
        default:
            break;
        }
    }
    if(optind >= argc){ /* missing input file */
        printf("No SAM File specified\n\n");
        print_help();
        return E_NO_FILE_SPECIFIED;
    }
    struct sam_file* sam = NULL;
    int err = parse_sam(&sam,argv[optind]);
    if(err == E_FILE_NOT_FOUND){
        printf("File %s not found\n",argv[optind]);
        return E_FILE_NOT_FOUND;
    }
    if(err != E_SUCCESS){
        printf("Unknown Error\n");
        return E_UNKNOWN;
    }

    struct cluster_list* list= NULL;
    err = create_clusters(&list,sam);
    if(err != E_SUCCESS){
        print_error(err);
        free_sam(sam);
        return err;
    }
    free_sam(sam);
    sort_clusters(list);
    
    
    free_clusters(list);

    return E_SUCCESS;
}



int print_help()
{
    printf(
        "Description:\n"
        "    cluster generates a list of main expression contigs based on\n"
        "    alignment data.\n"
        "Usage: miRA cluster [-g gapsize] [-w windowsize] [-f minreadnumber]\n"
        "                    [-l maxlength] [-o outputfile] [-h] <input SAM file> \n"
        );
    return E_SUCCESS;
}

int create_clusters(struct cluster_list** list,struct sam_file* sam){
    struct cluster_list* tmp_list = (struct cluster_list*)malloc(sizeof(struct cluster_list));
    if(tmp_list == NULL){
        return E_MALLOC_FAIL;
    }
    tmp_list->capacity = sam->n;
    tmp_list->clusters = (struct cluster*)malloc(tmp_list->capacity*sizeof(struct cluster));
    if(tmp_list->clusters == NULL){
        free(tmp_list);
        return E_MALLOC_FAIL;
    }
    size_t n =0;
    for (size_t i = 0; i < sam->n; i++)
    {
        int result = sam_to_cluster(tmp_list->clusters+n,sam->entries+i,n);
        if(result!=E_SUCCESS) continue;
        n++;
    }
    tmp_list->n = n;
    *list = tmp_list;
    return E_SUCCESS;
}
int sort_clusters(struct cluster_list* list){
    qsort(list->clusters,list->n,sizeof(struct cluster),compare_clusters);
    return E_SUCCESS;
}

int compare_clusters(const void* c1,const void* c2){
    struct cluster* cl1 = (struct cluster*)c1;
    struct cluster* cl2 = (struct cluster*)c2;
    int diff;
    diff = cl1->strand - cl2->strand;
    if(diff !=0) return diff;
    diff = strcmp(cl1->chrom,cl2->chrom);
    if(diff !=0) return diff;
    diff = cl1->chrom - cl2->chrom;
    return diff;
}


int sam_to_cluster(struct cluster* cluster,struct sam_entry* entry,long id)
{
    const char POSITIVE_STRAND_SYMBOL = '+';
    const char NEGATIVE_STRAND_SYMBOL = '-';

    cluster->id = id;
    if(entry->flag & REV_COMPLM){
        cluster->strand = NEGATIVE_STRAND_SYMBOL;
    } else {
        cluster->strand = POSITIVE_STRAND_SYMBOL;
    }
    int n = strlen(entry->rname)+1;
    cluster->chrom = (char*)malloc(n*sizeof(char));
    if(cluster->chrom == NULL) return E_MALLOC_FAIL;
    strcpy(cluster->chrom, entry->rname);

    cluster->start = entry->pos;
    cluster->end = entry->pos + strlen(entry->seq);
    cluster->readcount = 1;


    
    return E_SUCCESS;
}

int free_clusters(struct cluster_list* list){
    for(size_t i=0;i<list->n;i++){
        free_cluster(list->clusters+i);
    }
    free(list->clusters);
    free(list);
    return E_SUCCESS;
}
int free_cluster(struct cluster* c){
    free(c->chrom);
    return E_SUCCESS;
}





