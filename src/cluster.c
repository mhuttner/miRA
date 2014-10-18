#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "cluster.h"
#include "parse_sam.h"

int cluster(int argc,char** argv){
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
    }
    parse_sam(NULL,argv[optind]);


    return 0;
}



void print_help(){
    printf(
        "Description:\n"
        "    cluster generates a list of main expression contigs based on\n"
        "    alignment data.\n"
        "Usage: miRA cluster [-g gapsize] [-w windowsize] [-f minreadnumber]\n"
        "                    [-l maxlength] [-o outputfile] [-h] <input SAM file> \n"
        );
}





