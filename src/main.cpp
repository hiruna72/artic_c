//
// Created by hiruna72 on 2020-04-22.
//

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include "misc.h"
#include "error.h"
#include "assert.h"
#include "multiIntersectBed/multiIntersectBed.h"
#include "interface_artic_c.h"

#define ARTIC_C_VERSION "0.1"

//make the segmentation faults a bit cool
void sig_handler(int sig) {
    fprintf(stderr,"I regret to inform that a segmentation fault occurred. But at least "
                   "it is better than a wrong answer%s",
            ".");
    exit(EXIT_FAILURE);
}
int trim_main(int argc, char ** argv);
int mask_main(int argc, char ** argv);
int multiintersect_main(int argc, char* argv[]); // since bedtools is big, copied only the mutliintersect part and necessary utils

int print_usage(FILE *fp_help){

    fprintf(fp_help,"Usage: artic_c <command> [options]\n\n");
    fprintf(fp_help,"command:\n");
    fprintf(fp_help,"         trim    downsize dataset by trimming alignments\n");
    fprintf(fp_help,"         multiinter    merge bedfiles (copied from bedtools multiintersect)\n");
    fprintf(fp_help,"         mask    create a bed file containing bases with lower depth\n\n");

    if(fp_help==stderr){
        exit(EXIT_FAILURE);
    }
    else if(fp_help==stdout){
        exit(EXIT_SUCCESS);
    }
    else{
        assert(0);
    }


}


int init_artic(int argc, char **argv){

    double realtime0 = realtime();

    signal(SIGSEGV, sig_handler);

    int ret=1;

    if(argc<2){
        return print_usage(stderr);
    }
    if(strcmp(argv[1],"trim")==0){
        ret=trim_main(argc-1, argv+1);
    }
    else if(strcmp(argv[1],"mask")==0){
        ret=mask_main(argc-1, argv+1);
    }
    else if(strcmp(argv[1],"multiinter")==0){
        ret=multiintersect_main(argc-1, argv+1);
    }
    else if(strcmp(argv[1],"--version")==0 || strcmp(argv[1],"-V")==0){
        fprintf(stdout,"ARTIC_C %s\n",ARTIC_C_VERSION);
        exit(EXIT_SUCCESS);
    }
    else if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0){
        print_usage(stdout);
    }
    else{
        fprintf(stderr,"[artic_c] Unrecognised command %s\n",argv[1]);
        print_usage(stderr);
    }

    fprintf(stderr, "\n[%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) {
        fprintf(stderr, " %s", argv[i]);
    }

    fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU time: %.3f sec; Peak RAM: %.3f GB\n\n",
            __func__, realtime() - realtime0, cputime(),peakrss() / 1024.0 / 1024.0 / 1024.0);

    return ret;
}
