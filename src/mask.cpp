//
// Created by shan on 2020-04-22.
//
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <getopt.h>
static inline void strtok_null_check(char *tmp, int line_num){
    // TODO: the file no is passed. have to change
    if(tmp == NULL){
        fprintf(stderr,"Malformed depths file? An offending value is at line number %d\n",line_num);
        exit(EXIT_FAILURE);
    }
}

// TODO Add short description to OPTIONS
static const char *ARTIC_C_MASK_USAGE_MESSAGE =
        "Usage: artic_c mask [OPTIONS] -i [input file] -o [output bedfile]\n"
        "   -d  INT threshold depth\n"
        "   -v  Show Verbose logs\n";

int mask_main(int argc, char ** argv){

    char* ARGS_DEPTHS_FILE = NULL;
    char* ARGS_BED_FILE = NULL;
    int VERBORSITY = 0;
    int DEPTH_THRESHOLD = 20;

    const char* optstring = "vd:i:o:";
    if(argc == 1){
        fprintf (stderr, "%s", ARTIC_C_MASK_USAGE_MESSAGE);
        exit(EXIT_FAILURE);
    }
    int c;
    while ((c = getopt (argc, argv, optstring)) != -1)
        switch (c) {
            case 'd':
                DEPTH_THRESHOLD = atoi(optarg);
                break;
            case 'v':
                VERBORSITY = 1;
                break;
            case 'i':
                ARGS_DEPTHS_FILE = optarg;
                break;
            case 'o':
                ARGS_BED_FILE = optarg;
                break;
            default:
                fprintf (stderr, "%s", ARTIC_C_MASK_USAGE_MESSAGE);
                exit(EXIT_FAILURE);
        }

    if (ARGS_DEPTHS_FILE == NULL || ARGS_BED_FILE == NULL) {
        fprintf(stderr, "%s", ARTIC_C_MASK_USAGE_MESSAGE);
        exit(EXIT_FAILURE);
    }

    // open tab separated depths file
    FILE * fin = fopen(ARGS_DEPTHS_FILE, "r"); // read mode
    if (fin == NULL){
        fprintf(stderr,"Error while opening %s file.\n",ARGS_DEPTHS_FILE);
        exit(EXIT_FAILURE);
    }
    // open bed file
    FILE * fout = fopen(ARGS_BED_FILE, "w"); // write mode
    if (fout == NULL){
        fprintf(stderr,"Error while opening %s file.\n",ARGS_BED_FILE);
        exit(EXIT_FAILURE);
    }
    char *buf = NULL;
    size_t buf_size = 0;
    int line = 0;
    while ((getline(&buf, &buf_size, fin)) != -1) {
        //chromosome
        char * tmp = strtok(buf, "\t");
        strtok_null_check(tmp,line);
        char* chromosome = strdup(tmp);
        tmp = strtok(NULL, "\t");
        strtok_null_check(tmp,line);
        uint32_t position = atoi(tmp);
        tmp = strtok(NULL, "\t");
        strtok_null_check(tmp,line);
        uint32_t depth = atoi(tmp);
        if(depth < DEPTH_THRESHOLD){
            fprintf(fout,"%s\t%u\t%u\n",chromosome,position-1,position);
        }
        free(chromosome);
        line++;
    }


    fprintf(stderr,"hellow mask\n");
    return 0;

}

//    double realtime0 = realtime();
//
//    signal(SIGSEGV, sig_handler);
