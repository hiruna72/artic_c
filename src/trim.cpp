//
// Created by shan on 2020-04-02.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include "artic.h"
#include "htslib/sam.h"
#include "common.h"
#include "misc.h"

// default verbosity
int ARGS_VERBORSE = 0;

//long no_end_zero_entries = 0;
//long no_end_one_entries = 0;
std::set<std::string> pools;

// consumesReference lookup for if a CIGAR operation consumes the reference sequence
uint32_t consumeReference [] = {1,0,1,1,0,0,0,1};
// consumesQuery lookup for if a CIGAR operation consumes the query sequence
uint32_t consumesQuery [] = {1,1,0,0,1,0,0,1};



static inline void strtok_null_check(char *tmp, int line_num){
    // TODO: the file no is passed. have to change
    if(tmp == NULL){
        fprintf(stderr,"Malformed primer scheme file? An offending value is at line number %d\n",line_num);
        exit(EXIT_FAILURE);
    }
}

void merge_sites(bed_row* cannonical, std::vector<bed_row>::iterator alt){
    if(cannonical->direction_sign != alt->direction_sign){
        fprintf(stderr,"could not merge alt with different orientation to canonical\n");
        exit(EXIT_FAILURE);
    }
    //  merge the start/ends of the alt with the canonical to get the largest window possible
    if(cannonical->direction_sign == 1) {
        if (alt->start < cannonical->start) {
            cannonical->start = alt->start;
        }
        if (alt->end > cannonical->end) {
            cannonical->end = alt->end;
        }
    } else{
        if(alt->start > cannonical->start){
            cannonical->start = alt->start;
        }
        if(alt->end < cannonical->end){
            cannonical->end = alt->end;
        }
    }
}


std::vector<bed_row> read_bed_file(const char *bed_file_name) {
    /*Parses a bed file and collapses alts into canonical primer sites
    Parameters
    ----------
            fn : str
    The bedfile to parse
    Returns
    -------
            list
    A list of dictionaries, where each dictionary contains a row of the parsed bedfile.
            The available dictionary keys are - Primer_ID, direction, start, end
    */
    // create a vector that holds values to be returned. some values in ths vector will get updated but will not be deleted
    std::vector<bed_row> bed_file;
    std::vector<bed_row> alt_bed_file;
    char * alt_check = "_alt";
    char * left_check = "LEFT";
    char * right_check = "RIGHT";

    FILE * fin = fopen(bed_file_name, "r"); // read mode
    if (fin == NULL){
        fprintf(stderr,"Error while opening %s file.\n",bed_file_name);
        exit(EXIT_FAILURE);
    }
    char *buf = NULL;
    size_t buf_size = 0;
    int line = 0;
    int original_lines = 0;

    std::map<const char*,int> primerID_index;
    pools.insert(std::string("unmatched"));

    while ((getline(&buf, &buf_size, fin)) != -1) {
//        printf("%s", buf);
        bed_row* row;
        if(strstr(buf, alt_check) == NULL) {
            bed_file.emplace_back();
            row = &bed_file.back();
            row->original = 1;
        } else{
            alt_bed_file.emplace_back();
            row = &alt_bed_file.back();
            row->original = 0;
        }
        //chromosome
        char * tmp = strtok(buf, "\t");
        strtok_null_check(tmp,line);
        row->chromosome = strdup(tmp);
        tmp = strtok(NULL, "\t");
        strtok_null_check(tmp,line);
        row->start = atoi(tmp);
        tmp = strtok(NULL, "\t");
        strtok_null_check(tmp,line);
        row->end = atoi(tmp);
        tmp = strtok(NULL, "\t");
        strtok_null_check(tmp,line);
        row->Primer_ID = strdup(tmp);
        tmp = strtok(NULL, "\t");
        strtok_null_check(tmp,line);
        row->PoolName = strdup(tmp);
        if(strstr(row->Primer_ID, left_check) != NULL) {
            row->direction = strdup("+");
            row->direction_sign = 1;
        } else if (strstr(row->Primer_ID, right_check) != NULL){
            row->direction = strdup("-");
            row->direction_sign = -1;
        } else{
            fprintf(stderr,"LEFT/RIGHT must be specified in Primer ID - %s\n",row->Primer_ID);
            exit(EXIT_FAILURE);
        }
        line ++;
        if(row->original){
            if(primerID_index.count(row->Primer_ID)==0){
                primerID_index[row->Primer_ID] = original_lines;
                original_lines++;
            }else{
                // verify_integrity is used to prevent duplicate Primer_IDs being processed
                bed_file.pop_back();
            }
        }
        pools.insert(std::string(row->PoolName));

    }
    if(bed_file.empty() and alt_bed_file.empty()){
        fprintf(stderr,"primer scheme file is empty\n");
        exit(EXIT_FAILURE);
    }

    /*
     *  for _, row in alts.iterrows():
        primerID = row['Primer_ID'].split('_alt')[0]

        # check the bedFile if another version of this primer exists
        if primerID not in bedFile:

            # add to the bed file and continue
            bedFile[primerID] = row
            continue

        # otherwise, we've got a primer ID we've already seen so merge the alt
        mergedSite = merge_sites(bedFile[primerID], row)

        # update the bedFile
        bedFile[primerID] = mergedSite

     */
    // print map
//    for (auto it=primerID_index.begin(); it!=primerID_index.end(); ++it)
//        std::cout << it->first << " => " << it->second << '\n';

// print pools set
//    for (auto it=pools.begin(); it!=pools.end(); ++it)
//        std::cout << ' ' << *it;

    for (auto it=alt_bed_file.begin(); it!=alt_bed_file.end(); ++it){
        std::cout << it->Primer_ID << std::endl;
        char * key = strdup(it->Primer_ID);
        char * p = strstr(key,alt_check);
        *p = '\0';
        std::cout << key<< std::endl;
        if (primerID_index.count(key)){
            bed_file.emplace_back();
            int index = primerID_index.find(key)->second;
            merge_sites(&bed_file[index],it);
        } else{
            bed_file.emplace_back();
            bed_row* row = &bed_file.back();
            row->chromosome = strdup(it->chromosome);
            row->start = it->start;
            row->end = it->end;
            row->Primer_ID = strdup(it->Primer_ID);
            row->PoolName = strdup(it->PoolName);
            row->direction = strdup(it->direction);
            row->direction_sign = it->direction_sign;
            row->original = 1;
        }
    }

//    FILE* fout = fopen("bed_file_produced.txt","w");
//    for (auto it=bed_file.begin(); it!=bed_file.end(); ++it){
////        std::cout << it->chromosome << it->start << it->end << it->Primer_ID << it->PoolName << it->direction << std::endl;
//        fprintf(fout,"%s\t%d\t%d\t%s\t%s\t%s\n",it->chromosome,it->start,it->end,it->Primer_ID,it->PoolName,it->direction);
//    }
//    fclose(fout);
    return bed_file;
}

int find_primer(std::vector<bed_row> bed, uint32_t pos, int direction) {
    uint32_t min_distance = INT32_MAX;
    size_t index = 0;
    int i = 0;
    for(auto it=bed.begin();it!=bed.end();++it){
        if(it->direction_sign == direction and direction == 1){
            uint32_t distance = (it->start > pos) ? it->start - pos : pos - it->start;
            if(distance < min_distance){
                min_distance = distance;
                index = i;
            }
        }
        if(it->direction_sign == direction and direction == -1){
            uint32_t distance = (it->end > pos)? it->end - pos : pos - it->end;
            if( distance < min_distance){
                min_distance = distance;
                index = i;
            }
        }
        i++;
    }
    return index;
}

int trim(std::vector<uint32_t>*cigar,bam1_t *b, uint32_t end_pos,uint32_t primer_pos, int end) {

    uint32_t pos;
    if(end == 0){
        pos = b->core.pos;
        std::reverse(cigar->begin(),cigar->end()); // now it is easy to pop and append
    } else{
        pos = end_pos;
    }

    uint32_t eaten = 0;
    uint32_t iter = 0;
    while(1){
        iter++;
        uint32_t flag;
        uint32_t length;
        // chomp stuff until we reach pos
        if(cigar->size() >= 2){
            if(end){
                length = cigar->back();
                cigar->pop_back();
                flag = cigar->back();
                cigar->pop_back();
            }else{
                flag = cigar->back();
                cigar->pop_back();
                length = cigar->back();
                cigar->pop_back();
            }
        }else{
            if (ARGS_VERBORSE == 1) {
                fprintf(stderr,"Ran out of cigar during soft masking - completely masked read will be ignored\n");
            }
            break;
        }

// if the CIGAR operation consumes the reference sequence, increment/decrement the position by the CIGAR operation length
        if(consumeReference[flag]){
            if(end == 0){
                pos += length;
            } else{
                pos -= length;
            }
        }
// if the CIGAR operation consumes the query sequence, increment the number of CIGAR operations eaten by the CIGAR operation length
        if(consumesQuery[flag]){
            eaten += length;
        }
// stop processing the CIGAR if we've gone far enough to mask the primer
        if(end == 0 and pos >= primer_pos and flag == 0){
            break;
        }
        if(end == 1 and pos <= primer_pos and flag == 0){
            break;
        }
    }

    uint32_t extra = (pos > primer_pos)?pos-primer_pos:primer_pos-pos;
    if(extra){
        if(end){
            cigar->push_back(0);
            cigar->push_back(extra);
        }else{
            cigar->push_back(extra);
            cigar->push_back(0);
        }
        eaten -= extra;
    }

    // softmask the left primer
    if(end == 0){
        // update the position of the leftmost mappinng base
        b->core.pos = pos - extra;
        // if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if(cigar->back() == 2){
            while(1){
                uint32_t flag = cigar->back();
                if(flag != 2){
                    break;
                }
                cigar->pop_back();
                uint32_t length = cigar->back();
                cigar->pop_back();
                b->core.pos += length;
            }
        }
        // check the new CIGAR and replace the old one
        if(eaten <= 0){
            fprintf(stderr,"invalid cigar operation created - possibly due to INDEL in primer");
            std::reverse(cigar->begin(),cigar->end()); // reverse cigar again
            return 0;
        }
        cigar->push_back(eaten);
        cigar->push_back(4);

    } else{ // softmask the right primer
        // check the new CIGAR and replace the old one
        if(eaten <= 0){
            fprintf(stderr,"invalid cigar operation created - possibly due to INDEL in primer");
            return 0;
        }
        cigar->push_back(4);
        cigar->push_back(eaten);
    }

    // reverse cigar again
    if(end==0){
        std::reverse(cigar->begin(),cigar->end());
//        no_end_zero_entries++;
    }else{
//        no_end_one_entries++;
    }
    return 1;
}

// TODO Add short description to OPTIONS
static const char *ARTIC_C_USAGE_MESSAGE =
        "Usage: artic_c trim [OPTIONS] -b [bed file] -i [input samfile] -o [output samfile]\n"
        "   -n  INT\n"
        "   -r\n"
        "   -s\n"
        "   -g\n"
        "   -v  Show Verbose logs\n";

/////////////////////////////////////////////////////////////////////MAIN//////////////////////////////////////////////////

int trim_main(int argc, char ** argv){

//    double realtime0 = realtime();
//
//    signal(SIGSEGV, sig_handler);

    // args flags
    int ARGS_REMOVE_INCORRECT_PAIRS = 0;
    int ARGS_START = 0;
    int ARGS_NORMALISE = 200;
    int ARGS_NO_READ_GROUPS = 0;
    char* ARGS_SAMFILE_IN = NULL;
    char* ARGS_SAMFILE_OUT = NULL;
    char* ARGS_BEDFILE = NULL;

    const char* optstring = "grsvn:b:i:o:";

    if(argc == 1){
        fprintf (stderr, "%s", ARTIC_C_USAGE_MESSAGE);
        exit(EXIT_FAILURE);
    }
    int c;
    while ((c = getopt (argc, argv, optstring)) != -1)
        switch (c) {
            case 'r':
                ARGS_REMOVE_INCORRECT_PAIRS = 1;
                break;
            case 's':
                ARGS_START = 1;
                break;
            case 'n':
                ARGS_NORMALISE = atoi(optarg);
                break;
            case 'g':
                ARGS_NO_READ_GROUPS = 1;
                break;
            case 'v':
                ARGS_VERBORSE = 1;
                break;
            case 'b':
                ARGS_BEDFILE = optarg;
                break;
            case 'i':
                ARGS_SAMFILE_IN = optarg;
                break;
            case 'o':
                ARGS_SAMFILE_OUT = optarg;
                break;
            default:
                fprintf (stderr, "%s", ARTIC_C_USAGE_MESSAGE);
                exit(EXIT_FAILURE);
        }

    if (ARGS_BEDFILE == NULL || ARGS_SAMFILE_IN == NULL || ARGS_SAMFILE_OUT == NULL) {
        fprintf(stderr, "%s", ARTIC_C_USAGE_MESSAGE);
        exit(EXIT_FAILURE);
    }

    std::vector<bed_row> bed = read_bed_file(ARGS_BEDFILE);

    //these come from htslib/sam.h
    samFile *in = NULL;
    samFile *out = NULL;
    bam1_t *b= NULL;
    bam_hdr_t *header = NULL;

    //open the BAM file (though called sam_open is opens bam files too :P)
    in = sam_open(ARGS_SAMFILE_IN, "r");
    errorCheckNULL(in);

    fprintf(stderr,"SM out %s\n",ARGS_SAMFILE_OUT);
    out = sam_open(ARGS_SAMFILE_OUT, "w");      //for writing
    errorCheckNULL(out);

    //get the sam header.
    if ((header = sam_hdr_read(in)) == 0){
        fprintf(stderr,"No sam header?\n");
        exit(EXIT_FAILURE);
    }

    //write the SAM header
    int ret=sam_hdr_write(out,header);
    errorCheck(ret);

    b = bam_init1();

    //my structure for a read (see common.h)
//    struct alignedRead* read = (struct alignedRead*)malloc(sizeof(struct alignedRead));

    std::map<std::string,int[2]>counter;

//    long no_bam_entries = 0;
//    long no_completed_entries = 0;
//    long no_unmapped_entries = 0;
//    long no_supp_entries = 0;
//    long no_incorrect_entries = 0;
//    long no_before_entries = 0;
//    long no_p1_entries = 0;
//    long no_p2_entries = 0;

    //repeat until all reads in the file are retrieved
    while ( sam_read1(in, header, b) >= 0){
//        no_bam_entries++;

//        getRead(read, b);         //copy the current read to the myread structure. See common.c for information
        bam1_core_t *c = &(b->core);
        if(c->flag & BAM_FUNMAP){
            if (ARGS_VERBORSE == 1) {
                fprintf(stderr,"%s read is skipped because it is unmapped\n",bam_get_qname(b));
            }
//            no_unmapped_entries++;
            continue;
        }
        if(c->flag & BAM_FSUPPLEMENTARY){
            if (ARGS_VERBORSE == 1) {
                fprintf(stderr,"%s read is skipped because it is supplementary\n",bam_get_qname(b));
            }
//            no_supp_entries++;
            continue;
        }
        int p1 = find_primer(bed,c->pos,1);
        int p2 = find_primer(bed,bam_endpos(b),-1);

        char* primer_left = strdup(bed[p1].Primer_ID);
        char * p = strstr(primer_left,"LEFT");
        *p = '\0';

        char* primer_right = strdup(bed[p2].Primer_ID);
        p = strstr(primer_right,"RIGHT");
        *p = '\0';
        int correctly_paired = (strcmp(primer_left,primer_right)==0)?1:0; // return 0 if correctly paired

        if(ARGS_REMOVE_INCORRECT_PAIRS==1 and correctly_paired == 0){
//            no_incorrect_entries++;
            if (ARGS_VERBORSE == 1) {
                fprintf(stderr, "%s skipped as not correctly paired\n", bam_get_qname(b));
            }
            continue;
        }
        //report
        if (ARGS_VERBORSE == 1) {
            fprintf(stderr,"%s\t%ld\t%ld\t%s\t%s\n",bam_get_qname(b),c->pos,bam_endpos(b),bed[p1].Primer_ID,bed[p2].Primer_ID);
        }

        // if the alignment starts before the end of the primer, trim to that position
        uint32_t p1_position;
        uint32_t p2_position;
        if(ARGS_START){
            p1_position = bed[p1].start;
            p2_position = bed[p2].end;

        } else{
            p1_position = bed[p1].end;
            p2_position = bed[p2].start;
        }
//        no_before_entries++;

        std::vector<uint32_t> cigar;
        uint32_t *original_cigar = bam_get_cigar(b);
        for(auto i = 0; i < c->n_cigar; i++){
            uint32_t cigarFlag     = bam_cigar_op(original_cigar[i]);
            uint32_t cigarFlagLen  = bam_cigar_oplen(original_cigar[i]);
            cigar.push_back(cigarFlag);
            cigar.push_back(cigarFlagLen);
        }
        uint32_t end_pos = bam_endpos(b);

        if(c->pos < p1_position){
//            no_p1_entries++;
            int trim_success = trim(&cigar, b,end_pos, p1_position, 0);
            if(trim_success == 0){
                continue;
            }
        }

        if(end_pos > p2_position){
//            no_p2_entries++;
            int trim_success = trim(&cigar,b,end_pos, p2_position, 1);
            if(trim_success == 0){
                continue;
            }
        }

        // todo: exception handle
        if(ARGS_NORMALISE){
            char key[100];
            strcpy (key,bed[p1].Primer_ID);
            strcat(key,bed[p2].Primer_ID);
            int index = (c->flag & BAM_FREVERSE) ? 1:0;
            counter[key][index]++;
            if(counter[key][index] > ARGS_NORMALISE){
                continue;
            }
        }

        uint32_t new_cigar_len = cigar.size()/2;
        // check the the alignment still contains bases matching the reference
        for (auto i=0 ;i < new_cigar_len; i++){
            if(cigar[2*i] == 0){
                continue;
            }
        }

//        no_completed_entries++;

        // write back to b
        // we only have to update cigar in b with std::vector<> cigar and b->core.pos with read->pos
        uint32_t *new_cigar_ptr = (uint32_t*)malloc(sizeof(uint32_t)*new_cigar_len);
        uint32_t *new_cigar = new_cigar_ptr;
        for(auto it = cigar.begin(); it != cigar.end(); ++it) {
            uint32_t flag = *it;
            ++it;
            uint32_t length = *it;
            length = length << BAM_CIGAR_SHIFT;
            *new_cigar_ptr = (length | flag);
            new_cigar_ptr++;
        }

        //assert(new_cigar_len<=b->core.n_cigar);
       
        uint8_t *seq, *qual, *pointer;
 
        bam1_t *b_dup = bam_dup1(b);
        errorCheckNULL(b_dup);

        // allocate memory for the new CIGAR
        if (b->l_data + (new_cigar_len + 1) * 4 > b->m_data) { // not enough memory
            b_dup->m_data = b->l_data + new_cigar_len * 4;
            kroundup32(b_dup->m_data);
            b_dup->data = (uint8_t*)realloc(b_dup->data, b_dup->m_data);
            if (ARGS_VERBORSE == 1){
                fprintf(stderr,"reallocated\n");
            }
        }

        seq = bam_get_seq(b); 
        qual = bam_get_qual(b);
        int j = b->core.l_qseq;

        memcpy(bam_get_cigar(b_dup), new_cigar, new_cigar_len * 4); // set CIGAR
        pointer = b_dup->data + b_dup->core.l_qname + new_cigar_len * 4;
        memcpy(pointer, seq, (j+1)>>1); pointer += (j+1)>>1; // set SEQ 
        memcpy(pointer, qual, j); pointer += j; // set QUAL
        memcpy(pointer, bam_get_aux(b), bam_get_l_aux(b)); pointer += bam_get_l_aux(b); // set optional fields
        b_dup->core.n_cigar = new_cigar_len;
        // b->core.l_qseq = j; // update CIGAR length and query length
        b_dup->l_data = pointer - b_dup->data; // update record length

        b_dup->core.pos= c->pos;

//        TODO: modify header to include 'RG' tag and correct the code below
//        if(ARGS_NO_READ_GROUPS == 0){
//            if(correctly_paired){
//                const char* pool_name = bed[p1].PoolName;
//                int len = strlen(pool_name);
//                const uint8_t* p = reinterpret_cast<const uint8_t*>(pool_name);
//                bam_aux_append(b_dup,"RG",'Z',len,p);
//            }else{
//                const char* pool_name = "unmatched";
//                int len = strlen(pool_name);
//                const uint8_t* p = reinterpret_cast<const uint8_t*>(pool_name);
//                bam_aux_append(b_dup,"RG",'Z',len,p);
//            }
//        }

        ret=sam_write1(out,header,b_dup);
        errorCheck(ret);
        
        bam_destroy1(b_dup);

//        free(new_cigar_ptr);

    }

    //wrap up
//    free(read);
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out);

    std::cout << "C/C++" << std::endl;
//    std::cout << "no_of_bam_entries "<< no_bam_entries << std::endl;
//    std::cout << "no_unmapped_entries "<< no_unmapped_entries << std::endl;
//    std::cout << "no_supp_entries "<< no_supp_entries << std::endl;
//    std::cout << "no_incorrect_entries "<< no_incorrect_entries << std::endl;
//    std::cout << "no_before_entries "<< no_before_entries << std::endl;
//    std::cout << "no_p1_entries "<< no_p1_entries << std::endl;
//    std::cout << "no_end_zero_entries "<< no_end_zero_entries << std::endl;
//    std::cout << "no_p2_entries "<< no_p2_entries << std::endl;
//    std::cout << "no_end_one_entries "<< no_end_one_entries << std::endl;
//    std::cout << "no_completed_entries "<< no_completed_entries << std::endl;

//    fprintf(stderr, "[%s] Real time: %.3f sec; CPU time: %.3f sec; Peak RAM: %.3f GB\n\n",
//         __func__, realtime() - realtime0, cputime(),peakrss() / 1024.0 / 1024.0 / 1024.0);

    return 0;



}

// trim -n 200 -g -b ../test_files/nCoV-2019.bed -i ../test_files/SP1-mapped.bam -o trimmed.bam