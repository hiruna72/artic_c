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
#include <3rdparty/htslib/sam.h>
#include "artic.h"
#include "htslib/sam.h"
#include "common.h"

std::set<std::string> pools;


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
        perror("Error while opening output file.\n");
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

    FILE* fout = fopen("bed_file_produced.txt","w");
    for (auto it=bed_file.begin(); it!=bed_file.end(); ++it){
//        std::cout << it->chromosome << it->start << it->end << it->Primer_ID << it->PoolName << it->direction << std::endl;
        fprintf(fout,"%s\t%d\t%d\t%s\t%s\t%s\n",it->chromosome,it->start,it->end,it->Primer_ID,it->PoolName,it->direction);
    }
    fclose(fout);
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

int trim(std::vector<uint32_t>*cigar,alignedRead *read, int start_pos, int end) {
    int pos;
    if(end==0){
        pos = read->pos;
        std::reverse(cigar->begin(),cigar->end()); // now it is easy to pop and append
    } else{
        pos = read->end;
    }
//    fprintf(stderr,"pos = %d\n",pos);
    int eaten = 0;

//    fprintf(stderr,"before while.... in read %s\n",read->qname);
    int iter = 0;
    while(1){
        iter++;
        uint32_t flag;
        uint32_t length;
        // chomp stuff until we reach pos
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
//        fprintf(stderr,"flag = %d length = %d\n",flag,length);
        if(flag == 0){
            // match
            eaten += length;
            if(end == 0){
                pos += length;
            } else{
                pos -= length;
            }
        }
        if(flag == 1){
            // insertion to the ref
            eaten += length;
        }
        if(flag == 2){
            // deletion to the ref
            if(end == 0){
                pos += length;
            }else{
                pos -= length;
            }
        }
        if(flag == 4){
            eaten += length;
        }
        if(end == 0 and pos >= start_pos and flag == 0){
            break;
        }
        if(end == 1 and pos <= start_pos and flag == 0){
            break;
        }
    }
//    fprintf(stderr,"after while.... in read %s iterations %d\n",read->qname,iter);
    int extra = abs(pos-start_pos);
    if(extra){
        // if flag ==0 // this is always true
        cigar->push_back(0);
        cigar->push_back(extra);
        eaten -= extra;
    }
    if(end==0){
        read->pos = pos - extra;
    }
    if(end==0 and cigar->end()[-2]==2){
        std::reverse(cigar->begin(),cigar->end()); // reverse cigar again
        return 0;
    }

    if(end){
        cigar->push_back(4);
        cigar->push_back(eaten);
    }else{
        cigar->push_back(eaten);
        cigar->push_back(4);
    }

    // reverse cigar again
    if(end==0){
        std::reverse(cigar->begin(),cigar->end());
    }
//    fprintf(stderr,"before returning\n");
    return 1;
}

int check_still_matching_bases(alignedRead *read) {
    for (auto i=0 ;i < read->cigarLen; i++){
        if(read->cigarOps[2*i] == 0){
            return 1;
        }
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////MAIN//////////////////////////////////////////////////

int main(int argc, char ** argv){
    // args flags
    int ARGS_REMOVE_INCORRECT_PAIRS = 1;
    int ARGS_START = 1;
    int ARGS_VERBORSE = 1;
    int ARGS_NORMALISE = 200;
    char* ARGS_SAMFILE_IN = "SP1-mapped.bam";
    char* ARGS_SAMFILE_OUT = "trimmed.bam";
    char* ARGS_BEDFILE = "nCoV-2019.bed";

    std::vector<bed_row> bed = read_bed_file(ARGS_BEDFILE);

    //these come from htslib/sam.h
    samFile *in = NULL;
    samFile *out = NULL;
    bam1_t *b= NULL;
    bam_hdr_t *header = NULL;

    //open the BAM file (though called sam_open is opens bam files too :P)
    in = sam_open(ARGS_SAMFILE_IN, "r");
    errorCheckNULL(in);

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

    //print the chromosome names in the header
    //see the bam_hdr_t struct in htslib/sam.h for parsing the rest of the stuff in header
    int i;
    for(i=0; i< (header->n_targets); i++){
        printf("Chromosome ID %d = %s\n",i,(header->target_name[i]));
    }

    b = bam_init1();

    //my structure for a read (see common.h)
    struct alignedRead* read = (struct alignedRead*)malloc(sizeof(struct alignedRead));

    std::map<std::string,int[2]>counter;

    FILE* fout = fopen("read_names.txt","w");

    //repeat until all reads in the file are retrieved
    while ( sam_read1(in, header, b) >= 0){
        getRead(read, b);         //copy the current read to the myread structure. See common.c for information
//        fprintf(fout,"%s\t%d\t%d\n",read->qname,read->cigarLen, sizeof(read->cigarOps)/ sizeof(uint32_t));
//        continue;
        if(read->flag & BAM_FUNMAP){
//            fprintf(stderr,"%s read is skipped because it is unmapped\n",read->qname);
            continue;
        }
        if(read->flag & BAM_FSUPPLEMENTARY){
//            fprintf(stderr,"%s read is skipped because it is supplementary\n",read->qname);
            continue;
        }
        int p1 = find_primer(bed,read->pos,1);
        int p2 = find_primer(bed,read->end,-1);


        char* primer_left = strdup(bed[p1].Primer_ID);
        char * p = strstr(primer_left,"LEFT");
        *p = '\0';

        char* primer_right = strdup(bed[p2].Primer_ID);
        p = strstr(primer_right,"RIGHT");
        *p = '\0';
        int is_correctly_paired = strcmp(primer_left,primer_right); // return 0 if correctly paired

        if(ARGS_REMOVE_INCORRECT_PAIRS==1 and is_correctly_paired!=0){
            fprintf(fout,"%s\t%d\t%d\n",read->qname,read->cigarLen, sizeof(read->cigarOps)/ sizeof(uint32_t));
            continue;
        }
        std::vector<uint32_t> cigar(read->cigarOps, read->cigarOps + 2*read->cigarLen );
        //report
//        fprintf(stderr,"%s\t%s\t%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d",read->qname,read->pos,read->end,bed[p1].Primer_ID,bed[p2].Primer_ID,)

        // if the alignment starts before the end of the primer, trim to that position
        uint32_t primer_position;
        if(ARGS_START){
            primer_position = bed[p1].start;
        } else{
            primer_position = bed[p1].end;
        }
//        fprintf(stderr,"%s\t%d\t%d\t%d\t%d\n",read->qname,read->cigarLen, read->pos,primer_position,read->end);
        int both = 0;
        if(read->pos < primer_position){
            int trim_success = trim(&cigar, read, primer_position, 0);
//            break;
            both = 1;
            if(trim_success == 0){
                continue;
            }
        }
        else{
            if(ARGS_VERBORSE == 1){
//                fprintf(stderr,"ref end %d >= primer_position %d\n",read->end,primer_position);
            }
        }

        if(ARGS_START){
            primer_position = bed[p2].end;
        } else{
            primer_position = bed[p2].start;
        }
        if(read->end > primer_position){
            trim(&cigar,read, primer_position, 1);
            both = 2;
        }
//        if(both == 2){
//            fprintf(stderr,"doing two trims for a single record?");
//        }
        else{
            if(ARGS_VERBORSE){
//                fprintf(stderr,"ref end %d >= primer_position %d\n",read->end,primer_position);
            }
        }

        // todo: exception handle
        if(ARGS_NORMALISE){
            char key[100];
            strcpy (key,bed[p1].Primer_ID);
            strcat(key,bed[p2].Primer_ID);
            int index = (read->flag & BAM_FREVERSE) ? 1:0;
            counter[key][index]++;
//            fprintf(stderr,"bitwise op done. key=%s index= %d counter=%d\n",key,index,counter[key][index]);
            if(counter[key][index] > ARGS_NORMALISE){
                continue;
            }
        }

        if(check_still_matching_bases(read)==0){
            continue;
        }

        // write back to b
        // we only have to update cigar in b with std::vector<> cigar
        uint32_t *cigar_ptr = bam_get_cigar(b);
        uint32_t *new_cigar_ptr = (uint32_t*)malloc(sizeof(uint32_t)*cigar.size());

        for(auto it = cigar.begin(); it != cigar.end(); ++it) {
            uint32_t flag = *it;
            ++it;
            uint32_t length = *it;
            length = length << BAM_CIGAR_SHIFT;
            *new_cigar_ptr = (length | flag);
//            fprintf(stderr,"flag = %d length = %d new value = %d\n",flag,*it,*new_cigar_ptr);
            new_cigar_ptr++;
        }
        b->core.n_cigar = cigar.size(); // cause error >> E::sam_format1_append] Corrupted aux data for read 20a4e805-b12b-410a-8cdf-d4d0046187dc

        ret=sam_write1(out,header,b);
        errorCheck(ret);
//        free(new_cigar_ptr);

    }

    fclose(fout);

    //wrap up
    free(read);
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
    sam_close(out);

    std::cout << "hello world" << std::endl;

    return 0;



}