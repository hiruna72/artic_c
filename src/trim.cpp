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

int find_primer(std::vector<bed_row> bed, int pos, int direction) {
    int min_distance = INT32_MAX;
    size_t index = 0;
    int i = 0;
    for(auto it=bed.begin();it!=bed.end();++it){
        if(it->direction_sign == direction && direction == 1){
            if(abs(it->start-pos)<min_distance){
                min_distance = abs(it->start-pos);
                index = i;
            }
        }
        if(it->direction_sign == direction && direction == -1){
            if(abs(it->end-pos)<min_distance){
                min_distance = abs(it->end-pos);
                index = i;
            }
        }
        i++;
    }
    return index;
}

int trim(alignedRead *read, int start_pos, int end) {
    int pos;
    if(end==0){
        pos = read->pos;
    } else{
        pos = read->end;
    }
    int eaten = 0;
    int index = -1;
    if(end == 1){
        index = read->cigarLen;
    }
    while(1){
        uint32_t flag;
        uint32_t length;
        // chomp stuff until we reach pos
        if(end == 1){
            index--;
        }else{
            index++;
        }
        flag = read->cigarOps[2*index];
        length = read->cigarOps[2*index + 1];
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

    int extra = abs(pos-start_pos);
    

    return 0;
}

int main(int argc, char ** argv){
    // args flags
    int REMOVE_INCORRECT_PAIRS = 1;
    int START = 1;

    char * samfile = "samtools.bam";
    std::vector<bed_row> bed = read_bed_file("nCoV-2019.bed");

    //these come from htslib/sam.h
    samFile *in = NULL;
    bam1_t *b= NULL;
    bam_hdr_t *header = NULL;

    //open the BAM file (though called sam_open is opens bam files too :P)
    in = sam_open(samfile, "r");
    if(in == nullptr){
        fprintf(stderr,"could not open file %s",samfile);
    }

    //get the sam header.
    if ((header = sam_hdr_read(in)) == 0){
        fprintf(stderr,"No sam header?\n");
        exit(EXIT_FAILURE);
    }

    //print the chromosome names in the header
    //see the bam_hdr_t struct in htslib/sam.h for parsing the rest of the stuff in header
    int i;
    for(i=0; i< (header->n_targets); i++){
        printf("Chromosome ID %d = %s\n",i,(header->target_name[i]));
    }

    b = bam_init1();

    //my structure for a read (see common.h)
    struct alignedRead* read = (struct alignedRead*)malloc(sizeof(struct alignedRead));

    //repeat until all reads in the file are retrieved
    while ( sam_read1(in, header, b) >= 0){
        getRead(read, b);         //copy the current read to the myread structure. See common.c for information
        if(read->flag & BAM_FUNMAP == 0){
            fprintf(stderr,"%s read is skipped because it is unmapped",read->qname);
            continue;
        }
        if(read->flag & BAM_FSUPPLEMENTARY == 0){
            fprintf(stderr,"%s read is skipped because it is supplementary",read->qname);
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

        if(REMOVE_INCORRECT_PAIRS==1 and is_correctly_paired!=0){
            continue;
        }
        //report
//        fprintf(stderr,"%s\t%s\t%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d",read->qname,read->pos,read->end,bed[p1].Primer_ID,bed[p2].Primer_ID,)

        // if the alignment starts before the end of the primer, trim to that position
        int primer_position;
        if(START){
            primer_position = bed[p1].start;
        } else{
            primer_position = bed[p1].end;
        }

        if(read->pos < primer_position){
            int trim_success = trim(read,primer_position,0);
            if(!trim_success){
                continue;
            }
        }


//        printRead(myread,header);   //print data in  myread structure. See common.c for information
    }


    //wrap up
    free(myread);
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);


    return 0;

//    std::cout << "hello world" << std::endl;


}