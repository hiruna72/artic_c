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

int main(int argc, char ** argv){
//    std::vector<bed_row> bed_file = read_bed_file("nCoV-2019.bed");

    samFile *fp_in = hts_open("samtools.bam","r"); //open bam file
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment

    char *chrom = "chr22";
    int locus = 0;
    int comp ;

    printf("%s\t%d\n", chrom, locus);

    //header parse
    //uint32_t *tar = bamHdr->text ;
    //uint32_t *tarlen = bamHdr->target_len ;

    //printf("%d\n",tar);

    while(sam_read1(fp_in,bamHdr,aln) > 0){

        int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
        char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
        uint32_t len = aln->core.l_qseq; //length of the read.

        uint8_t *q = bam_get_seq(aln); //quality string
        uint32_t q2 = aln->core.qual ; //mapping quality


        char *qseq = (char *)malloc(len);

        for(int i=0; i< len ; i++){
            qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
        }

        //printf("%s\t%d\t%d\t%s\t%s\t%d\n",chr,pos,len,qseq,q,q2);

        if(strcmp(chrom, chr) == 0){

            if(locus > pos+len){
                printf("%s\t%d\t%d\t%s\t%s\t%d\n",chr,pos,len,qseq,q,q2);
            }
        }
    }

    bam_destroy1(aln);
    sam_close(fp_in);

    return 0;

    std::cout << "hello world" << std::endl;


}