//
// Created by shan on 2020-04-02.
//

#ifndef ARTIC_C_ARTIC_H
#define ARTIC_C_ARTIC_H

// bed file values
typedef struct {
    char* chromosome;
    int start;
    int end;
    char* Primer_ID;
    char* PoolName;
    char* direction;
    int direction_sign;
    int original;
} bed_row;

#endif //ARTIC_C_ARTIC_H
