#ifndef MEMORY_H_
#define MEMORY_H_

#include "Location.h"

struct Memory {

    Location addr;
    int U; // work period (of this level)
    int Q; // colony size
    int age = 0; // time step % U

    bool countSig[8] = {0,0,0,0,0,0,0,0}; // N,W,E,S,NW,NE,SW,SE
    bool n_countSig[8] = {0,0,0,0,0,0,0,0};
    bool flipSig[4] = {0,0,0,0}; // N,W,E,S
    bool n_flipSig[4] = {0,0,0,0};
    int count[9] = {0,0,0,0,0,0,0,0,0}; // N,W,E,S,NW,NE,SW,SE,C
};

#endif
