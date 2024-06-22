#ifndef CA_H_
#define CA_H_

#include "Cell.h"
#include "Location.h"
#include "ToricCode.h"

class CA {

    private:
        int L;
        Cell*** cells;
        Location** corrections;

    public:
        CA(int L, int U, double fC, double fN);
        virtual ~CA();
        void reset();
        Cell* getCell(int i, int j);
        Location** step(bool** syndromes);

};

#endif
