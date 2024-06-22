#ifndef CELL_H_
#define CELL_H_

#include "Memory.h"
#include "Location.h"

class Cell {
    private:
        Location addr; // level-0 address (within k=0 hierarchy)
        Cell** neighbors; // 8 nearest neighbor cells
        Memory** memory; // one memory per hierarchy level
        bool* syndromes; // anyon presence (N,W,E,S,NW,NE,SW,SE,C)

        int d; // max. hierarchy level
        double fN; // threshold for count of neighbor signals
        double fC; // threshold for count of own syndrome
    public:
        Cell(int row, int col, int Q, int U, int d, double fC, double fN);
        virtual ~Cell();
        void reset();

        void acquire(); // Get data from neighbors
        void update(); // move signal data from temp to actual (or broadcast)
        Location rule(); // apply local rule to actual data (or higher-level rule)

        void setNeighbors(Cell**); // assign neighbor cells
        void setSyndrome(bool syndrome); // set current center syndrome (i.e. anyon presence)
        Memory* getMemory(int k); // get k-th level memory of this cell

        Location harringtonRule(Location addr, bool* syndromes);
};

#endif
