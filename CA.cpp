#include "CA.h"
#include "Cell.h"
#include "Location.h"

#include <cmath>
#include <cassert>

CA::CA(int L, int U, double fC, double fN) {
    this->L = L; // linear size of lattice
    int Q = 3; // colony size (hard-coded)

    float df = std::log2(L)/std::log2(Q); // calc. hierarchy level
    assert( (::ceilf(df) == df) || (::floorf(df) == df) ); // assure L and Q are compatible
    int d = int(df); // hierarchy level

    // create cells
    this->cells = new Cell**[L];
    for (int i=0; i<L; i++) {
        this->cells[i] = new Cell*[L];
        for (int j=0; j<L; j++) {
            this->cells[i][j] = new Cell(i,j,Q,U,d,fC,fN);
        }
    }

    // output of global rule: LxL corrections
    this->corrections = new Location*[L];
    for (int i=0; i<L; i++) {
        this->corrections[i] = new Location[L];
    }

    // assign neighbors
    for (int i=0; i<L; i++) {
        for (int j=0; j<L; j++) {
            Cell** neighbors = new Cell*[8];

            int i_ = (i-1 == -1) ? L-1 : i-1;
            int j_ = (j-1 == -1) ? L-1 : j-1;
            int ip = (i+1) % L;
            int jp = (j+1) % L;

            neighbors[Location::N]  = this->cells[i_][j];
            neighbors[Location::W]  = this->cells[i][j_];
            neighbors[Location::E]  = this->cells[i][jp];
            neighbors[Location::S]  = this->cells[ip][j];
            neighbors[Location::NW] = this->cells[i_][j_];
            neighbors[Location::NE] = this->cells[i_][jp];
            neighbors[Location::SW] = this->cells[ip][j_];
            neighbors[Location::SE] = this->cells[ip][jp];

            this->cells[i][j]->setNeighbors(neighbors);
        }
    }
}

CA::~CA() {
    for(int i=0; i<this->L; i++) {
        delete[] this->cells[i];
    }
}

void CA::reset() {
    for(int i=0; i<this->L; i++) {
        for(int j=0; j<this->L; j++) {
            this->corrections[i][j] = Location::None;
            this->cells[i][j]->reset();
        }
    }
}

Cell* CA::getCell(int i, int j) {
    return this->cells[i][j];
}

Location** CA::step(bool** syndromes) {

    // 1. Measure syndrome, assign to cells
    for (int i=0; i<this->L; i++) {
        for (int j=0; j<this->L; j++) {
            this->cells[i][j]->setSyndrome(syndromes[i][j]);
        }
    }

    // 2. Copy neighbor data
    for (int i=0; i<this->L; i++) {
        for (int j=0; j<this->L; j++) {
            this->cells[i][j]->acquire();
        }
    }

    // 3. Synchronous update: temp->actual
    for (int i=0; i<this->L; i++) {
        for (int j=0; j<this->L; j++) {
            this->cells[i][j]->update();
        }
    }

    // 4. Perform (synchronized) local rule
    for (int i=0; i<this->L; i++) {
        for (int j=0; j<this->L; j++) {
            this->corrections[i][j] = this->cells[i][j]->rule();
        }
    }

    return this->corrections;

}
