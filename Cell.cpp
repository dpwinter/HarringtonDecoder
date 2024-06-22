#include "Cell.h"
#include "Location.h"
#include "Memory.h"

#include <cmath>

Cell::Cell(int row, int col, int Q, int U, int d, double fC, double fN) {

    this->d = d;
	this->fC = fC;
	this->fN = fN;

    this->neighbors = new Cell*[8];
    this->syndromes = new bool[9];
    this->addr = locFromCoords(row, col);

	this->memory = new Memory*[d-1];
    for(int k=0; k<d; k++) {

        int offset = int((std::pow(Q,k) - 1) / 2.0);
        float krow = (row - offset) / std::pow(Q,k);
        float kcol = (col - offset) / std::pow(Q,k);
        Location kaddr = (int(krow) == krow && int(kcol) == kcol) ? locFromCoords(int(krow) % Q, (int(kcol) % Q)) : Location::None;
    
        if(k==0) {
            this->addr = kaddr; // assign level-0 address
        } else {
            this->memory[k-1] = new Memory {kaddr, int(std::pow(U,k)), int(std::pow(Q,k))};
        }
    }

    this->reset();
}

Cell::~Cell() {
    delete[] this->neighbors;
    delete[] this->syndromes;

	for(int i=0; i<this->d-1; i++) {
        delete[] this->memory[i];
    }
}

void Cell::reset() {
    for (int i=0; i<9; i++) {
        this->syndromes[i] = 0;
    }

	for (int k=0; k<this->d-1; k++) {
        this->memory[k]->age = 0;

        for (int i=0; i<4; i++) { 
            this->memory[k]->flipSig[i] = 0;
            this->memory[k]->n_flipSig[i] = 0;
        }

        for (int i=0; i<8; i++) {
            this->memory[k]->countSig[i] = 0;
            this->memory[k]->n_countSig[i] = 0;
        }

        for (int i=0; i<9; i++) {
            this->memory[k]->count[i] = 0;
        }
        
    }
}

void Cell::setNeighbors(Cell** neighbors) {
    this->neighbors = neighbors;
}

void Cell::setSyndrome(bool syndrome) {
    this->syndromes[Location::C] = syndrome;
}

Memory* Cell::getMemory(int k) {
    return this->memory[k];
}

void Cell::acquire() {

    // copy neighbor center syndrome
    for (int i=0; i<8; i++) {
        this->syndromes[i] = this->neighbors[i]->syndromes[Location::C];
    }

	// propagate signals (from opposite neighbor in same direction)
    for (int k=0; k<this->d-1; k++) {
        for (int i=0; i<8; i++)
            this->memory[k]->n_countSig[i] = this->neighbors[oppositeLoc(Location(i))]->memory[k]->countSig[i];
        for (int i=0; i<4; i++)
            this->memory[k]->n_flipSig[i] = this->neighbors[oppositeLoc(Location(i))]->memory[k]->flipSig[i];

    }
}

void Cell::update() {
    for (int k=0; k<this->d-1; k++) {
        this->memory[k]->age = (this->memory[k]->age + 1) % this->memory[k]->U; // increment age

        if (this->memory[k]->addr != Location::None) { // hierarchy representatives
            for (int i=0; i<8; i++)
                this->memory[k]->countSig[i] = this->syndromes[Location::C]; // broadcast

			// update count array
			this->memory[k]->count[Location::C] += this->syndromes[Location::C];
			for (int i=0; i<8; i++)
				this->memory[k]->count[i] += this->memory[k]->n_countSig[oppositeLoc(Location(i))]; // direction it came from

        } else { // copy signals for all non-representatives
            for (int i=0; i<8; i++)
                this->memory[k]->countSig[i] = this->memory[k]->n_countSig[i];
            for (int i=0; i<4; i++)
                this->memory[k]->flipSig[i] = this->memory[k]->n_flipSig[i];
        }
    }
}

Location Cell::rule() {
	for (int k=0; k<this->d-1; k++) {
        if (this->memory[k]->age == 0 && this->memory[k]->addr != Location::None) { // at t=U -> decide flipSig

            bool* syndromes = new bool[9]; // k-level syndrome..
            for (int i=0; i<9; i++) {
                double f = (i == Location::C) ? this->fC : this->fN;
                syndromes[i] = this->memory[k]->count[i] >= f * this->memory[k]->U; // ..determined from k-level count
                this->memory[k]->count[i] = 0; // reset count
            }

            Location dir = this->harringtonRule(this->memory[k]->addr, syndromes); // higher-level rule
            if (dir != Location::None) // emit flipSig
				this->memory[k]->flipSig[dir] = 1;
        }
        else if (this->memory[k]->age == this->memory[k]->Q) { // at t=U+Q, do correction chain, if applicable

            for (int i=0; i<4; i++) {
                if (this->memory[k]->flipSig[i]) {
					this->memory[k]->flipSig[i] = 0;
					return Location(i); // issue correction in direction of (first) flipSig
                } 
            }
        }
    }
    return this->harringtonRule(this->addr, this->syndromes);
}

Location Cell::harringtonRule(Location addr, bool* syndromes) {

    // Immediate exit
    if(addr == Location::C || syndromes[Location::C] == 0) {
        return Location::None;
    }
    // W border
    if(addr == Location::NW || addr == Location::W || addr == Location::SW) {
        if(syndromes[Location::W] || syndromes[Location::NW] || syndromes[Location::SW]) { 
            return Location::W; 
        }
    }
    // S border
    if(addr == Location::S || addr == Location::SW || addr == Location::SE) {
        if(syndromes[Location::S] || syndromes[Location::SW] || syndromes[Location::SE]) { 
            return Location::S; 
        }
    }

    // SW quadrant
    if(addr == Location::SW) {
        if(syndromes[Location::S] || syndromes[Location::W]) {
            return Location::None;
        }
        else if(syndromes[Location::N]) {
            return Location::N;
        }
        else if(syndromes[Location::E]) {
            return Location::E;
        }
        else if(syndromes[Location::SW]) {
            return Location::None;
        }
        else if(syndromes[Location::NW]) {
            return Location::N;
        }
        else if(syndromes[Location::SE]) {
            return Location::E;
        }
        else {
            return Location::E;
        }
    }

    // W corridor
    if(addr == Location::W) {
        if(syndromes[Location::S] || syndromes[Location::W] || syndromes[Location::N]) {
            return Location::None;
        }
        else if(syndromes[Location::E]) {
            return Location::E;
        }
        else if(syndromes[Location::SW] || syndromes[Location::NW]) {
            return Location::None;
        }
        else {
            return Location::E;
        }
    }

    // NW quadrant
    if(addr == Location::NW) {
        if(syndromes[Location::W] || syndromes[Location::N]) {
            return Location::None;
        }
        else if(syndromes[Location::E]) {
            return Location::E;
        }
        else if(syndromes[Location::S]) {
            return Location::S;
        }
        else if(syndromes[Location::NW]) {
            return Location::None;
        }
        else if(syndromes[Location::NE]) {
            return Location::E;
        }
        else if(syndromes[Location::SW]) {
            return Location::S;
        }
        else {
            return Location::E;
        }
    }

    // N corridor
    if(addr == Location::N) {
        if(syndromes[Location::W] || syndromes[Location::N] || syndromes[Location::E]) {
            return Location::None;
        }
        else if(syndromes[Location::S]) {
            return Location::S;
        }
        else if(syndromes[Location::NW] || syndromes[Location::NE]) {
            return Location::None;
        }
        else {
            return Location::S;
        }
    }

    // NE quadrant
    if(addr == Location::NE) {
        if(syndromes[Location::N] || syndromes[Location::E]) {
            return Location::None;
        }
        else if(syndromes[Location::S]) {
            return Location::S;
        }
        else if(syndromes[Location::W]) {
            return Location::W;
        }
        else if(syndromes[Location::NE]) {
            return Location::None;
        }
        else if(syndromes[Location::SE]) {
            return Location::S;
        }
        else if(syndromes[Location::NW]) {
            return Location::W;
        }
        else {
            return Location::W;
        }
    }

    // E corridor
    if(addr == Location::E) {
        if(syndromes[Location::N] || syndromes[Location::E] || syndromes[Location::S]) {
            return Location::None;
        }
        else if(syndromes[Location::W]) {
            return Location::W;
        }
        else if(syndromes[Location::NE] || syndromes[Location::SE]) {
            return Location::None;
        }
        else {
            return Location::W;
        }
    }

    // SE quadrant
    if(addr == Location::SE) {

        if(syndromes[Location::E] || syndromes[Location::S]) {
            return Location::None;
        }
        else if(syndromes[Location::W]) {
            return Location::W;
        }
        else if(syndromes[Location::N]) {
            return Location::N;
        }
        else if(syndromes[Location::SE]) {
            return Location::None;
        }
        else if(syndromes[Location::SW]) {
            return Location::W;
        }
        else if(syndromes[Location::NE]) {
            return Location::N;
        }
        else {
            return Location::W;
        }
    }

    // S corridor
    if(addr == Location::S) {
        if(syndromes[Location::E] || syndromes[Location::S] || syndromes[Location::W]) {
            return Location::None;
        }
        else if(syndromes[Location::N]) {
            return Location::N;
        }
        else if(syndromes[Location::SE] || syndromes[Location::SW]) {
            return Location::None;
        }
        else {
            return Location::N;
        }
    }

    return Location::None;
}
