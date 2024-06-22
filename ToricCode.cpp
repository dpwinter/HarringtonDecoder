#include "ToricCode.h"
#include "Location.h"

#include <iostream>


ToricCode::ToricCode(int L) {
    this->L = L;
    this->qubits = new bool**[L];
    this->stabs = new bool*[L];

    for (int i=0; i<L; i++) {
        this->qubits[i] = new bool*[L];
        this->stabs[i] = new bool[L];
        for (int j=0; j<L; j++) {
            this->qubits[i][j] = new bool[2]; // (N,W)-unit cell
        }
    }
    
    this->reset();
}

ToricCode::~ToricCode() {
    for (int i = 0; i<this->L; i++) {
        delete[] this->stabs[i];
        for (int j=0; j<this->L; j++) {
            delete[] this->qubits[i][j];
        }
    }
}

void ToricCode::reset() {
    for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
            this->stabs[i][j] = 0;
            this->qubits[i][j][0] = 0;
            this->qubits[i][j][1] = 0;
		}
	}
}

bool ToricCode::getStab(int i, int j) { // plaquette operator
    bool ret = 0;

    ret ^= this->qubits[i][j][0]; // N
    ret ^= this->qubits[i][j][1]; // W
    ret ^= this->qubits[i][(j+1)%this->L][1];  // E
    ret ^= this->qubits[(i+1)%this->L][j][0];  // S

    return ret;
}

bool** ToricCode::getSyndromes() {
    for (int i = 0; i<this->L; i++) {
        for (int j = 0; j<this->L; j++) {
            this->stabs[i][j] = this->getStab(i,j);
        }
    }
    return this->stabs;
}

bool ToricCode::getQubit(int i, int j, int k) {
    return this->qubits[i][j][k];
}

void ToricCode::noise(double p) {
    for (int i = 0; i<this->L; i++) {
        for (int j = 0; j<this->L; j++) {
            for (int k=0; k<2; k++) {
                double r = this->randDist(this->randGen);
                if(r <= p) {
                    this->qubits[i][j][k] ^= 1;
                }
            }
        }
    }
}

void ToricCode::flip(int i, int j, int loc) {

    switch(loc) {
        case Location::N: 
            this->qubits[i][j][0] ^= 1;
            break;
        case Location::W:
            this->qubits[i][j][1] ^= 1;
            break;
        case Location::E:
            this->qubits[i][(j+1)%this->L][1] ^= 1;
            break;
        case Location::S:
            this->qubits[(i+1)%this->L][j][0] ^= 1;
            break;
    }
}

bool ToricCode::hasLogErr() {
    int paritySumRows = 0;
    int paritySumCols = 0;

    for(int i=0; i<this->L; i++) {
        bool parityRows = 0;
        bool parityCols = 0;

        for(int j=0; j<this->L; j++) {
            parityRows ^= this->qubits[i][j][0];
            parityCols ^= this->qubits[j][i][1];
        }
        
        paritySumRows += int(parityRows);
        paritySumCols += int(parityCols);
    }
    
    return (paritySumRows > this->L/2) || (paritySumCols > this->L/2);
}

