/*
 *  The Mx_xxx_Matrix class is based on mx.cpp from Robert Edgar.
 *
 *  Created by Sarah Westcott on 3/25/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "mxmatrix.h"

/******************************************************************************/
MxFloatMatrix::MxFloatMatrix(unsigned row, unsigned col, float init) {
    
    unsigned N = max(row, col);
	N += 16;

    matrix.resize(N);

    for(unsigned i = 0; i < N; i++) {
        matrix[i].resize(N, init);
    }
}
/******************************************************************************/
void MxFloatMatrix::resize(unsigned row, unsigned col) {
    
    unsigned N = max(row, col);
	N += 16;

    // resize neccesary
    if (matrix.size() < N) {
        matrix.resize(N);

        for(unsigned i = 0; i < N; i++) {
            matrix[i].resize(N, 0.0);
        }
    }
}
/******************************************************************************/
MxFloatMatrix::~MxFloatMatrix(){}
/******************************************************************************/
MxByteMatrix::MxByteMatrix(unsigned row, unsigned col, Byte init) {
    
    unsigned N = max(row, col);
	N += 16;

    matrix.resize(N);

    for(unsigned i = 0; i < N; i++) {
        matrix[i].resize(N, init);
    }
}
/******************************************************************************/
void MxByteMatrix::resize(unsigned row, unsigned col) {
    
    unsigned N = max(row, col);
	N += 16;

    // resize neccesary
    if (matrix.size() < N) {
        matrix.resize(N);

        for(unsigned i = 0; i < N; i++) {
            matrix[i].resize(N, 0);
        }
    }
}
/******************************************************************************/
MxByteMatrix::~MxByteMatrix(){}
/******************************************************************************/
