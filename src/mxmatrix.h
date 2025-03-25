#ifndef MXMATRIX_H
#define MXMATRIX_H

/*
 *  The Mx_xxx_Matrix class is based on mx.h from Robert Edgar.
 *
 *  Created by Sarah Westcott on 3/25/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 */

#include "uchime.h"
#include "seqdb.h"

const float MINUS_INFINITY = -9e9f;

/******************************************************************************/
class MxFloatMatrix {
    public:

        MxFloatMatrix(unsigned row = 0, unsigned col = 0, float init = 0.0);
        ~MxFloatMatrix();

        void resize(unsigned row = 0, unsigned col = 0);

        vector<vector<float> > matrix;

    private:

        

};
/******************************************************************************/

class MxByteMatrix {
    public:

        MxByteMatrix(unsigned RowCount = 0, unsigned ColCount = 0, Byte init = 0);
        ~MxByteMatrix();

        void resize(unsigned RowCount = 0, unsigned ColCount = 0);

        vector<vector<Byte> > matrix;

    private:

        

};
/******************************************************************************/

#endif
