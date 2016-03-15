#ifndef _HAM_HOG_H_
#define _HAM_HOG_H_

#include <stdio.h>
typedef unsigned char uchar;

typedef struct {
    int cols;
    int rows;
    uchar *d;
}Mat;

typedef struct {
    size_t gradOfs;
    size_t qangleOfs;
    int histOfs[4];
    float histWeights[4];
    float gradWeight;
}PixData;

typedef struct {
    int histOfs;
    int imgOfsX;
    int imgOfsY;
}BlockData;



#endif

