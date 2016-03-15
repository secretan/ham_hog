#include "hog.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

// Window paramenter
int window_width = 64;
int window_height = 128;
int window_stride = 8;

// Block paramenter
int block_width = 16;
int block_height = 16;
int block_stride = 8;

// Cell paramenter
int cell_width = 8;
int cell_height = 8;

//
PixData *pixData = NULL;
BlockData *blockData = NULL;

//
Mat *grad = NULL;
Mat *angle = NULL;

void hog_init()
{
    int w = 0, h = 0;

    // compute weight start
    float *weight = (float *)malloc(block_width*block_height*sizeof(float));
    float win_sigma = (block_width+block_height)/8;
    float scale = 1.0/(win_sigma*win_sigma*2);
    {
        float *dw = (float *)malloc(block_width*sizeof(float));
        float *dh = (float *)malloc(block_height*sizeof(float));

        for (w = 0; w < block_width; w++)
        {
            dw[w] = dw[w]*(w - block_width/2);
        }
        for (h = 0; h < block_height; h++)
        {
            dh[h] = dh[w]*(h - block_height/2);
        }

        for (h = 0; h < block_height; h++)
        {
            for (w = 0; w < block_width; w++)
            {
                weight[h*block_width+w] = exp(-(dw[w]+dh[h])*scale);
            } // for (w = 0; w < block_width; w++)
        }   // for (h = 0; h < block_height; h++)
    }
    // compute weight end

    // compute Pixel Data start
    int rawBlockSize = block_height*block_width;
    pixData = (PixData*)malloc(block_height*block_width*sizeof(PixData)*3);
    memset(pixData, 0, rawBlockSize*sizeof(PixData)*3);
    int count1 = 0,count2 = 0,count4=0;
    int i,j;
    const char ncells_height = 2;
    const char ncells_width = 2;
    const char nbins = 9;
    for (h = 0; h < block_height; h++)
    {
        for (w = 0; w < block_width; w++)
        {
            PixData* data = 0;
            i = h;j = w;
            float cellX = (j+0.5f)/cell_width - 0.5f;
            float cellY = (i+0.5f)/cell_height - 0.5f;
            int icellX0 = (int)floor(cellX);
            int icellY0 = (int)floor(cellY);
            int icellX1 = icellX0 + 1, icellY1 = icellY0 + 1;
            cellX -= icellX0;
            cellY -= icellY0;

            if( (unsigned)icellX0 < (unsigned)ncells_width &&
               (unsigned)icellX1 < (unsigned)ncells_width )
            {
                if( (unsigned)icellY0 < (unsigned)ncells_height &&
                   (unsigned)icellY1 < (unsigned)ncells_height )
                {
                    data = &pixData[rawBlockSize*2 + (count4++)];
                    data->histOfs[0] = (icellX0*ncells_height + icellY0)*nbins;
                    data->histWeights[0] = (1.f - cellX)*(1.f - cellY);
                    data->histOfs[1] = (icellX1*ncells_height + icellY0)*nbins;
                    data->histWeights[1] = cellX*(1.f - cellY);
                    data->histOfs[2] = (icellX0*ncells_height + icellY1)*nbins;
                    data->histWeights[2] = (1.f - cellX)*cellY;
                    data->histOfs[3] = (icellX1*ncells_height + icellY1)*nbins;
                    data->histWeights[3] = cellX*cellY;
                }
                else
                {
                    data = &pixData[rawBlockSize + (count2++)];
                    if( (unsigned)icellY0 < (unsigned)ncells_height )
                    {
                        icellY1 = icellY0;
                        cellY = 1.f - cellY;
                    }
                    data->histOfs[0] = (icellX0*ncells_height + icellY1)*nbins;
                    data->histWeights[0] = (1.f - cellX)*cellY;
                    data->histOfs[1] = (icellX1*ncells_height + icellY1)*nbins;
                    data->histWeights[1] = cellX*cellY;
                    data->histOfs[2] = data->histOfs[3] = 0;
                    data->histWeights[2] = data->histWeights[3] = 0;
                }
            }
            else
            {
                if( (unsigned)icellX0 < (unsigned)ncells_width )
                {
                    icellX1 = icellX0;
                    cellX = 1.f - cellX;
                }

                if( (unsigned)icellY0 < (unsigned)ncells_height &&
                   (unsigned)icellY1 < (unsigned)ncells_height )
                {
                    data = &pixData[rawBlockSize + (count2++)];
                    data->histOfs[0] = (icellX1*ncells_height + icellY0)*nbins;
                    data->histWeights[0] = cellX*(1.f - cellY);
                    data->histOfs[1] = (icellX1*ncells_height + icellY1)*nbins;
                    data->histWeights[1] = cellX*cellY;
                    data->histOfs[2] = data->histOfs[3] = 0;
                    data->histWeights[2] = data->histWeights[3] = 0;
                }
                else
                {
                    data = &pixData[count1++];
                    if( (unsigned)icellY0 < (unsigned)ncells_height )
                    {
                        icellY1 = icellY0;
                        cellY = 1.f - cellY;
                    }
                    data->histOfs[0] = (icellX1*ncells_height + icellY1)*nbins;
                    data->histWeights[0] = cellX*cellY;
                    data->histOfs[1] = data->histOfs[2] = data->histOfs[3] = 0;
                    data->histWeights[1] = data->histWeights[2] = data->histWeights[3] = 0;
                }
            }
            data->gradOfs = (block_width*i + j);
            data->qangleOfs = (block_width*i + j);
            data->gradWeight = weight[i*block_width+j];
        } // for (w = 0; w < block_width; w++)
    }   // for (h = 0; h < block_height; h++)
#if 0
    if (count1+count2+count4 == block_width*block_height)
    {
        for (i = 0; i < count2; i++)
        {
            //pixData[i+count1] = pixData[i+block_width*block_height];
        }
        for (i = 0; i < count4; i++)
        {
            //pixData[i+count2+count1] = pixData[i+rawBlockSize*2];
        }
        count2 += count1;
        count4 += count2;
    }
#endif
    // compute Pixel Data end

    // compute Block Data Offset start
    char nblockw = (window_width-block_width)/block_stride + 1;
    char nblockh = (window_height-block_height)/block_stride + 1;
    blockData = (BlockData *)malloc(nblockw*nblockh*sizeof(BlockData));
    memset(blockData, 0, nblockw*nblockh*sizeof(BlockData));
    for (h = 0; h < nblockh; h++)
    {
        for (w = 0; w < nblockw; w++)
        {
            blockData[h*nblockw+w].histOfs = (h*nblockw+w)*4*9;
            blockData[h*nblockw+w].imgOfsX = w*block_stride;
            blockData[h*nblockw+w].imgOfsY = h*block_stride;
        }
    }
    // compute Block Data Offset end
    free(weight);
}

static void compute_gradient(Mat *_img, Mat *_grad, Mat *_angle)
{
    int i;
    int h, w;

    // initilize gradient & angle
    _grad->cols = _img->cols*2;
    _grad->rows = _img->rows;
    _grad->d = (uchar *)malloc(_grad->cols*_grad->rows*sizeof(uchar));
    memset(_grad->d, 0, _grad->cols*_grad->rows*sizeof(uchar));
    _angle->cols = _img->cols*2;
    _angle->rows = _img->rows;
    _angle->d = (uchar *)malloc(_angle->cols*_angle->rows*sizeof(uchar));
    memset(_angle->d, 0, _angle->cols*_angle->rows*sizeof(uchar));

    // gamma transfer
    char gamma_trans_flag = 1;
    float *_lut = (float *)malloc(256*sizeof(float));
    if (gamma_trans_flag)
    {
        for (i = 0; i < 256; i++)
        {
            _lut[i] = sqrt((float)i);
        }
    }
    else
    {
        for (i = 0; i < 256; i++)
        {
            _lut[i] = (float)i;
        }
    }
    // computing
    float angleScale = 9.0/(3.1415);
    float *dbuf = (float*)malloc(_img->cols*4*sizeof(float));
    int lwidth = _img->cols;
    float mag = 0.0;
    float tangle = 0.0;
    int angleOfs = 0;


    for (h = 1; h < _img->rows-1; h++)
    {
        uchar *prevh = &_img->d[(h-1)*_img->cols];
        uchar *nexth = &_img->d[(h+1)*_img->cols];
        // compute one line start
        for (w = 1; w < _img->cols-1; w++)
        {
            //
            dbuf[w] = _lut[_img->d[h*lwidth+w+1]]-_lut[_img->d[h*lwidth+w-1]];
            dbuf[w+lwidth] = _lut[nexth[w]] - _lut[prevh[w]];
            //
            mag = (float)sqrt(dbuf[w]*dbuf[w]+ dbuf[w+lwidth]*dbuf[w+lwidth]);
            tangle = (float)atanf(dbuf[w+lwidth]/(dbuf[w+lwidth]+0.001));
            tangle = -0.5+tangle*angleScale;
            angleOfs = (int)floorf(tangle);
            tangle -= angleOfs;
            //
            _grad->d[(h*lwidth+w)*2+0] = mag*(1-tangle);
            _grad->d[(h*lwidth+w)*2+1] = mag*tangle;
            //
            if (angleOfs < 0)
            {
                angleOfs += 9;
            }
            else if (angleOfs > 9)
            {
                angleOfs -= 9;
            }
            _angle->d[(h*lwidth+w)*2+0] = angleOfs;
            angleOfs += 1;
            if (angleOfs >= 9)
            {
                _angle->d[(h*lwidth+w)*2+1] = 0;
            }
            else
            {
                _angle->d[(h*lwidth+w)*2+1] = angleOfs;
            }
        }
        // compute one line end
    }

    free(dbuf);
}

void hog_prepare(Mat *_mat)
{
    //
    grad = (Mat *)malloc(sizeof(Mat));
    angle = (Mat *)malloc(sizeof(Mat));
    // compute gradient
    compute_gradient(_mat, grad, angle);
}

void hog_detect()
{

}

void hog_destroy()
{
    free(grad->d);
    free(grad);
    free(angle->d);
    free(angle);
    free(pixData);
    free(blockData);
}

