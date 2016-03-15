#include <string.h>
#include <stdlib.h>
#include "hog.h"

int main(int argc, char **argv)
{
    FILE *fp = NULL;
    uchar *pic_buf = NULL;
    Mat lmat;

    if (argc < 2)
    {
        printf("Please type:%s <PictureFile>\n", argv[0]);
        return -1;
    }
    lmat.rows = 480;
    lmat.cols = 640;
    // get image buffer
    fp = fopen(argv[1], "r");
    lmat.d = (uchar *)malloc(640*480*2*sizeof(uchar));
    fread(pic_buf, 640*480*2,1, fp);
    fclose(fp);

    //
    printf("test point0\n");
    hog_init();
    printf("test point1\n");
    hog_prepare(&lmat);
    printf("test point2\n");

    hog_detect();
    printf("test point3\n");
    hog_destroy();
    printf("test point4\n");
    // free image buffer
    free(pic_buf);

    return 0;
}
