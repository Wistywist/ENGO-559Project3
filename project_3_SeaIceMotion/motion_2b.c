/* motion20.c */
/* Last modified 6/19/96 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define NARGS 		5
#define NPARMS 		12
#define TABLESIZE 	1024
#define COR_THRESHOLD 	0.6
#define BUFFER_LEN	1024
#define STRING_LEN	64

#define round(x)	((int)(((x) < 0.) ? (x) - 0.5 : (x) + 0.5))
#define pixel(image,xsize,ysize,i,j) (((i) < 0 || (i) >= (xsize) || \
			             (j) < 0 || (j) >= (ysize)) ? 0 : \
				     (image)[j][i])

/* For debugging */
#define printd_int(var) fprintf(stderr, #var ": %d\t", var)
#define printd_double(var) fprintf(stderr, #var ": %f\t", var)
#define newline() fprintf(stderr, "\n")
#define got_here(i) fprintf(stderr, "Got here %d.\n", i);

typedef enum _ImageType
{
  int1b,
  int2bh,
  int2bl,
} ImageType;

void usage();
void find_target(short **image1, short **image2, int xsize, int ysize,
		 int i_0, int j_0, int box_width, int box_height, 
		 int skip_val, int sub_val,
                 int abs_range, int skip_range, int norm_range, int sub_range,
		 double *i_offset, double *j_offset, double *coor_coeff);
void ft_skip(short **image1, short **image2, int i_0, int j_0,
	     int box_width, int box_height, int skip_val, int skip_range,
	     double *i_offset, double *j_offset, double *coor_coeff);
void ft_norm(short **image1, short **image2, 
	     int i_0, int j_0, int i_start, int j_start,
	     int box_width, int box_height, int norm_range, int abs_range,
	     double *i_offset, double *j_offset, double *coor_coeff);
void ft_sub(short **image1, short **image2, int xsize, int ysize, 
	    int i_0, int j_0, int i_start, int j_start, 
	    int box_width, int box_height,
	    int sub_val, int sub_range, int abs_range,
	    double *i_offset, double *j_offset, double *coor_coeff);
double correlation(short **image1, short **image2, int i1_0, int j1_0, 
		   int i2_0, int j2_0, int box_width, int box_height);
int good_target(short **image, int i0, int j0, int box_width, int box_height);
short **expand(short **image, int xsize, int ysize, 
	       int i0, int j0, int in, int jn, int factor);
void freemem_2d(void *array, int jn);
void getparms(FILE *parmfile, char *parmfilename, int *xsize, int *ysize, 
	      int *box_width, int *box_height, int *skip_val, int *sub_val, 
	      int* abs_range, int *skip_range, int *sub_range, int *norm_range,
	      int *overlap, int *mask_min, int *mask_max, 
	      double *delta_t, double *resolution, 
	      ImageType *image_type, int *buffer_output, 
	      char **map_filename, char **points_filename);
FILE *open_read(char *filename);
FILE *open_write(char *filename);
void *getmem(int nbytes);
void *getmem_2d(int nrows, int nbytes);
void readimage_uint1b(short **image, int xsize, int ysize, FILE *file, 
		      char *filename);
void readimage_int2bh(short **image, int xsize, int ysize, FILE *file,
		      char *filename);
void readimage_int2bl(short **image, int xsize, int ysize, FILE *file,
		      char *filename);
void write_image(short **image, int xsize, int ysize, char *filename);
void readmap(unsigned char **map, int xsize, int ysize, FILE *file,
		      char *filename);
void output_result(FILE *outfile, double i_center, double j_center,
		   double i_offset, double j_offset, double conversion,
		   double coor_coeff);

/*int **mult, *square;*/

int main(int argc, char *argv[])
{
  FILE *mapfile, *parmfile, *infile1, *infile2, *outfile;
  ImageType image_type;
  int xsize, ysize;
  int box_height, box_width;
  int abs_range, skip_range, norm_range, sub_range;
  int skip_val, sub_val;
  int mask_min, mask_max;
  int overlap;
  double delta_t, resolution;
  int buffer_output;

  short **image1, **image2;

  int nsearches_x, nsearches_y; 
  double **i_offset, **j_offset, **coor_coeff;
  double i_off, j_off;
  double coeff;

  char *map_filename, *points_filename;
  unsigned char **map;

  double i_center, j_center;
  double conversion;

  int i_0, j_0;

  int i, j;

  /* Usage */

  if (argc != NARGS)
    usage();

  /* Open files */

  parmfile = open_read(argv[1]);
  infile1 = open_read(argv[2]);
  infile2 = open_read(argv[3]);
  outfile = open_write(argv[4]);

  /* Read values from parmfile */

  getparms(parmfile, argv[1], &xsize, &ysize, 
	   &box_width, &box_height,
	   &skip_val, &sub_val, 
	   &abs_range, &skip_range, &sub_range, &norm_range,
	   &overlap, &mask_min, &mask_max, &delta_t, &resolution, 
	   &image_type, &buffer_output, &map_filename, &points_filename);

  if (map_filename)
    mapfile = open_read(map_filename);

  conversion = (1000*resolution)/(36*delta_t);

  /* Allocate space for images */

  if (map_filename)
    map = getmem_2d(ysize, xsize*sizeof(unsigned char));
  image1 = getmem_2d(ysize, xsize*sizeof(short));
  image2 = getmem_2d(ysize, xsize*sizeof(short));

  /* Read in images */

  if (map_filename)
    readmap(map, xsize, ysize, mapfile, map_filename);
  switch(image_type) {
    case int1b:
      readimage_uint1b(image1, xsize, ysize, infile1, argv[2]);
      readimage_uint1b(image2, xsize, ysize, infile2, argv[3]);
      break;
    case int2bh:
      readimage_int2bh(image1, xsize, ysize, infile1, argv[2]);
      readimage_int2bh(image2, xsize, ysize, infile2, argv[3]);
      break;
    case int2bl:
      readimage_int2bl(image1, xsize, ysize, infile1, argv[2]);
      readimage_int2bl(image2, xsize, ysize, infile2, argv[3]);
  }
     
  /* Mask out pixels specified in the map */

  if (map_filename)
    for (j=0;j<ysize;j++)
      for (i=0;i<xsize;i++) 
	if (map[j][i]) {
	  image1[j][i] = 0;
	  image2[j][i] = 0;
        }

  /* Mask out values less than min and greater than max */

  for (j=0;j<ysize;j++)
    for (i=0;i<xsize;i++) {
      if (image1[j][i] < mask_min || image1[j][i] > mask_max)
        image1[j][i] = 0;
      if (image2[j][i] < mask_min || image2[j][i] > mask_max)
        image2[j][i] = 0;
    }

  /* Initialize lookup tables */ 

  /*mult = getmem_2d(TABLESIZE, TABLESIZE*sizeof(int));
  for (j=0;j<TABLESIZE;j++)
    for (i=0;i<TABLESIZE;i++)
      mult[j][i] = i*j;*/

 /* square = (int *)getmem(TABLESIZE*sizeof(int));
  for (i=0;i<TABLESIZE;i++)
    square[i] = i*i;*/

  /* Determine number of searches in each direction */ 

  nsearches_x = ((xsize - 2*abs_range - box_width)/
		(box_width - overlap)) + 1;
  nsearches_y = ((ysize - 2*abs_range - box_height)/
		(box_height - overlap)) + 1;

  /* Output file header */

  fprintf(outfile, "%s\t%s\n", argv[2], argv[3]);
  fprintf(outfile, "%d\t%d\t%d\t%d\t%f\n", nsearches_x, nsearches_y,
	  xsize, ysize, 1/conversion);

  /* Allocate space for offsets, and correlation coefficients */

  if (buffer_output) {
    i_offset = getmem_2d(nsearches_y, nsearches_x*sizeof(double));
    j_offset = getmem_2d(nsearches_y, nsearches_x*sizeof(double));
    coor_coeff = getmem_2d(nsearches_y, nsearches_x*sizeof(double));
  }

  /* Run Correlations */

  for (j=0;j<nsearches_y;j++)
    for (i=0;i<nsearches_x;i++) {
      i_0 = (box_width - overlap)*i + abs_range;
      j_0 = (box_height - overlap)*j + abs_range;
      find_target(image1, image2, xsize, ysize, i_0, j_0, 
		  box_height, box_width,
		  skip_val, sub_val, abs_range, skip_range, 
		  norm_range, sub_range, &i_off, &j_off, &coeff);
      if (buffer_output) {
        i_offset[j][i] = i_off;
        j_offset[j][i] = j_off;
        coor_coeff[j][i] = coeff;
      }
      else {
	i_center = i_0 + box_width/2. - 0.5;
	j_center = j_0 + box_height/2. - 0.5;
	output_result(outfile, i_center, j_center, i_off, j_off, 
		     conversion, coeff);
        fflush(outfile);
      }
    }

  if (buffer_output) 
    for (j=0;j<nsearches_y;j++)
      for (i=0;i<nsearches_x;i++) {
	i_center = (box_width - overlap)*i + abs_range + box_width/2. - 0.5;
	j_center = (box_height - overlap)*j + abs_range + box_height/2. - 0.5; 
	output_result(outfile, i_center, j_center, 
		     i_offset[j][i], j_offset[j][i], 
		     conversion, coor_coeff[j][i]);
      }

  return 0;
}

void usage()
{
  fprintf(stderr, "usage: motion <parmfile> <infile1> <infile2> <outfile>\n");
  exit(1);
}

void find_target(short **image1, short **image2, int xsize, int ysize,
		 int i_0, int j_0, int box_width, int box_height, 
		 int skip_val, int sub_val,
                 int abs_range, int skip_range, int norm_range, int sub_range,
		 double *i_offset, double *j_offset, double *coor_coeff)
{
  double i_off, j_off, coeff;

  i_off = 0.;
  j_off = 0.;
  coeff = -2.0;

  if (good_target(image1, i_0, j_0, box_width, box_height)) {
    if (skip_range) 
      ft_skip(image1, image2, i_0, j_0, box_width, box_height, 
	      skip_val, skip_range, &i_off, &j_off, &coeff);

    if (norm_range) 
      ft_norm(image1, image2, i_0, j_0, round(i_0+i_off), round(j_0+j_off),
	      box_width, box_height, norm_range, abs_range, 
	      &i_off, &j_off, &coeff);

    if (sub_range) 
      ft_sub(image1, image2, xsize, ysize, i_0, j_0, 
	     round(i_0+i_off), round(j_0+j_off), box_width, box_height, 
	     sub_val, sub_range, abs_range, 
	     &i_off, &j_off, &coeff);
  }

  *i_offset = i_off; 
  *j_offset = j_off; 
  *coor_coeff = coeff;
}

void ft_skip(short **image1, short **image2, int i_0, int j_0,
	     int box_width, int box_height, int skip_val, int skip_range,
	     double *i_offset, double *j_offset, double *coor_coeff)
{
  int best_i, best_j;
  double coeff, best_coeff;
  int i, j;

/*** May want to change best_i & j initial for when coeff = -2.0 ***/
  best_i = 0;
  best_j = 0;
  best_coeff = -2.0;

  for (j = j_0 - skip_range*skip_val;
       j < j_0 + skip_range*skip_val + skip_val;
       j += skip_val)
    for (i = i_0 - skip_range*skip_val;
	 i < i_0 + skip_range*skip_val + skip_val;
	 i += skip_val) {
      coeff = correlation(image1, image2, i_0, j_0, i, j, 
			  box_width, box_height);
      if (coeff > best_coeff) {
	best_coeff = coeff;
	best_i = i;
	best_j = j;
      }
    }

  *i_offset = best_i - i_0;
  *j_offset = best_j - j_0;
  *coor_coeff = best_coeff;
}

void ft_norm(short **image1, short **image2, 
	     int i_0, int j_0, int i_start, int j_start,
	     int box_width, int box_height, int norm_range, int abs_range,
	     double *i_offset, double *j_offset, double *coor_coeff)
{
  int best_i, best_j;
  double coeff, best_coeff;
  int i, j;

  best_i = 0;
  best_j = 0;
  best_coeff = -2.0;

  for (j=j_start-norm_range;j<j_start+norm_range+1;j++)
    for (i=i_start-norm_range;i<i_start+norm_range+1;i++) {
      if (abs(i-i_0) <= abs_range && abs(j-j_0) <= abs_range) {
        coeff = correlation(image1, image2, i_0, j_0, i, j, 
                            box_width, box_height);
        if (coeff > best_coeff) {
  	  best_coeff = coeff;
	  best_i = i;
	  best_j = j;
        }
      }
    }
    

  *i_offset = best_i - i_0;
  *j_offset = best_j - j_0;
  *coor_coeff = best_coeff;
}

void ft_sub(short **image1, short **image2, int xsize, int ysize, 
	    int i_0, int j_0, int i_start, int j_start, 
	    int box_width, int box_height,
	    int sub_val, int sub_range, int abs_range,
	    double *i_offset, double *j_offset, double *coor_coeff)
{
  short **image1_sub, **image2_sub;
  int best_i, best_j;
  int sw_i0_norm, sw_j0_norm, sw_in_norm, sw_jn_norm;
  int i_min_sub, j_min_sub, i_n_sub, j_n_sub;
  double coeff, best_coeff;
  int i, j;

  best_i = 0;
  best_j = 0;
  best_coeff = -2.0;

  sw_i0_norm = i_start - (sub_range + sub_val - 1)/sub_val;
  sw_j0_norm = j_start - (sub_range + sub_val - 1)/sub_val;
  sw_in_norm = 2*(sub_range + sub_val - 1)/sub_val + box_width;
  sw_jn_norm = 2*(sub_range + sub_val - 1)/sub_val + box_height;

  i_min_sub = (sub_val - sub_range%sub_val)%sub_val;
  j_min_sub = (sub_val - sub_range%sub_val)%sub_val;
  i_n_sub = 2*sub_range + 1;
  j_n_sub = 2*sub_range + 1;

  image1_sub = expand(image1, xsize, ysize, i_0, j_0, 
		      box_width, box_height, sub_val);
  image2_sub = expand(image2, xsize, ysize, sw_i0_norm, sw_j0_norm, 
		      sw_in_norm, sw_jn_norm, sub_val);
/*write_image(image1_sub, 37, 37, "outfile1");
*/

  for (j=j_min_sub;j<j_min_sub+j_n_sub;j++)
    for (i=i_min_sub;i<i_min_sub+i_n_sub;i++)
      if (abs(i/sub_val + sw_i0_norm - i_0) <= abs_range && 
	  abs(j/sub_val + sw_j0_norm - j_0) <= abs_range) {
        coeff = correlation(image1_sub, image2_sub, 0, 0, 
			    i, j, sub_val*(box_width-1)+1, 
			    sub_val*(box_height-1)+1);
/* fputc((char)(coeff*100 + 100), stderr);
*/
        if (coeff > best_coeff) {
	  best_coeff = coeff;
	  best_i = i;
	  best_j = j;
        }
      }

  freemem_2d(image1_sub, sub_val*(box_height-1)+1);
  freemem_2d(image2_sub, sub_val*(sw_jn_norm-1)+1);

  *i_offset = (double)best_i/sub_val + sw_i0_norm - i_0;
  *j_offset = (double)best_j/sub_val + sw_j0_norm - j_0;
  *coor_coeff = best_coeff;
}

int good_target(short **image, int i0, int j0, int box_width, int box_height)
{
  int good_pixels;
  int i, j;

  good_pixels = 0;

  for (j=j0;j<j0+box_height;j++)
    for (i=i0;i<i0+box_height;i++)
      if (image[j][i])
	good_pixels++;

  if (good_pixels < box_width*box_height*COR_THRESHOLD)
    return 0;
  else
    return 1;
}

short **expand(short **image, int xsize, int ysize, 
	       int i0, int j0, int in, int jn, int factor)
{
  double **tmp_image;
  short **new_image;
  int i, j;

  if (factor < 1)
    return NULL;

  tmp_image = getmem_2d(factor*(jn-1)+1, (factor*(in-1)+1)*sizeof(double)); 
/*for (j=0;j<factor*(jn-1)+1;j++)
for (i=0;i<factor*(in-1)+1;i++)
tmp_image[j][i] = 0;
*/

  /* Fill temporary array with starting pixels */
  for (j=0;j<jn;j++)
    for (i=0;i<in;i++)
      tmp_image[factor*j][factor*i] = pixel(image, xsize, ysize, i+i0, j+j0);

  /* Smooth pixels horizontally */
  for (j=0;j<factor*(jn-1)+1;j+=factor)
    for (i=0;i<factor*(in-1)+1;i++)
      if (i%factor) 
        tmp_image[j][i] = ((factor - i%factor)*tmp_image[j][i-i%factor]
 		          + (i%factor)*tmp_image[j][i-i%factor+factor])
			  /factor;

  /* Smooth pixels vertically */
  for (j=0;j<factor*(jn-1)+1;j++)
    for (i=0;i<factor*(in-1)+1;i++)
      if (j%factor)
        tmp_image[j][i] = ((factor - j%factor)*tmp_image[j-j%factor][i]
			  + (j%factor)*tmp_image[j-j%factor+factor][i])
			  /factor;

  /* Make final array */ 
  new_image = getmem_2d(factor*(jn-1)+1, (factor*(in-1)+1)*sizeof(short));
  for (j=0;j<factor*(jn-1)+1;j++)
    for (i=0;i<factor*(in-1)+1;i++)
      new_image[j][i] = round(tmp_image[j][i]);

  freemem_2d(tmp_image, factor*(jn-1)+1);

  return new_image;
}

void freemem_2d(void *array, int jn)
{
  int j;

  for (j=0;j<jn;j++)
    free(((void **)(array))[j]);
  free(array);
}

double correlation(short **image1, short **image2, int i1_0, int j1_0,
		   int i2_0, int j2_0, int box_width, int box_height)
{
  double sigma_i1, sigma_i2;
  double sigma_i12, sigma_i22;
  double sigma_i1i2;
  double numerator, denomenator;
  int n;
  int i, j;

  sigma_i1 = 0.;
  sigma_i2 = 0.;
  sigma_i12 = 0.;
  sigma_i22 = 0.;
  sigma_i1i2 = 0.;
  n = 0;

  for (j=0;j<box_height;j++)
    for (i=0;i<box_width;i++) 
      if (image1[j1_0+j][i1_0+i] && image2[j2_0+j][i2_0+i]) {
	n++;
	sigma_i1 += image1[j1_0+j][i1_0+i];
	sigma_i2 += image2[j2_0+j][i2_0+i];
	sigma_i12 += image1[j1_0+j][i1_0+i]*image1[j1_0+j][i1_0+i];
	sigma_i22 += image2[j2_0+j][i2_0+i]*image2[j2_0+j][i2_0+i];
	sigma_i1i2 += image1[j1_0+j][i1_0+i]*image2[j2_0+j][i2_0+i];
      }
  
  if (n < box_width*box_height*COR_THRESHOLD)
    return -2.0;

  denomenator = sqrt((n*sigma_i12-sigma_i1*sigma_i1)*
	             (n*sigma_i22-sigma_i2*sigma_i2));

  if (denomenator == 0.)
    return 2.0;

  numerator = n*sigma_i1i2-sigma_i1*sigma_i2;

  return numerator/denomenator;
}

void getparms(FILE *parmfile, char *parmfilename, int *xsize, int *ysize, 
	      int *box_width, int *box_height, int *skip_val, int *sub_val, 
	      int* abs_range, int *skip_range, int *sub_range, int *norm_range,
	      int *overlap, int *mask_min, int *mask_max, 
	      double *delta_t, double *resolution, 
	      ImageType *image_type, int *buffer_output, 
	      char **map_filename, char **points_filename)
{
  char buffer[BUFFER_LEN];
  int rng, skip_v, ps_rng, sub_v;
  double sub_rng;
  char it_string[STRING_LEN], buffer_string[STRING_LEN]; 
  char map_string[STRING_LEN], points_string[STRING_LEN];
  int errors;

  sub_rng = 0.;		/* Important.  Gets rid of core dump (FE). */
  errors = 0;

  fgets(buffer, BUFFER_LEN, parmfile);
  if (sscanf(buffer, "%d,%d\n", xsize, ysize) != 2) {
    fprintf(stderr, "Error in parameter file: %s, line 1:\n", parmfilename);
    fprintf(stderr, "\t%s", buffer);
    fprintf(stderr, "Correct Syntax: \n\t<xsize>, <ysize> [comments]\n");
    errors = 1;
  }

  fgets(buffer, BUFFER_LEN, parmfile);
  if (sscanf(buffer, "%d,%d\n", box_width, box_height) != 2) {
    fprintf(stderr, "Error in parameter file: %s, line 2:\n", parmfilename);
    fprintf(stderr, "\t%s", buffer);
    fprintf(stderr, "Correct Syntax: \n\t<box width>, <box height> ");
    fprintf(stderr, "[comments]\n");
    errors = 1;
  }

  fgets(buffer, BUFFER_LEN, parmfile);
  if (sscanf(buffer, "%d\n", &rng) != 1) {
    fprintf(stderr, "Error in parameter file: %s, line 3:\n", parmfilename);
    fprintf(stderr, "\t%s", buffer);
    fprintf(stderr, "Correct Syntax: \n\t<range> [comments]\n");
    errors = 1;
  }

  fgets(buffer, BUFFER_LEN, parmfile);
  if (sscanf(buffer, "%d,%d\n", &skip_v, &ps_rng) != 2) {
    fprintf(stderr, "Error in parameter file: %s, line 4:\n", parmfilename);
    fprintf(stderr, "\t%s", buffer);
    fprintf(stderr, "Correct Syntax: \n\t<skip value>, <post skip range> ");
    fprintf(stderr, "[comments]\n");
    errors = 1;
  }

  fgets(buffer, BUFFER_LEN, parmfile);
  if (sscanf(buffer, "%d,%lf\n", &sub_v, &sub_rng) != 2) {
    fprintf(stderr, "Error in parameter file: %s, line 5:\n", parmfilename);
    fprintf(stderr, "\t%s", buffer);
    fprintf(stderr, "Correct Syntax: \n\t<sub value>, <sub range> ");
    fprintf(stderr, "[comments]\n");
    errors = 1;
  }


  *skip_val = skip_v;
  *abs_range = rng;
  if (skip_v) {
      *skip_range = rng/skip_v;
    *norm_range = ps_rng;
  }
  else {
    *skip_range = 0;
    *norm_range = rng;
  }
  
  *sub_val = sub_v;
  *sub_range = round(sub_v*sub_rng); 

  fgets(buffer, BUFFER_LEN, parmfile);
  if (sscanf(buffer, "%d\n", overlap) != 1) {
    fprintf(stderr, "Error in parameter file: %s, line 6:\n", parmfilename);
    fprintf(stderr, "\t%s", buffer);
    fprintf(stderr, "Correct Syntax: \n\t<overlap> [comments]\n");
    errors = 1;
  }

  fgets(buffer, BUFFER_LEN, parmfile);
  if (sscanf(buffer, "%d,%d\n", mask_min, mask_max) != 2) {
    fprintf(stderr, "Error in parameter file: %s, line 7:\n", parmfilename);
    fprintf(stderr, "\t%s", buffer);
    fprintf(stderr, "Correct Syntax: \n\t<min>, <max> [comments]\n");
    errors = 1;
  }

  fgets(buffer, BUFFER_LEN, parmfile);
  if (sscanf(buffer, "%lf\n", delta_t) != 1) {
    fprintf(stderr, "Error in parameter file: %s, line 8:\n", parmfilename);
    fprintf(stderr, "\t%s", buffer);
    fprintf(stderr, "Correct Syntax: \n\t<time difference> [comments]\n");
    errors = 1;
  }

  fgets(buffer, BUFFER_LEN, parmfile);
  if (sscanf(buffer, "%lf\n", resolution) != 1) {
    fprintf(stderr, "Error in parameter file: %s, line 9:\n", parmfilename);
    fprintf(stderr, "\t%s", buffer);
    fprintf(stderr, "Correct Syntax: \n\t<resolution> [comments]\n");
    errors = 1;
  }

  fgets(buffer, BUFFER_LEN, parmfile);
  sscanf(buffer, "%s\n", it_string);

  if (!strcmp(it_string, "1"))
    *image_type = int1b;
  else if (!strcmp(it_string, "2") || !strcmp(it_string, "2h"))
    *image_type = int2bh;
  else if (!strcmp(it_string, "2l"))
    *image_type = int2bl;

  fgets(buffer, BUFFER_LEN, parmfile);
  sscanf(buffer, "%s\n", buffer_string);

  if (!strcasecmp(buffer_string, "yes") || !strcasecmp(buffer_string, "y"))
    *buffer_output = 1;
  else if (!strcasecmp(buffer_string, "no") || 
	   !strcasecmp(buffer_string, "n"))
    *buffer_output = 0;

  fgets(buffer, BUFFER_LEN, parmfile);
  sscanf(buffer, "%s\n", map_string);

  if (!strcmp(map_string, "none"))
    *map_filename = NULL;
  else
    *map_filename = map_string;
  
  fgets(buffer, BUFFER_LEN, parmfile);
  sscanf(buffer, "%s\n", points_string);

  if (!strcmp(points_string, "none"))
    *points_filename = NULL;
  else
    *points_filename = points_string;

  if (errors)
    exit(1);
}

FILE *open_read(char *filename)
{
  FILE *file;

  if (!(file = fopen(filename, "rb"))) {
    perror(filename);
    exit(1);
  }

  return file;
}

FILE *open_write(char *filename)
{
  FILE *file;

  if (!(file = fopen(filename, "wb"))) {
    perror(filename);
    exit(1);
  }

  return file;
}

void *getmem(int nbytes)
{
  void *mem;

  if (!(mem = (void *)malloc(nbytes))) {
    fprintf(stderr, "Out of memory.\n");
    exit(1);
  }

  return mem;
}

void *getmem_2d(int nrows, int nbytes) 
{
  void *mem;
  int j;

  if (!(mem = (void *)malloc(nrows*sizeof(void *)))) {
    fprintf(stderr, "Out of memory.\n");
    exit(1);
  }

  for (j=0;j<nrows;j++)
    ((void **)mem)[j] = getmem(nbytes);

  return mem;
}

void readimage_uint1b(short **image, int xsize, int ysize, FILE *file, 
		      char *filename)
{
  unsigned char *buffer;
  int i, j;

  buffer = getmem(xsize);
  for (j=0;j<ysize;j++) {
    if (fread(buffer, 1, xsize, file) != xsize) {
      fprintf(stderr, "%s: Unexexpected EOF.\n", filename);
      exit(1);
    }
    for (i=0;i<xsize;i++)
      image[j][i] = buffer[i];
  }
  free(buffer);
}

void readimage_int2bh(short **image, int xsize, int ysize, FILE *file,
		      char *filename)
{
  char *buffer;
  int i, j;

  buffer = getmem(2*xsize);
  for (j=0;j<ysize;j++) {
    if (fread(buffer, 1, 2*xsize, file) != 2*xsize) {
      fprintf(stderr, "%s: Unexexpected EOF.\n", filename);
      exit(1);
    }
    for (i=0;i<xsize;i++)
      image[j][i] = buffer[2*i]<<8 | (unsigned char)buffer[2*i+1];
  }
  free(buffer);
}

void readimage_int2bl(short **image, int xsize, int ysize, FILE *file,
		      char *filename)
{
  char *buffer;
  int i, j;

  buffer = getmem(2*xsize);
  for (j=0;j<ysize;j++) {
    if (fread(buffer, 1, 2*xsize, file) != 2*xsize) {
      fprintf(stderr, "%s: Unexexpected EOF.\n", filename);
      exit(1);
    }
    for (i=0;i<xsize;i++)
      image[j][i] = buffer[2*i+1]<<8 | (unsigned char)buffer[2*i];
  }
  free(buffer);
}

void write_image(short **image, int xsize, int ysize, char *filename)
{
  FILE *outfile;
  int j;

  outfile = open_write(filename);

  for (j=0;j<ysize;j++)
    fwrite(image[j], sizeof(short), xsize, outfile); 

  fclose(outfile);
}

void readmap(unsigned char **map, int xsize, int ysize, FILE *file,
	     char *filename)
{
  unsigned char *buffer;
  int i, j;

  buffer = getmem(xsize*sizeof(unsigned char));
  for (j=0;j<ysize;j++) {
    if (fread(buffer, sizeof(unsigned char), xsize, file) != xsize) {
      fprintf(stderr, "%s: Unexexpected EOF.\n", filename);
      exit(1);
    }
    for (i=0;i<xsize;i++)
      map[j][i] = buffer[i];
  }
  free(buffer);
}

void output_result(FILE *outfile, double i_center, double j_center,
		   double i_offset, double j_offset, double conversion,
		   double coor_coeff)
{

  if (coor_coeff == -2.0)
    fprintf(outfile, "%6.1f\t%6.1f\t%6.2f\t%6.2f\t%6.2f\n",
	    i_center, j_center, 0., 0., 0.);
  else
    fprintf(outfile, "%6.1f\t%6.1f\t%6.2f\t%6.2f\t%6.2f\n", i_center, j_center,
	    i_offset*conversion, -j_offset*conversion, coor_coeff);
}
