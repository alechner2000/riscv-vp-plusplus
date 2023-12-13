#include "irq.h"
#include "filter.h"

void camera_irq_handler(void);
void StartCamera(void);
void StopCamera(void);
void data_in(unsigned char *image, unsigned int i);
void data_out(unsigned char *edge, unsigned int i);

static volatile char *const CAMERA_START_ADDR = (char *const)0x51000000;
static volatile char *const CAMERA_END_ADDR   = (char *const)0x51ffffff;
static volatile unsigned char *const CAMERA_FRAME_BUFFER_ADDR =
			(unsigned char *const)(CAMERA_START_ADDR + 0x000000);
static volatile uint32_t *const CAMERA_CAPTURE_INTERVAL_REG_ADDR =
			(uint32_t *const)(CAMERA_START_ADDR + 0xff0000);
static volatile uint32_t *const CAMERA_WIDTH_REG_ADDR  =
			(uint32_t *const)(CAMERA_START_ADDR + 0xff0004);
static volatile uint32_t *const CAMERA_HEIGHT_REG_ADDR =
			(uint32_t *const)(CAMERA_START_ADDR + 0xff0008);
const uint32_t CAMERA_IRQ_NUMBER = 5;
volatile unsigned int frames_captured = 0;

/* Gaussian kernel (computed at beginning, then constant) */
float Gaussian_kernel[WINSIZE] = {0.0};
int Gaussian_kernel_center = 0;

int main(void)
{
   unsigned char image[SIZE];
   unsigned char edge[SIZE];
   unsigned int i;
   int windowsize; /* Dimension of the gaussian kernel. */

   StartCamera(); /* turn on the camera device */

   /****************************************************************************
   * Create a 1-dimensional gaussian smoothing kernel.
   ****************************************************************************/
   if(VERBOSE) printf("Computing the gaussian smoothing kernel.\n");
   make_gaussian_kernel(SIGMA, Gaussian_kernel, &windowsize);
   Gaussian_kernel_center = windowsize / 2;

   for(i=0; i<IMG_NUM; i++)
   {
      /*************************************************************************
      * Input a frame.
      *************************************************************************/
      printf("Input image frame %u.\n", i+1);
      data_in(image, i);

      /*************************************************************************
      * Perform the edge detection. All of the work takes place here.
      *************************************************************************/
      if(VERBOSE) printf("Starting Canny edge detection.\n");
      canny(image, ROWS, COLS, SIGMA, TLOW, THIGH, edge);

      /*************************************************************************
      * Output a frame.
      *************************************************************************/
      printf("Output image frame %u.\n", i+1);
      data_out(edge, i);
   }

   StopCamera(); /* turn on the camera device */

   return(0); /* exit cleanly */
}

void camera_irq_handler(void)
{
    frames_captured++;
}

void StartCamera(void) /* configure and turn on the camera */
{
  // TODO
}

void StopCamera(void)
{
  // TODO
}


void data_in(unsigned char *image, unsigned int i)
{
   /*************************************************************************
   * Grab an image from the camera
   *************************************************************************/
  // TODO
}

void data_out(unsigned char *edge, unsigned int i)
{
   char outfilename[128]; /* Name of the output "edge" image */
   unsigned int n;

   /*************************************************************************
   * Write out an edge image to a file.
   *************************************************************************/
   n = i % AVAIL_IMG;
   sprintf(outfilename, IMG_OUT, n+1);
   if(VERBOSE) printf("Writing the edge image %s.\n", outfilename);
   if(write_pgm_image(outfilename, edge, ROWS, COLS, "", 255) == 0){
      fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
      exit(1);
   }
}
