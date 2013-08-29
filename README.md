PoissonNoiseReduction
=====================

Poisson Noise Reduction example

This sample program is presented at the PDPTA 2013.

  "Poisson Observed Image Restoration using a Latent Variational Approximation with Gaussian MRF"
  H.Shouno and M.Okada

The presentation slides can see at

  http://www.slideshare.net/HAL9801/slideshelf

In this sample, prior implementation is just different from the proceeding.
The priro shown in the slide is correct implementation.

The implementation language is R with several libraries, "Matrix", "png".
"MyPoisson2D.R" is my local library source.
When you start on the R, set working directory to the repositry directory,
and you can run with

> source( 'Poisson2D13.R' )

If you'd like to use with command-line of several shells, please type

$ R --vanilla < Poisson2D13.R

After running, you can obtain several result files such like "Pois2D13_02_20_01.RData"
The numbers in the filename mean "Pois2D13_[min]_[max]_[trial].RData". 
Thus the example means, the first trial of restoration whose original contrast is [02, 20].

If you'd like see the image, please load the data in the R interpreter.

> load( 'Pois2D13_02_20_01.RData' )

The image data is assigned as followings

y0: Original Image (with 1D vector format)
y1: Poisson noised Image (with 1D vector format)
est: Restored Image (with 1D vector format)

The original image size is 64x64 gray scale image.
Hence, you can see the image as

> image( matrix( y1, 64, 64 ) )

The image should be shown in rotated. If you'd like to correct, please try the rotation the matrix
as t(matrix(y1,64,64)[64:1,]).
