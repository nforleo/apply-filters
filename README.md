# apply-filters
C++ program to apply filters (sequentially, parallel, and CUDA)

Note: need a GPU processor to run conv2d.cu (replace 'seq'

To compile program:

g++ -std=c++11 -Wall -pthread -conv2d.cc -o prog -ltiff

To run program (example):

./prog 'images/sunflower.tiff' 'final.tiff' 'filters/edge-detection.txt' seq

