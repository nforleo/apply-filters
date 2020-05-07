#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <math.h>
#include <sstream>
#include <cassert>
#include "tiffio.h"

typedef struct {
	int n_threads;
	uint32 *raster;
	uint32 w;
	uint32 h;
	const float *filter;
	int f_len;
	int thread_id;
} params;

void *par_filter_img(void *args) {
	params *param = (params *) args;
	
	// extract variables
	uint32 w = param->w;
	uint32 h = param->h;
	int f_len = param->f_len;
	int thread_id = param->thread_id;
	uint32 n_threads = param->n_threads;
	int start_ind = thread_id*(h / n_threads);
	const float *filter = param->filter;
	uint32 *raster = param->raster;
	
	int par_height = (h / n_threads) + ((h / n_threads)*thread_id);

	int f_width = sqrt(f_len);
	
	// define copy of image
	uint32 *copy = new uint32[w*h*sizeof(uint32)];
	for(int i = 0; i < w*h; i++) {
		copy[i] = raster[i];
	}
	
	bool isEnd = false;
	int last_ind = w*h - 1;
	int filter_wall = floor(f_width / 2);
	uint32 last_upd = last_ind-(w*filter_wall+filter_wall);
	
	// determine where to start processing on each thread
	int row_ind;
	if(thread_id == 0) {
		row_ind = start_ind + filter_wall;
	} else {
		row_ind = start_ind;
	}
	int col_ind = filter_wall;
	
	// Defaults alpha value to max
	uint32 alpha = 255;
	
	for(uint32 i = row_ind; i < par_height; i++) {
		for(uint32 j = col_ind; j < w; j++) {
			float r_sum = 0;
			float g_sum = 0;
			float b_sum = 0;
			for(int a = 0; a < f_width; a++) {
				for(int b = 0; b < f_width; b++) {
					int curr_ind = ((i+a-filter_wall)*w+(j+b-filter_wall));
					int filter_ind = a*f_width+b;
					if(!isEnd) {
						r_sum += TIFFGetR(copy[curr_ind]) * filter[filter_ind];
						g_sum += TIFFGetG(copy[curr_ind]) * filter[filter_ind];
						b_sum += TIFFGetB(copy[curr_ind]) * filter[filter_ind];
						if(curr_ind == last_ind) {
							isEnd = true;
						}
					}
				}	
			}
			uint32 ij = i*w+j;
		
			if(ij <= last_upd) {
				if(ij !=  (i+filter_wall)*w-filter_wall) {
					// Ugly and redundant, but prevents total refactoring after stupid mistake
					int roundedR = r_sum; //round(r_sum);
					int roundedG = g_sum; //round(g_sum);
					int roundedB = b_sum; //round(b_sum);
					// Make sure RGB values don't go above 255 or below 0
					if(roundedR > 255) {
						roundedR = 255;
					} else if (roundedR < 0) {
						roundedR = 0;
					}
					if(roundedG > 255) {
						roundedG = 255;
					} else if(roundedG < 0) {
						roundedG = 0;
					}
					
					if(roundedB > 255) {
						roundedB = 255;
					} else if(roundedB < 0) {
						roundedB = 0;
					}
					
					//alpha = TIFFGetA(raster[ij]);
					raster[ij] = alpha << 24 | ((uint32) roundedB) << 16 | ((uint32) roundedG) << 8 | ((uint32) roundedR);
					//raster[ij] = 0xFF000000 | ((uint32) roundedB) << 16 | ((uint32) roundedG) << 8 | ((uint32) roundedR);
				}
			}
		}
	}
	
}

// saves TIFF file from data in `raster`
void save_tiff(const char *fname, uint32 *raster, uint32 w, uint32 h) {
    TIFF *tif = TIFFOpen(fname, "w");
    if (! raster) {
        throw std::runtime_error("Could not open output file");
    }
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 4);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFWriteEncodedStrip(tif, 0, raster, w*h*4);
    TIFFClose(tif);
}

// loads image data from `fname` (allocating dynamic memory)
// *w and *h are updated with the image dimensions
// raster is a matrix flattened into an array using row-major order
// every uint32 in the array is 4 bytes, enconding 8-bit packed ABGR
// A: transparency attribute (can be ignored)
// B: blue pixel
// G: green pixel
// R: red pixel
uint32 *load_tiff(const char *fname, uint32 *w, uint32 *h) {
    TIFF *tif = TIFFOpen(fname, "r");
    if (! tif) {
        throw std::runtime_error("Could not open input file");
    }
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, h);
    uint32 *raster = (uint32 *) _TIFFmalloc(*w * *h * sizeof (uint32));
    if (! raster) {
        TIFFClose(tif);
        throw std::runtime_error("Memory allocation error");
    }
    if (! TIFFReadRGBAImageOriented(tif, *w, *h, raster, ORIENTATION_TOPLEFT, 0)) {
        TIFFClose(tif);
        throw std::runtime_error("Could not read raster from TIFF image");
    }
    TIFFClose(tif);
    return raster;
}

void filter_image_seq(uint32 *raster, uint32 w, uint32 h, const float *filter, int f_len) {
    // to get RGB values from a pixel, you can either use bitwise masks
    // or rely on the following macros:
    // TIFFGetR(raster[i]) red
    // TIFFGetG(raster[i]) green
    // TIFFGetB(raster[i]) blue 
    // TIFFGetA(raster[i]) this value should be ignored
    //
    // to modify RGB values from a pixel, you can use bitwise shifts or masks
    // each pixel stores values in the order ABGR
    // raster[i] = 0xFF0000FF;  // assigns Red to a pixel i
    // raster[i] = 0xFF00FF00;  // assigns Green to a pixel i
    // raster[i] = 0xFFFF0000;  // assigns Blue to a pixel i
    // note that the first byte is always FF (alpha = 255)
    //
    // TODO: here you will filter the image in raster
    //
	
	/**
	*	Find the width of the filter
	*	Used to cycle through filter array
	*/
	int f_width = sqrt(f_len);
	
	/**
	*	Create copy of image
	*	Used to calcualte new pixel value without changing pixel
	* 	values of image
	*/
	uint32 *copy = new uint32[w*h*sizeof(uint32)];

	for(int i = 0; i < w*h; i++) {
		copy[i] = raster[i];
	}
	

	
	/**
	*	Determine if we need to continue with the image processing
	*/
	bool isEnd = false;
	int last_ind = w*h - 1;
	int filter_wall = floor(f_width / 2);
	uint32 last_upd = last_ind-(w*filter_wall+filter_wall);
	
	// Default alpha value to max
	uint32 alpha = 255;
	
	for(uint32 i = filter_wall; i < h; i++) {
		for(uint32 j = filter_wall; j < w; j++) {
			float r_sum = 0;
			float g_sum = 0;
			float b_sum = 0;
			for(int a = 0; a < f_width; a++) {
				for(int b = 0; b < f_width; b++) {
					int curr_ind = ((i+a-filter_wall)*w+(j+b-filter_wall));
					int filter_ind = a*f_width+b;
					if(!isEnd) {
						r_sum += TIFFGetR(copy[curr_ind]) * filter[filter_ind];
						g_sum += TIFFGetG(copy[curr_ind]) * filter[filter_ind];
						b_sum += TIFFGetB(copy[curr_ind]) * filter[filter_ind];
						if(curr_ind == last_ind) {
							isEnd = true;
						}
					}
				}	
			}
			uint32 ij = i*w+j;
		
			if(ij <= last_upd) {
				if(ij !=  (i+filter_wall)*w-filter_wall) {
					// Ugly and redundant, but prevents total refactoring after stupid mistake
					int roundedR = r_sum; //round(r_sum);
					int roundedG = g_sum; //round(g_sum);
					int roundedB = b_sum; //round(b_sum);
					
					if(roundedR > 255) {
						roundedR = 255;
					} else if (roundedR < 0) {
						roundedR = 0;
					}
					if(roundedG > 255) {
						roundedG = 255;
					} else if(roundedG < 0) {
						roundedG = 0;
					}
					
					if(roundedB > 255) {
						roundedB = 255;
					} else if(roundedB < 0) {
						roundedB = 0;
					}
					
					raster[ij] = alpha << 24 | ((uint32) roundedB) << 16 | ((uint32) roundedG) << 8 | ((uint32) roundedR);
				}
			}
		}
	}
}

void filter_image_par(uint32 *raster, uint32 w, uint32 h, const float *filter, int f_len, int n_threads) {
    //
    // TODO: here you will filter the image in raster using threads
    //
	
	// params list
	/*
	typedef struct {
		int n_threads;
		uint32 *raster;
		uint32 w;
		uint32 h;
		const float *filter;
		int f_len;
		int thread_id;
	} params;
	*/
	
	
	pthread_t *threads = new pthread_t[n_threads];
	
	// Define param data and create threads
	for(int i = 0; i < n_threads; i++) {
		params *data = new params;
		data->n_threads = n_threads;
		data->raster = raster;
		data->w = w;
		data->h = h;
		data->filter = filter;
		data->f_len = f_len;
		data->thread_id = i;
		
		int ret = pthread_create(&threads[i], NULL, par_filter_img, (void *)data);
		
		assert(! ret);
	}
    
	// wait until each threads finishes 
    for (int i = 0 ; i < n_threads ; i ++) {
        pthread_join(threads[i], NULL);
    }
	
}

float *load_filter(const char *fname, int *n) {
    std::ifstream myfile(fname);
    if (! myfile) {
        throw std::runtime_error("Could not open filter file");
    }
    myfile >> *n;
    float *filter = new float[*n];
    for (int i = 0 ; i < *n ; i++) myfile >> filter[i];
    myfile.close();
    return filter;
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cout << "Usage:\t./filter <in_fname> <out_fname> <filter_fname> <algo> <n_threads>" << std::endl;
        std::cout << "<in_fname> path to the input image" << std::endl;
        std::cout << "<out_fname> path to the output image" << std::endl;
        std::cout << "<filter_fname> path to the filter file" << std::endl;
        std::cout << "<algo> whether to use the sequential (seq) or parallel algorithm (par)" << std::endl;
        std::cout << "<n_threads> number of threads to use" << std::endl;
        return 0;
    }

    uint32 width, height;
    int n_threads = std::stoi(argv[5]);

    // loads the filter
    int f_len;
    float *filter = load_filter(argv[3], &f_len);
    
    // loads image bytes from file name supplied as a command line argument
    // this function allocates memory dynamically
    uint32 *image = load_tiff(argv[1], &width, &height);
	
    // measure time of the algorithm
    auto start = std::chrono::high_resolution_clock::now();
    if (! std::strcmp(argv[4], "seq")) {
        // call the sequential implementation
        filter_image_seq(image, width, height, filter, f_len);
    } else if (! std::strcmp(argv[4], "par")) {
        // TODO: call the parallel implementation
		filter_image_par(image, width, height, filter, f_len, n_threads);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << diff.count() << std::endl;   

    // save new file with filtered image
    save_tiff(argv[2], image, width, height);

    // frees memory allocated by load_filter and load_tiff
    delete [] filter;
    _TIFFfree(image);

    return 0;
}