//Input-------------------------------------------------------------------------------------------------
#define WindowDimension 3						// this is the dimension of the window. 
#define PatchSigma 0.1                          // this is h squared , mentioned in the report
#define Sigma 0.05   							// this is the sigma squared , mentioned in the report
#define FILENAME "images/rasp_noise.csv"		// path to the csv of the image you want to use as input
//End of input------------------------------------------------------------------------------------------



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define sharedSize 24

void printFloatArray(float* arr, int x, int y){
	int i;
	int j;

	for(i=0; i<x; i++){
		for(j=0; j<y; j++){
			printf("%f ", arr[i*y + j]);
		}
		printf("\n");
	}
}
//fast exp , sacrifising some accuracy
__device__ float expf_fast(float a) {
  union { float f; int x; } u;
  u.x = (int) (12102203 * a + 1064866805);
  return u.f;
}



//Calculate the mean of each pixel , using the moving-window technique
__global__ void means( float* pixels , int *windowSize, int* iconLine, float *means)
{
	

	int xIndex = threadIdx.x+ blockIdx.x* blockDim.x;
	int* index = &xIndex;
	int i, j;
	int lim = (*windowSize)/2;
	int dimension = *windowSize;

	float result[WindowDimension*WindowDimension];
	
	for(i=-lim; i<=lim; i++){
		for(j=-lim; j<=lim; j++){
			int location = (*index) + i +j*(*iconLine);

			bool outLeftorRight = (location/(*iconLine) != ((*index)/(*iconLine) + j) ) || location<0;
			bool outUpOrDown = location/(*iconLine) <0 || location/(*iconLine) >= (*iconLine);
			bool ifResult = !(outLeftorRight || outUpOrDown);
			int resultIndex = (j+lim)*(*windowSize) + (i+lim);

			//if the window is completely inside the image
			if(ifResult){ 
				result[resultIndex] = pixels[location];
			}
			
			//if part of the window is outside of the image
			else{
				

				location = (*index) - i -j*(*iconLine);
				
				outLeftorRight =location/(*iconLine) != ((*index)/(*iconLine) - j) || location<0;
				outUpOrDown = location/(*iconLine) <0 || location/(*iconLine) >= (*iconLine);

				if(!outLeftorRight && !outUpOrDown){
					result[resultIndex] = pixels[location];
					continue;
				}

				location = (*index) -i +j*(*iconLine);
				outLeftorRight =location/(*iconLine) != ((*index)/(*iconLine) + j) || location<0;
				outUpOrDown = location/(*iconLine) <0 || location/(*iconLine) >= (*iconLine);

				if(!outLeftorRight && !outUpOrDown){
					result[resultIndex] = pixels[location];
					continue;
				}

				location = (*index) +i -j*(*iconLine);
				result[resultIndex] = pixels[location];
			}

		}	
	}
	float mean=0;
	
	float patchSigma = PatchSigma;
	float tmp = 0;
	for(i=0; i<*windowSize*(*windowSize); i++){
		int x = i%dimension - dimension/2;
		int y = i/dimension - dimension/2;

		float fx = (float)x;
		float fy = (float)y;
		float arithmitis = fx*fx + fy*fy;
		float paronomastis = 2*M_PI*patchSigma;
		
		mean = mean + result[i]*expf_fast(-arithmitis/paronomastis)*0.5;
		tmp +=expf_fast(-arithmitis/paronomastis)*0.5;
	

	}
	
	mean = mean/tmp;
	means[xIndex] = mean; 

}

//denoise the image using the formulas mentioned in the report
__global__ void denoise_shared(float* pixels, float* sigma, int* imageDimension,int* windowDimension,float* means, float* result){
	int xIndex = threadIdx.x+ blockIdx.x* blockDim.x;
	int thread= threadIdx.x;
	int windowSize;
	int imageSize;
	imageSize = *imageDimension*(*imageDimension);
	windowSize = *windowDimension*(*windowDimension); 
	
	__shared__ float shared_pixels[sharedSize];
	__shared__ float shared_means[sharedSize];
	float mean1 = means[xIndex];
	int i=0;
	float sumW = 0 ;
	float sumP = 0;

	for(i=0; i<imageSize; i= i + sharedSize){
		if(i%sharedSize==0){
			if(threadIdx.x<sharedSize)
				shared_pixels[threadIdx.x] = pixels[threadIdx.x+i];
				

		}
		__syncthreads();
		if(i%24==0){
			if(threadIdx.x<sharedSize)
				shared_means[threadIdx.x] = means[threadIdx.x+i];
				

		}
		__syncthreads();
		
		int j =0; 
		for(j=0; j<sharedSize; j++){
			float mean2 = shared_means[j];
			float tmp = (mean1 - mean2)*(mean1-mean2)*(-1);
			float weight = exp(tmp/(*sigma));

			sumW = sumW + weight;
			sumP = sumP + weight*shared_pixels[j];
		}
	}
	sumP = sumP/sumW;
	result[xIndex] = sumP;
	
}

//read the csv and put it in a float array
float* readCSVfile(char* filename, int* dimension){
	FILE *file;
	file = fopen(filename, "r");
	int local_dimension = 256; 
	float* result = (float*)malloc(local_dimension*local_dimension*sizeof(float));  // dont forget to free the memory afterwards
	if(file == NULL){
		printf("The file could not be opened/ does not exits");
		exit(1);
	}

	char* line =NULL;
	size_t length = 0 ;
	ssize_t read;
	const char delimeters[] = ", ";

	char* number ;
	int outer_counter = 0;
	while ((read = getline(&line, &length, file) != -1)){
		
		number = strtok(line, delimeters);
		float fnumber = atof(number);
		int counter = 0;

		while( number !=NULL){
			result[counter + outer_counter*local_dimension] = fnumber;
			counter++;
		

			number = strtok(NULL, delimeters);
			if(number != NULL)	
				fnumber = atof(number);
		}

		if(counter != local_dimension){
			local_dimension = counter;
			
			result =  (float*)realloc(result, counter*counter*sizeof(float));
			
		}


		outer_counter++;
		

	}

	fclose(file);
	if(line)
		free(line);

	*dimension = local_dimension;
	return result;
}



//save a float array to csv
void floatToCSV(char* filename, int dimension, float* arr){
	int i;
	FILE* fileWriter;
	fileWriter = fopen(filename, "a");
	if(fileWriter == NULL){
		printf("something went wrong when saving");
		exit(1);
	}
	for(i=0; i<dimension*dimension ; i++){
	
		if(i!=0 && i%dimension==0)
			fputs("\n", fileWriter);
		if(arr[i]!=arr[i])
			arr[i]=0;
		
		char stringFl[10];
		sprintf(stringFl, "%f", arr[i]);
		fputs( stringFl, fileWriter );
		fputs( ",", fileWriter );
	}
	fclose(fileWriter);
}
int main(void){

	
	float sigma = Sigma;
	int windowDimension = WindowDimension;
	int dimension ; 
	float *image ;
	image = readCSVfile(FILENAME, &dimension);

	
	
	//device pointers 
	float *d_image;
	float *d_sigma;
	int* d_windowDimension;
	int* d_imageDimension;
	float* d_weights;
	float* d_finalImage;
	float* d_means;

	cudaMalloc((void**)&d_means, dimension*dimension*sizeof(float));
	cudaMalloc((void**)&d_finalImage, dimension*dimension*sizeof(float));

	cudaMalloc((void**)&d_image, dimension*dimension*sizeof(float));
	cudaMemcpy(d_image, image, dimension*dimension*sizeof(float), cudaMemcpyHostToDevice);
	
	cudaMalloc((void**)&d_sigma, sizeof(float));
	cudaMemcpy(d_sigma	, &sigma, sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_windowDimension, sizeof(int));
	cudaMemcpy(d_windowDimension, &windowDimension, sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_imageDimension, sizeof(int));
	cudaMemcpy(d_imageDimension, &dimension, sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_weights, dimension*dimension*sizeof(float));


	

	float* finalPixel = (float*)malloc(sizeof(float)*dimension*dimension); //final result 
	clock_t start, end;
    double cpu_time_used;

    start = clock();
    //the whole algorithm gets executed here
    means<<<dimension, dimension>>>(d_image, d_windowDimension, d_imageDimension, d_means);
    denoise_shared<<<dimension , dimension>>>(d_image, d_sigma, d_imageDimension, d_windowDimension,d_means, d_finalImage);
    //-------------------------------------
    cudaDeviceSynchronize();
    cudaError_t error = cudaGetLastError();
  	if(error != cudaSuccess)
  	{
    // print the CUDA error message and exit
    	printf("CUDA error: %s\n", cudaGetErrorString(error));
   		exit(-1);
  	}
    cudaMemcpy(finalPixel, d_finalImage, sizeof(float)*dimension*dimension, cudaMemcpyDeviceToHost);
	

	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("program took %f seconds to execute \n", cpu_time_used);

    floatToCSV("imageAfter.csv", dimension, finalPixel);
    
	cudaFree(d_weights);
	cudaFree(d_image);
	cudaFree(d_sigma);
	cudaFree(d_windowDimension);
	cudaFree(d_imageDimension); 






	return 0; 

}