//Input-------------------------------------------------------------------------------------------------
#define WindowDimension 3						// this is the dimension of the window. 
#define PatchSigma 0.01                         // this is h squared , mentioned in the report
#define Sigma 0.05   							// this is the sigma squared , mentioned in the report
#define FILENAME "images/rasp_noise.csv"		// path to the csv of the image you want to use as input
//End of input------------------------------------------------------------------------------------------





#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//fast exp , sacrifising some accuracy
float expf_fastCPU(float a) {
  union { float f; int x; } u;
  u.x = (int) (12102203 * a + 1064866805);
  return u.f;
}

//Calculate the mean of each pixel , using the moving-window technique
void meansCPU(float* pixels , int *windowSize, int* iconLine, float *means){
	int xIndex = 0;
	for(xIndex = 0; xIndex<*iconLine*(*iconLine); xIndex++){


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
			
			mean = mean + result[i]*expf_fastCPU(-arithmitis/paronomastis)*0.5;
			tmp = tmp+expf_fastCPU(-arithmitis/paronomastis)*0.5;
		
		}
		
		mean = mean/tmp;
		means[xIndex] = mean; 
	}
}
//denoise the image using the formulas mentioned in the report
void denoiseCPU(float* pixels, float* sigma, int* imageDimension,int* windowDimension,float* means, float* result){
	int xIndex =0;
	for(xIndex=0; xIndex<*imageDimension*(*imageDimension); xIndex++){
		int windowSize;
		int imageSize;
		imageSize = *imageDimension*(*imageDimension);
		windowSize = *windowDimension*(*windowDimension); 
		float mean1 = means[xIndex];
		int i=0;
		float sumW = 0 ;
		float sumP = 0;
		for(i=0; i<imageSize; i++){
			float mean2 = means[i];
		
			float tmp = (mean1 - mean2)*(mean1-mean2)*(-1);
			float weight = exp(tmp/(*sigma));

			sumW = sumW + weight;
			sumP = sumP + weight*pixels[i];

		}
	sumP = sumP/sumW;
	result[xIndex] = sumP;
	}
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

	float* means = (float*)malloc(dimension*dimension*sizeof(float));

	
	
	
	


	

	float* finalPixel = (float*)malloc(sizeof(float)*dimension*dimension);
	clock_t start, end;
    double cpu_time_used;
    start = clock();

    //the whole algorithm gets executed here
    meansCPU(image, &windowDimension, &dimension, means);

    denoiseCPU(image, &sigma, &dimension, &windowDimension, means, finalPixel);
	//-------------------------------------

	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("program took %f seconds to execute \n", cpu_time_used);

    floatToCSV("imageAfter.csv", dimension, finalPixel);
    
	








	

	return 0; 

}