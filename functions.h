#ifndef FUNCTIONS_H
#define FUNCTIONS_H


// you will need to implement these yourself for the assignment
//double* readWavFile(int *arraySize, int *channels, char *filename);	
//void readWavFileHeader(int *channels, int *numSamples, FILE *inputFile);

void writeWavFile(double *outputArray, int outputArraySize, int channels, char *filename);
void writeWavFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile);

void readWavFileHeader(int *channels, int *numberSamples, FILE *inputFile);
double* readWavFile(int *arraySize, int *channels, char *inputFilename);

bool isPowerOfTwo(int n);
int  nextPowTwo ( int  x );
double* convolution(double* inputArray, int inputArraySize, double* IRArray, int IRArraySize, int* outputArraySize);

size_t fwriteIntLSB(int data, FILE *stream);
int freadIntLSB(FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);
short int freadShortLSB(FILE *stream);

#endif