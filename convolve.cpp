/*
	Produces a single tone, and writes the result in .wav file format to the filename given as user input
	Note that the specified filename should have a .wav ending
	
	Usage: ./testtone outputFile
	
	Contains functions for writing .wav files that can be used for the assignment
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>	// includes sin
#include <string>
#include <fstream>
#include <iostream>
#include "functions.h"
#include <complex>
#include <valarray>

// CONSTANTS ******************************

#define PI              	3.14159265358979

// Frequency of tone to be created (A = 440 Hz)
#define TONE_FREQUENCY		440

// Duration of tone to be created (3 seconds for now)
#define SECONDS				3

//  Standard sample rate in Hz
#define SAMPLE_RATE     	44100.0

//  Standard sample size in bits
#define BITS_PER_SAMPLE		16

// Standard sample size in bytes		
#define BYTES_PER_SAMPLE	(BITS_PER_SAMPLE/8)

// Rescaling factor to convert between 16-bit shorts and doubles between -1 and 1
#define MAX_SHORT_VALUE		32768

// Number of channels
#define MONOPHONIC			1
#define STEREOPHONIC		2

// Offset of the fmt chunk in the WAV header
#define FMT_OFFSET			12

using namespace std;


typedef complex<double> Complex;
typedef valarray<Complex> CArray;


int main(int argc, char **argv) 
{
	char *outputFileName, *inputFileName, *IRFileName;

	if (argc != 4) 
    {
		printf("Please enter a filename\n");
		exit(-1);
	}

	inputFileName = argv[1];
	IRFileName = argv[2];
    outputFileName = argv[3];


    int inputChannels, inputArraySize;
    double *inputArray;
    
    printf("\nReading input from input file %s....\n", inputFileName);
    inputArray = readWavFile(&inputArraySize, &inputChannels, inputFileName);

    int IRChannels, IRArraySize;
    double *IRArray;

    printf("\nReading input from IR file %s....\n", IRFileName);
    IRArray = readWavFile(&IRArraySize, &IRChannels, IRFileName);
    

    int outputChannels = 1;

    if (IRChannels == 2)
    {
        inputChannels = 2;
        outputChannels = 2;
        inputArray = convertToStereo(inputArray, &inputArraySize, inputArraySize);
    }

    
   
    int outputArraySize = inputArraySize + IRArraySize - 1;

    double *outputArray = convolution(inputArray, inputArraySize, IRArray, IRArraySize, &outputArraySize);

	printf("Writing result to file %s...\n", outputFileName);
    cout << outputChannels << endl;
    writeWavFile(outputArray, outputArraySize, outputChannels, outputFileName);
    
    printf("Finished");

}

//The code for the following fft and ifft functions were taken from https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B. 

void fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;
 
    // divide
    CArray even = x[slice(0, N/2, 2)];
    CArray  odd = x[slice(1, N/2, 2)];
 
    // conquer
    fft(even);
    fft(odd);
 
    // combine
    double toMultiply = -2 * PI * 1/N;
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = polar(1.0, toMultiply*k) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

void ifft(CArray& x)
{
    // conjugate the complex numbers
    x = x.apply(conj);
 
    // forward fft
    fft( x );
 
    // conjugate the complex numbers again
    x = x.apply(conj);
 
    // scale the numbers
    x /= x.size();
}
\

double* convolution(double* inputArray, int inputArraySize, double* IRArray, int IRArraySize, int* outputArraySize)
{
    double *outputArray = new double[*outputArraySize];
    int complexArraySize = *outputArraySize;

    if (!isPowerOfTwo(complexArraySize))
    {
        complexArraySize = nextPowTwo(complexArraySize);
    }

    double* inputPadded = new double[complexArraySize];
    double* IRPadded = new double[complexArraySize];

   for (int i = 0; i < inputArraySize; i++)
    {
        inputPadded[i] = inputArray[i];
    }

    for (int i = inputArraySize; i < complexArraySize; i++ )
    {
        inputPadded[i] = 0;
    }

    for (int i = 0; i < IRArraySize; i++)
    {
        IRPadded[i] = IRArray[i];
    }

    for (int i = IRArraySize; i < complexArraySize; i++ )
    {
        IRPadded[i] = 0;
    }

    double* inputComplex = new double[complexArraySize*2];
    double* IRComplex= new double[complexArraySize*2];

    for (int i = 0; i < complexArraySize; i++) 
    {
        inputComplex[i*2] = inputPadded[i];
        inputComplex[i*2+1] = 0;
    }

    for (int i = 0; i < complexArraySize; i++) {
        IRComplex[i*2] = IRPadded[i];
        IRComplex[i*2+1] = 0;
    }


    CArray inputData((Complex*)inputComplex, complexArraySize);
    CArray IRData((Complex*)IRComplex, complexArraySize);

    fft(inputData);
    fft(IRData);

    double* outputData = new double[*outputArraySize];
    int i;
    for ( i = 0; i < inputData.size() - 1; i+=2) 
    {
        inputData[i] = inputData[i]*IRData[i];
        inputData[i+1] = inputData[i+1]*IRData[i+1];
    }
    if (i == (inputData.size() - 1))
    {
        inputData[i] = inputData[i]*IRData[i];
    }

    ifft(inputData);
    int j;
    for (j = 0; j < *outputArraySize - 1; j+=2) 
    {

        outputData[j] = inputData[j].real();
        outputData[j+1] = inputData[j+1].real();

    }

    if (j == ( *outputArraySize - 1))
    {
            outputData[j] = inputData[j].real();
    }

    return outputData;
}

int  nextPowTwo ( int  x )
{
    int  c =  0 ;
    while ( x > 0 )
    {
        c++ ;
        x  =  x  >>  1 ;
    }
    return  (1 << c) ;
}

//Taken from https://www.geeksforgeeks.org/program-to-find-whether-a-given-number-is-power-of-2/
// Function to check if x is power of 2
bool isPowerOfTwo(int n)
{
   if(n==0)
   return false;
 
   return (ceil(log2(n)) == floor(log2(n)));
}

double* convertToStereo(double *inputArray, int *outputArraySize, int inputArraySize) 
{
    *outputArraySize = 2 * inputArraySize;
    double *outputArray = new double [*outputArraySize];
    for (int i = 0; i < *outputArraySize-1; i+=2)
    {
        outputArray[i] = inputArray[i/2];
        outputArray[i+1] = inputArray[i/2];
    }

    return outputArray;
}


//The code for readWavFileHeader along with readWavFile was taken straight from the Bardia's tutorial

void readWavFileHeader(int *channels, int *numSamples, FILE *inputFile)
{
    int sampleRate;
    int bytesPerSecond;
    int dataChunkSize;
    unsigned char buffer[64];

    fread(buffer, sizeof(unsigned char), FMT_OFFSET, inputFile);

    freadIntLSB(inputFile);
    int fmtSize = freadIntLSB(inputFile);
    freadShortLSB(inputFile);

    *channels = freadShortLSB(inputFile);
    sampleRate = freadIntLSB(inputFile);
    bytesPerSecond = freadIntLSB(inputFile);

    int frameSize = freadShortLSB(inputFile);
    int bitRate = freadShortLSB(inputFile);

    if (bitRate != BITS_PER_SAMPLE)
    {
        printf("Error: bit rate of provided WAV file is not 16. Exiting.");
        exit(-1);
    }
    if (sampleRate != SAMPLE_RATE)
    {
        printf("Error: sample rate of provided WAV file is not 44.1 KHz. Exiting");
        exit(-1);
    }
    fread(buffer, sizeof(unsigned char), fmtSize - 12, inputFile);
    dataChunkSize = freadIntLSB(inputFile);
    printf("Data chunk size: %d\n", dataChunkSize);
    *numSamples = dataChunkSize / (BYTES_PER_SAMPLE * (*channels));
}

double *readWavFile(int *arraySize, int *channels, char *filename)
{
    double *array;
    FILE *inputFileStream = fopen(filename, "rb");
    if (inputFileStream == NULL)
    {
        printf("File %s could not be opened for reading\n", filename);
        exit(-1);
    }
    int numSamples;
    readWavFileHeader(channels, &numSamples, inputFileStream);
    printf("Channels: %d\n", *channels);
    printf("Number of samples: %d\n", numSamples);
    if (numSamples <= 0)
    {
        printf("The file %s doesn't contain any samples. Exiting the program.\n", filename);
        exit(0);
    }
    *arraySize = numSamples * (*channels);
    array = new double[*arraySize];
    short *intArray = new short[*arraySize];
    int count = fread(intArray, BYTES_PER_SAMPLE, numSamples, inputFileStream);
   
    int largest = 0;
    for (int i = 0; i < *arraySize; i++)
    {
        if (intArray[i] > largest)
        {
            largest = intArray[i];
        }
    }
    for (int i = 0; i < *arraySize; i++)
    {
        array[i] = ((double)intArray[i]) / largest;
    }
   
    delete[] intArray;
    return array;
} 


/*
Writes the header for a WAV file with the given attributes to 
 the provided filestream
*/

void writeWavFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile) 
{
    // Note: channels is not currently used. You will need to add this functionality
	// yourself for the bonus part of the assignment
	
	/*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = numberSamples * BYTES_PER_SAMPLE;
	
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
	
    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * BYTES_PER_SAMPLE;
	
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
      
    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */

    if (channels == 1)
    {
        fwriteShortLSB((short) MONOPHONIC, outputFile);
    }
    else
    {
        fwriteShortLSB((short) STEREOPHONIC, outputFile); 
    }
    

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}


/*
Creates a WAV file with the contents of the provided outputArray as the samples, and writes
it to the given filename
 */

void writeWavFile(double *outputArray, int outputArraySize, int channels, char *filename) {
    // Note: channels is not currently used. You will need to add this functionality
	// yourself for the bonus part of the assignment

  //open a binary output file stream for writing
    FILE *outputFileStream = fopen(filename, "wb");
    if (outputFileStream == NULL) {
      printf("File %s cannot be opened for writing\n", filename);
        return;
    }

    //create an 16-bit integer array to hold rescaled samples
    short *intArray = new short[outputArraySize];

    //find the largest entry and uses that to rescale all other
    // doubles to be in the range (-1, 1) to prevent 16-bit integer overflow
    double largestDouble = 1;
    for (int i=0; i< outputArraySize; i++) {
		if (abs(outputArray[i]) > largestDouble) {
			largestDouble = abs(outputArray[i]);
		}
    }

    for (int i=0; i<outputArraySize; i++) {
		intArray[i] = (short) ((outputArray[i]/largestDouble)*MAX_SHORT_VALUE);
    }
	
    int numSamples = outputArraySize;

	// actual file writing
    writeWavFileHeader(channels, numSamples, SAMPLE_RATE, outputFileStream);
    fwrite(intArray, sizeof(short), outputArraySize, outputFileStream);
    
    //clear memory from heap
    delete[] intArray;
}


//writes an integer to the provided stream in little-endian form
size_t fwriteIntLSB(int data, FILE *stream) {
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}


//reads an integer from the provided stream in little-endian form
int freadIntLSB(FILE *stream) {
    unsigned char array[4];

    fread(array, sizeof(unsigned char), 4, stream);
    
    int data;
    data = array[0] | (array[1] << 8) | (array[2] << 16) | (array[3] << 24);

    return data;
}


//writes a short integer to the provided stream in little-endian form
size_t fwriteShortLSB(short int data, FILE *stream) {
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}


//reads a short integer from the provided stream in little-endian form
short int freadShortLSB(FILE *stream) {
    unsigned char array[2];

    fread(array, sizeof(unsigned char), 2, stream);
    
    int data;
    data = array[0] | (array[1] << 8);

    return data;
}


