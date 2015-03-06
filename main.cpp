// Programmers: Majed Sahli, Tariq Alturkestani
// Date: January 2014 
// Main driver for APlug class 


#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string>

#include "APlug.h"


// for use with getopt(3) 
extern char *optarg;
extern int optind;
extern int optopt;
extern int opterr;

using namespace std;

static void usage(char *prog, int status)
{
	if (status == EXIT_SUCCESS)
	{
		printf("Usage: %s [-h] [-t n] [-s n] [-f n]\n",
		prog);
		printf("    -h      help\n");
		printf("    -t n    number of tasks\n");
		printf("    -s n    sample size\n");
		printf("    -f s    sample filename\n");
		printf("    -m s    miminum number of workers\n");
		printf("    -x s    maximum number of workers\n");
		printf("    -p s    steps to take from min to max workers\n");
	}
	else
	{
		fprintf(stderr, "%s: Try '%s -h' for usage information.\n", prog, prog);
	}

	exit(status);
}


int main(int argc, char**argv)
{
	APlug myPlug;					// Cloud Scalibility plug
	int ch;  // for use with getopt(3) 
	unsigned int numOfTasks = 0;
	unsigned int sampleSize = 0;
	unsigned int dataRead  = 0;
	int minWorkers = 1 ;
	int maxWorkers = MAX_NUM_WORKERS;
	int step = 15;
		
	double timeInFile; // used as a buf
	double sumTest =0;
	string sourceFilename ;
	bool isFilenamePassed = false;
	FILE * fd;
	
	while ((ch = getopt(argc, argv, ":ht:s:f:m:x:p:")) != -1)
	{
		switch (ch) 
		{
			case 'h':
				usage(argv[0], EXIT_SUCCESS);
				break;
			case 't':
				numOfTasks = strtol(optarg, NULL, 10);
				break;
			case 's':
				sampleSize = strtol(optarg, NULL, 10);
				break;
			case 'f':
				sourceFilename = string(optarg);
				isFilenamePassed = true;
				break;
			case 'm':
				minWorkers = atoi(optarg);
				break;
			case 'x':
				maxWorkers = atoi(optarg);
				break;
			case 'p':
				step = atoi(optarg);
				break;	
			case '?':
				printf("%s: invalid option '%c'\n", argv[0], optopt);
				usage(argv[0], EXIT_FAILURE);
				break;
			case ':':
				printf("%s: invalid option '%c' (missing argument)\n", argv[0], optopt);
				usage(argv[0], EXIT_FAILURE);
				break;
			default:
				usage(argv[0], EXIT_FAILURE);
				break;
		}
	}
	
	if ( !isFilenamePassed)
	{
		cout << "Source filename was not passed!\n";
		usage(argv[0], EXIT_FAILURE);
	}
	    
	// read file 
	fd = fopen(sourceFilename.c_str(),"r");
	if (!fd)
	{
		cout << "File: " << sourceFilename << " could not be opend" << endl;
		usage(argv[0], EXIT_FAILURE);
	}
	
	if (sampleSize == 0)
	{
		cout << "Sample size was not entered\n"
			 << "Sample size will be determined by the count of numbers in the source file\n";
					 
		while (fscanf(fd, "%lf ", &timeInFile) != EOF)
		{
			myPlug.pushTime(timeInFile);
			++sampleSize;
			sumTest +=timeInFile;
		}
	}
	else
	{
		while (fscanf(fd, "%lf ", &timeInFile) != EOF)
		{
			myPlug.pushTime(timeInFile);
			++dataRead;
			sumTest +=timeInFile;
		}
		
		if ( dataRead > sampleSize)
		{
			cout << "File: " << sourceFilename << " contains more data than expected\n"
				 << "Sample size is: " << sampleSize << " and the file contains: "
				 << dataRead << "\nSample size is reset" << endl;
			sampleSize = dataRead;
		}
		else if ( dataRead < sampleSize)
		{
			cout << "File: " << sourceFilename << " contains less data than expected\n"
				 << "Sample size is: " << sampleSize << " and the file contains: "
				 << dataRead << "\nSample size is reset" << endl;
			sampleSize = dataRead;
		}
	}
		
	fclose(fd);
	
	if ( numOfTasks == 0)
	{
		cout << "Number of tasks was not entered\n"
		     <<	"Number of tasks will be set to sample size\n";
		numOfTasks = sampleSize;	
	}
	else if ( numOfTasks < sampleSize )
    {
		cout << "Total number of tasks cannot be less than sample size!\n"
			 <<	"Number of tasks will be set to sample size\n";
		numOfTasks = sampleSize;	
	}	
	
	if (minWorkers < 0 || maxWorkers < 0 || step < 1)
	{
		cout << "Minimum/Maximum number of workers and step number can't be negative\n"
			<< "Setting them to default: min = 1, max = 5000, step = 15\n";
		minWorkers = 1;
		maxWorkers = MAX_NUM_WORKERS;
		step = 15;
		
	}  
	if (minWorkers >= maxWorkers)
	{
		cout << "Minimum number of workers can't be more than the maximum\n"
			 << "setting min to 1 and max to min...\n";
		maxWorkers = minWorkers;
		minWorkers = 1; 
		
	}
	
	
	// note the following lines would only be executed if all the variables are set
    cout << "Number of Tasks: " << numOfTasks << endl;
    cout << "Sample Size:     " << sampleSize << endl;
    cout << "Filename:        " << sourceFilename << endl;
    cout << "Getting recommendation for [" << minWorkers << ", " 
			<< (minWorkers == 1? step : minWorkers + step)<< ", " 
			<< (minWorkers == 1? 2*step : minWorkers + 2*step)<< ", " 
			<< (minWorkers == 1? 3*step : minWorkers + 3*step)<< ", " 
			<< ", ... , " << maxWorkers << ")\n" ;
    
	myPlug.setSampleSizeAndNumOfTasks(sampleSize, numOfTasks);
	myPlug.getRecommendation(minWorkers,maxWorkers, step );
	
	return 0;
}
