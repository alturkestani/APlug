// Programmers: Majed Sahli, Tariq Alturkestani
// Date: January 2014 
// Adaptive Parallelism for Scientific Applications (APlug) class implementation 

#include "APlug.h"


APlug::APlug()
{
	init();
}

APlug::APlug(unsigned int sampleSize)
{
	init();
	this->sampleSize = sampleSize;
	this->numOfTasks = sampleSize;
	this->upRate = 1.0;
}

APlug::APlug(unsigned int numOfTasks, unsigned int sampleSize)
{
	init();
	this->numOfTasks = numOfTasks;
	this->sampleSize = sampleSize;
	this->upRate = (double)numOfTasks/(double)sampleSize;
}

void APlug::init()
{
	this->numOfTasks = 0;
	this->sampleSize = 0;
	this->sum   = 0;
	this->mean  = 0;
	this->stdev = 0;
	this->min   = 0;
	this->max   = 0;
	this->upRate= 0;
	this->isSorted = false;
	this->isComputed = false;
	this->isStatSet = false;
	this->recommendedNumOfWorkers = 0;
	this->numOfBuckets = 0;
	this->bucketsGranularity =2;
	this->estSerialTimeSec = 0.0;
}
	
// setters 		
void APlug::setBucketsGranularity(unsigned int bucketsGranularity)
{
	this->bucketsGranularity = bucketsGranularity;
}

void APlug::setNumOfTasks(unsigned int numOfTasks)
{
	this->numOfTasks = numOfTasks;
}

void APlug::setSampleSize(unsigned int sampleSize)
{
	this->sampleSize = sampleSize;
}

void APlug::setSampleSizeAndNumOfTasks(unsigned int sampleSize, unsigned int numOfTasks)
{
	this->numOfTasks = numOfTasks;
	this->sampleSize = sampleSize;
	this->upRate = (double)numOfTasks/(double)sampleSize;
};

void APlug::loadSamples(double sampleTimes[], int numOfSamples)
{
	for ( int i = 0; i < numOfSamples; ++i)
	{
		this->sampleTaskTimes.push_back(sampleTimes[i]);
	}
}

void APlug::setStats()
{
	if (this->isStatSet)
	{
		return;
	}
	else
	{
		this->isStatSet = true;
	}
	// sum	
	this->sum=0;
	for(unsigned int i=0;i < this->getSampleSize();i++) 
	{
		this->sum += this->sampleTaskTimes[i];
	}
	
	// mean
	this->mean = this->sum/this->sampleSize;

	//stdev
	double E=0;
	double inverse = 1.0 / static_cast<double>(this->sampleSize);
	for(unsigned int i=0;i<getSampleSize();i++) 
	{
		E += pow(static_cast<double>(sampleTaskTimes[i]) - mean, 2);
	}
	this->stdev = sqrt(inverse * E);
	
	// set min and max assumes that the sampleTaskTimes array is sorted
	if (isSorted)
	{
		this->min = this->sampleTaskTimes[this->sampleSize-1];
		this->max = this->sampleTaskTimes[0];
	}
	else // sort first (we need it sorted)
	{
		sort( sampleTaskTimes.begin(), sampleTaskTimes.end(), greater<double>());
		this->min = this->sampleTaskTimes[this->sampleSize-1];
		this->max = this->sampleTaskTimes[0];		
	}
	
}

// getters 
unsigned int APlug::getNumOfTasks()
{
	return this->numOfTasks;
}

unsigned int APlug::getSampleSize()
{
	return this->sampleSize;
}

unsigned int APlug::getMaxWorkers()
{
	return this->getNumOfTasks();
}
	
double APlug::getSum(void)
{
	if ( this->sum == 0 ) 
	{
		this->setStats();
	}
	return this->sum;
}

double APlug::getMean(void)
{
	if ( this->mean == 0 ) 
	{
		this->setStats();
	}
	return this->mean;
}

double APlug::getStdev(void)
{
	if ( this->stdev == 0 ) 
	{
		this->setStats();
	}
	return this->stdev;
}
	
double APlug::getMin(void)
{
	if ( this->min == 0 ) 
	{
		this->setStats();
	}
	return this->min;
}

double APlug::getMax(void)
{
	if ( this->max == 0 ) 
	{
		this->setStats();
	}
	return this->max;
}

double APlug::getUpRate(void)
{
	return this->upRate;
}

unsigned int APlug::getBucketsGranularity(void)
{	
	return this->bucketsGranularity;
}

double  APlug::getEstSerialTimeSecHistogram()
{
	return this->estSerialTimeSec;
}

double  APlug::getEstSerialTimeSecGamma()
{
	return this->estSerialTimeSecUsingGamma;
}

unsigned int APlug::getNumOfBuckets(void)
{
	return this->numOfBuckets;
}


unsigned int APlug::getRecommendedNumOfWorkers(void)
{
	
	return this->recommendedNumOfWorkers;
	
}
	
// other functions

int APlug::changeDecomposition(int last)
// users need to implement this function 
{
	cout << "changeDecomposition() is not implemented" << endl; 
	return last;
}


void APlug::sortTimes()
{
	sort( sampleTaskTimes.begin(), sampleTaskTimes.end(), greater<double>());
	sampleTaskTimes.resize(getSampleSize(),-1); // fill error with negative number (sanity check)
	this->isSorted = true;
}
 
void APlug::pushTime(double time)
{
	this->sampleTaskTimes.push_back(time);
}


void APlug::buildHistogram(void)
{
	this->numOfBuckets = ceil(getBucketsGranularity()*(getMax()-getMin())/getStdev());
	
	this->bucketsCount.resize(numOfBuckets,0);
	this->bucketsTime.resize(numOfBuckets,0);
	this->bucketsTimeAVG.resize(numOfBuckets,0);
	
	double rangeInBucket = (getMax()-getMin())/numOfBuckets;
	
	
	this->bucketsTime[0]= getMin();
	
	for(unsigned int i=1;i<this->numOfBuckets;i++)
	{
		bucketsTime[i]=bucketsTime[i-1]+rangeInBucket;
	}
		
	this->bucketsTimeAVG[0] = (this->bucketsTime[0] + this->bucketsTime[1] ) /2.0;
	
	for(unsigned int i=1;i<this->numOfBuckets;i++)
	{
		this->bucketsTimeAVG[i]=this->bucketsTimeAVG[i-1]+rangeInBucket;
	}
	
	
	for(unsigned int i=getSampleSize()-1; true ;--i)
	{
		for ( unsigned int j =0; j <this->numOfBuckets-1; j++)
		{
			if(this->sampleTaskTimes[i] >= bucketsTime[j] && this->sampleTaskTimes[i] < bucketsTime[j+1])
			{
				this->bucketsCount[j]++;
				break;
			}
		} 
		if( this->sampleTaskTimes[i] >=bucketsTime[this->numOfBuckets-1] && this->sampleTaskTimes[i] <= getMax())
		{
			this->bucketsCount[this->numOfBuckets-1]++;
		}
		if ( i==0)
		{
			break;
		}
	}
}

void APlug::makeTasks(void)
{
	
	for (unsigned int i = 0; i <this->getNumOfBuckets();++i)
	{
		for (unsigned int j = 0; j  < ceil(this->bucketsCount[i] * getUpRate()) ; ++j)
		{
			estAllTasksTime.push_back(this->bucketsTimeAVG[i]);
		}
	} 
	
	// set estimated serial time 
	for(unsigned int i=0;i<estAllTasksTime.size();i++)
	{
		this->estSerialTimeSec += estAllTasksTime[i];
	}
	
	// now randomize the vector
	srand(time(0));
	random_shuffle (estAllTasksTime.begin(), estAllTasksTime.end()); 
		
}

void APlug::makeTasksUsingGamma(void)
{
	//estAllTasksTimeGamma.push_back()
	
	//~ for (unsigned int i = 0; i <this->getNumOfBuckets();++i)
	//~ {
		//~ for (unsigned int j = 0; j  < ceil(this->bucketsCount[i] * getUpRate()) ; ++j)
		//~ {
			//~ estAllTasksTime.push_back(this->bucketsTimeAVG[i]);
		//~ }
	//~ } 
	
	double alpha;
	double beta;
	double maxRandom;
	double ratio;
	
	vector<double> rangeCDFs;
	
	
	//A = (m*m) / (std * std)
	//B = m / A
	//sample mean: 0.0722703
	//sample stdev: 0.0428602
	alpha = this->getMean()*this->getMean() / (this->getStdev() * this->getStdev());
	beta = (this->getMean()/alpha);
	
	volatile double cdf; 
	int i;
	for (  i = 0;  ; ++i)
	{
		cdf = incompletegamma(alpha,beta*i);
		if ( cdf >= 0.999999)
		{
			break; 
		}
	}
	maxRandom = i;
	ratio = this->getMax()/(double)i;
	
	FILE * gammaFD = fopen("GammaDistribution.csv","w");
	fprintf(gammaFD,"x = time (sec), P(X <= x) i.e. CDF\n");
	for (  i = 0; i < maxRandom  ; ++i)
	{
		cdf = incompletegamma(alpha,beta*i);
		fprintf(gammaFD,"%.15f, %.30f\n", ratio*i, cdf);
	}
	
	fclose(gammaFD);
	
	double mid;
	int localNumOfTasks = this->getNumOfTasks();
	if ( localNumOfTasks % 10 != 0)
	{
		localNumOfTasks +=  (10-(localNumOfTasks % 10 ));
	}
	
	//printf("num of tasks [%d]\n", localNumOfTasks);
	FILE * PDF = fopen("numOfTaskPerTimeGamma.csv", "w");
	fprintf(PDF, "numOfTasks, probability, totalTime\n");
	int j;
	int lnumOfTasks= 0;
	double ltime;
	for (  i = 0; i < maxRandom  ; i+= GAMMA_INTERVAL)
	{
		lnumOfTasks =0;
		mid = incompletegamma(alpha,beta*(i+ GAMMA_INTERVAL)) - incompletegamma(alpha,beta*i);
		//~ ltime = ratio*(i+GAMMA_INTERVAL/2.0);
		//~ ltime = ratio*(2*i+GAMMA_INTERVAL)/2.0;
		ltime = ratio*(i+GAMMA_INTERVAL);
		//ltime += this->getMin();
		for ( j = 0; j < ceil(mid*localNumOfTasks); ++j)
		{
			estAllTasksTimeGamma.push_back(ltime);
			++lnumOfTasks;
		}
		fprintf(PDF, "%d, %.20f, %.15f\n", lnumOfTasks, mid, ltime);
	}
	fclose(PDF);	
	// set estimated serial time 
	for(unsigned int i=0;i<estAllTasksTimeGamma.size();i++)
	{
		this->estSerialTimeSecUsingGamma += estAllTasksTimeGamma[i];
	}
	
	// now randomize the vector
	srand(time(0));
	random_shuffle (estAllTasksTimeGamma.begin(), estAllTasksTimeGamma.end()); 
		
}



struct record makeRecord(unsigned int numOfWorkers,	double tasksPerWorker,
						double parallelTime, double speedUpEfficiency)
{
	struct record aRecord;
	aRecord.numOfWorkers = numOfWorkers;
	aRecord.tasksPerWorker = tasksPerWorker;
	aRecord.parallelTime = parallelTime;
	aRecord.speedUpEfficiency = speedUpEfficiency;
	return aRecord;
}

double APlug::getTaskHistogram()
{
	static unsigned int i = 0;
	double ret;
	
	if ( i < this->estAllTasksTime.size())
	{
		ret = this->estAllTasksTime[i++];
	}
	else
	{
		i = 0; // counter reset 
		ret = -1;
	}
	
	return ret;	
}

double APlug::getTaskGamma()
{
	static unsigned int i = 0;
	double ret;
	
	if ( i < this->estAllTasksTimeGamma.size())
	{
		ret = this->estAllTasksTimeGamma[i++];
	}
	else
	{
		i = 0; // counter reset 
		ret = -1;
	}
	
	return ret;	
}

void APlug::speedupEstimationHistogram(int minWorkers , int maxWorkers, int step )
{
	this->recommendedNumOfWorkers = 0;
	double mostEfficient = 0;
	double efficiency;
	double parallelTime;
	//~ for ( int workers = 50; workers < MAX_NUM_WORKERS; workers+=50)
	//~ for ( int workers = 1; workers < MAX_NUM_WORKERS; ((workers != 1)? workers *= 2: workers = step))  
        for ( int workers = minWorkers; workers <= maxWorkers; workers += step)
	{		
		double core[workers];
		double coreSum[workers];
		
		std::fill_n(core, workers, 0);
		std::fill_n(coreSum, workers, 0);
		cout << "Speedup estimation for: " << workers << " workers" << endl;
		
		double curWork = 0; 
		while(true)
		{
			curWork = getTaskHistogram();
			if ( curWork != -1) // if there is a task to schedule then go ahead and schedule it
			{
				int i;
				for ( i = 0 ; i < workers ; i++) // look for an empty core
				{
					if ( core[i] == 0) // empty core has been found 
					{
						break;
					}
				}
				
				if ( i < workers) //if there is an empty slot then schedule 
				{
					core[i] = curWork;
					coreSum[i] += curWork;
				}
				else // if not then clean up the cores and find a new place for curTask
				{
					i = 0;
					double min = core[i];
					for ( int j = 0 ; j < workers ; j++)
					{
						if ( core[j] < min)
						{
							min = core[j];
							i = j;  // index for the new empty slot
						}
					}
										
					for ( int j = 0 ; j < workers ; j++)
					{
						core[j] -= min;
					}
					
					core[i] = curWork;
					coreSum[i] += curWork;
				}
			}
			else // no more job? what are you doing here!! 
			{
				break;
			}
		};
		
		// find parallel time by finding the maximum time among cores 
		parallelTime = coreSum[0];
		for ( int i = 0; i < workers ; ++i)
		{
			if ( coreSum[i] > parallelTime)
			{
				parallelTime = coreSum[i];
			}
		}
		
		
		efficiency = getEstSerialTimeSecHistogram()/(parallelTime*workers);
		this->estimations[workers] = makeRecord(workers, getNumOfTasks()/(double)workers, 
				parallelTime,efficiency);
				
		if ( efficiency > mostEfficient)
		{
			mostEfficient = efficiency;
			this->recommendedNumOfWorkers = workers;
		} 		
	
	
	}

}

void APlug::speedupEstimationGamma(int minWorkers , int maxWorkers, int step )
{
	this->recommendedNumOfWorkers = 0;
	double mostEfficient = 0;
	double efficiency;
	double parallelTime;
	//~ for ( int workers = 50; workers < MAX_NUM_WORKERS; workers+=50)
	//~ for ( int workers = 1; workers < MAX_NUM_WORKERS; ((workers != 1)? workers *= 2: workers = step))  
        for ( int workers = minWorkers; workers <= maxWorkers; workers += step)
	{		
		double core[workers];
		double coreSum[workers];
		
		std::fill_n(core, workers, 0);
		std::fill_n(coreSum, workers, 0);
		cout << "Speedup estimation for: " << workers << " workers" << endl;
		
		double curWork = 0; 
		while(true)
		{
			curWork = getTaskGamma();
			if ( curWork != -1) // if there is a task to schedule then go ahead and schedule it
			{
				int i;
				for ( i = 0 ; i < workers ; i++) // look for an empty core
				{
					if ( core[i] == 0) // empty core has been found 
					{
						break;
					}
				}
				
				if ( i < workers) //if there is an empty slot then schedule 
				{
					core[i] = curWork;
					coreSum[i] += curWork;
				}
				else // if not then clean up the cores and find a new place for curTask
				{
					i = 0;
					double min = core[i];
					for ( int j = 0 ; j < workers ; j++)
					{
						if ( core[j] < min)
						{
							min = core[j];
							i = j;  // index for the new empty slot
						}
					}
										
					for ( int j = 0 ; j < workers ; j++)
					{
						core[j] -= min;
					}
					
					core[i] = curWork;
					coreSum[i] += curWork;
				}
			}
			else // no more job? what are you doing here!! 
			{
				break;
			}
		};
		
		// find parallel time by finding the maximum time among cores 
		parallelTime = coreSum[0];
		for ( int i = 0; i < workers ; ++i)
		{
			if ( coreSum[i] > parallelTime)
			{
				parallelTime = coreSum[i];
			}
		}
		
		
		efficiency = getEstSerialTimeSecGamma()/(parallelTime*workers);
		this->estimationsGamma[workers] = makeRecord(workers, getNumOfTasks()/(double)workers, 
				parallelTime,efficiency);
				
		if ( efficiency > mostEfficient)
		{
			mostEfficient = efficiency;
			this->recommendedNumOfWorkers = workers;
		} 		
		
	}

}

void APlug::recommendationByHistogramToFile(void)
{
	FILE * fd;
	fd = fopen("recommendationsByHistogram.csv", "w");
	if (!fd)
	{
		cout << "File not created!\n";
	}
	
	fprintf(fd, "NumberOfWorkers, TasksPerWorker, ParallelTime, SpeedupEfficiency\n");
	
	for ( this->it= estimations.begin(); it!=estimations.end(); ++it)
	{
		fprintf(fd, "%u, %lf, %lf, %lf\n", it->first, it->second.tasksPerWorker, 
													it->second.parallelTime,
													it->second.speedUpEfficiency);
	}	
	
	fclose(fd);
}

void APlug::recommendationByGammaToFile(void)
{
	FILE * fd;
	fd = fopen("recommendationsByGamma.csv", "w");
	if (!fd)
	{
		cout << "File not created!\n";
	}
	
	fprintf(fd, "NumberOfWorkers, TasksPerWorker, ParallelTime, SpeedupEfficiency\n");
	
	for ( this->it= estimationsGamma.begin(); it!=estimationsGamma.end(); ++it)
	{
		fprintf(fd, "%u, %lf, %lf, %lf\n", it->first, it->second.tasksPerWorker, 
													it->second.parallelTime,
													it->second.speedUpEfficiency);
	}	
	
	fclose(fd);
}

void APlug::histogramToFile(void)
{
	FILE * fd;
	fd = fopen("histogram.csv", "w");
	if (!fd)
	{
		cout << "File not created!\n";
	}
	
	fprintf(fd, "BucketNumber, TasksInBucket, StartTimeOfBucket, MiddleOfBucket\n");
	
	for (unsigned int i = 0; i< this->getNumOfBuckets() ;++i)
	{
		fprintf(fd, "%d, %lf, %lf, %lf\n", i, this->bucketsCount[i],
										this->bucketsTime[i], this->bucketsTimeAVG[i]);
	}	
	
	fclose(fd);
}

void APlug::compute(int minWorkers  , int maxWorkers , int step)
{
	this->isComputed = true;
	
	sortTimes();
	cout << "Sample times are sorted\n";
	// set sum, mean, avg, min and max
	setStats();
	cout << "Stats are set\n";
	// build histogram
	buildHistogram();
	cout << "Histogram is built\n";
	// making taks array using historgam and making sure to shuffle it 
	makeTasks();
	cout << "Tasks are made using the historgam\n";
	// making taks array and making sure to shuffle it 
	makeTasksUsingGamma();
	cout << "Tasks are made using Gamma distribution\n";
	// speedup estimation
	speedupEstimationHistogram(minWorkers, maxWorkers, step);
	cout << "Speedup estimation using histogram is done\n";
	// speedup estimation
	speedupEstimationGamma(minWorkers, maxWorkers, step);
	cout << "Speedup estimation using gamma is done\n";
	// write to csv file
	recommendationByHistogramToFile();
	cout << "Recommendation using historgram written to a file\n";
	// write to csv file
	recommendationByGammaToFile();
	cout << "Recommendation using gamma written to a file\n";
	histogramToFile();
	cout << "Histograms written to a file\n";
	
}

void APlug::getRecommendation(int minWorkers  , int maxWorkers , int step )
{
	
	if ( !this->isComputed)
	{
		this->compute(minWorkers, maxWorkers, step);
	}
	
	cout << "Number of tasks:  " << getNumOfTasks() << endl;
    cout << "Sample size:      " << getSampleSize() << endl;
	cout << "sample sum:       " << getSum()  << endl;
	cout << "sample min:       " << getMin()  << endl;
	cout << "sample max:       " << getMax()  << endl;
	cout << "sample mean:      " << getMean() << endl;
	cout << "sample stdev:     " << getStdev()<< endl;
	cout << "Estimated serial time (histogram): " << getEstSerialTimeSecHistogram() << " sec" << endl;	
	cout << "Estimated serial time (gamma): " << getEstSerialTimeSecGamma() << " sec" << endl;	
	cout << "Suggested max workers: " << getMaxWorkers() << endl; 
	cout << "Recommended number of workers: " << getRecommendedNumOfWorkers() << endl;
	cout << "You can also find/plot other options in: \n"
	     << "\"recommendationsByHistogram.csv\" \"recommendationsByGamma.csv\" \n\"histogram.csv\" \"GammaDistribution.csv\" " << endl;
			
			
}

/*************************************************************************
Copyright (c) Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/
///////////////////////////////////////
double APlug::getParallelTime(int numOfWorkers)
{
	
	if ( this->estimations.find(numOfWorkers) == this->estimations.end() ) 
	{
		this->compute(numOfWorkers, numOfWorkers, 1);
	}
	
	return this->estimations[numOfWorkers].parallelTime;
}

double APlug::getSpeedUpEff(int numOfWorkers)
{
	if ( this->estimations.find(numOfWorkers) == this->estimations.end() ) 
	{
		this->compute(numOfWorkers, numOfWorkers, 1);
	}
	
	return this->estimations[numOfWorkers].speedUpEfficiency;
	
}

bool APlug::ae_fp_less(double v1, double v2)
{
    /* IEEE-strict floating point comparison */
    volatile double x = v1;
    volatile double y = v2;
    return x<y;
}

bool APlug::ae_fp_less_eq(double v1, double v2)
{
    /* IEEE-strict floating point comparison */
    volatile double x = v1;
    volatile double y = v2;
    return x<=y;
}

bool APlug::ae_fp_greater(double v1, double v2)
{
    /* IEEE-strict floating point comparison */
    volatile double x = v1;
    volatile double y = v2;
    return x>y;
}

bool APlug::ae_fp_greater_eq(double v1, double v2)
{
    /* IEEE-strict floating point comparison */
    volatile double x = v1;
    volatile double y = v2;
    return x>=y;
}
bool APlug::ae_fp_eq(double v1, double v2)
{
    /* IEEE-strict floating point comparison */
    volatile double x = v1;
    volatile double y = v2;
    return x==y;
}


/*************************************************************************
Natural logarithm of gamma function

Input parameters:
    X       -   argument

Result:
    logarithm of the absolute value of the Gamma(X).

Output parameters:
    SgnGam  -   sign(Gamma(X))

Domain:
    0 < X < 2.55e305
    -2.55e305 < X < 0, X is not an integer.

ACCURACY:
arithmetic      domain        # trials     peak         rms
   IEEE    0, 3                 28000     5.4e-16     1.1e-16
   IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
The error criterion was relative when the function magnitude
was greater than one but absolute when it was less than one.

The following test used the relative error criterion, though
at certain points the relative error could be much higher than
indicated.
   IEEE    -200, -4             10000     4.8e-16     1.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
*************************************************************************/
double APlug::lngamma(double x, double* sgngam)
{

    double a;
    double b;
    double c;
    double p;
    double q;
    double u;
    double w;
    double z;
    long long int i;
    double logpi;
    double ls2pi;
    double tmp;
    double result;

    *sgngam = 0;

    *sgngam = 1;
    logpi = 1.14472988584940017414;
    ls2pi = 0.91893853320467274178;
    if( ae_fp_less(x,-34.0) )
    {
        q = -x;
        w = lngamma(q, &tmp);
        p = floor(q);
        i = floor(p+.5);
        if( i%2==0 )
        {
            *sgngam = -1;
        }
        else
        {
            *sgngam = 1;
        }
        z = q-p;
        if( ae_fp_greater(z,0.5) )
        {
            p = p+1;
            z = p-q;
        }
        z = q*sin(PI*z);
        result = logpi-log(z)-w;
        return result;
    }
    if( ae_fp_less(x,13) )
    {
        z = 1;
        p = 0;
        u = x;
        while(ae_fp_greater_eq(u,3))
        {
            p = p-1;
            u = x+p;
            z = z*u;
        }
        while(ae_fp_less(u,2))
        {
            z = z/u;
            p = p+1;
            u = x+p;
        }
        if( ae_fp_less(z,0) )
        {
            *sgngam = -1;
            z = -z;
        }
        else
        {
            *sgngam = 1;
        }
        if( ae_fp_eq(u,2) )
        {
            result = log(z);
            return result;
        }
        p = p-2;
        x = x+p;
        b = -1378.25152569120859100;
        b = -38801.6315134637840924+x*b;
        b = -331612.992738871184744+x*b;
        b = -1162370.97492762307383+x*b;
        b = -1721737.00820839662146+x*b;
        b = -853555.664245765465627+x*b;
        c = 1;
        c = -351.815701436523470549+x*c;
        c = -17064.2106651881159223+x*c;
        c = -220528.590553854454839+x*c;
        c = -1139334.44367982507207+x*c;
        c = -2532523.07177582951285+x*c;
        c = -2018891.41433532773231+x*c;
        p = x*b/c;
        result = log(z)+p;
        return result;
    }
    q = (x-0.5)*log(x)-x+ls2pi;
    if( ae_fp_greater(x,100000000) )
    {
        result = q;
        return result;
    }
    p = 1/(x*x);
    if( ae_fp_greater_eq(x,1000.0) )
    {
        q = q+((7.9365079365079365079365*0.0001*p-2.7777777777777777777778*0.001)*p+0.0833333333333333333333)/x;
    }
    else
    {
        a = 8.11614167470508450300*0.0001;
        a = -5.95061904284301438324*0.0001+p*a;
        a = 7.93650340457716943945*0.0001+p*a;
        a = -2.77777777730099687205*0.001+p*a;
        a = 8.33333333333331927722*0.01+p*a;
        q = q+a/x;
    }
    result = q;
    return result;

}



double APlug::gammafunc_gammastirf(double x)
{
    double y;
    double w;
    double v;
    double stir;
    double result;


    w = 1/x;
    stir = 7.87311395793093628397E-4;
    stir = -2.29549961613378126380E-4+w*stir;
    stir = -2.68132617805781232825E-3+w*stir;
    stir = 3.47222221605458667310E-3+w*stir;
    stir = 8.33333333333482257126E-2+w*stir;
    w = 1+w*stir;
    y = exp(x);
    if( ae_fp_greater(x,143.01608) )
    {
        v = pow(x, 0.5*x-0.25);
        y = v*(v/y);
    }
    else
    {
        y = pow(x, x-0.5)/y;
    }
    result = 2.50662827463100050242*y*w;
    return result;
}


/*************************************************************************
Gamma function

Input parameters:
    X   -   argument

Domain:
    0 < X < 171.6
    -170 < X < 0, X is not an integer.

Relative error:
 arithmetic   domain     # trials      peak         rms
    IEEE    -170,-33      20000       2.3e-15     3.3e-16
    IEEE     -33,  33     20000       9.4e-16     2.2e-16
    IEEE      33, 171.6   20000       2.3e-15     3.2e-16

Cephes Math Library Release 2.8:  June, 2000
Original copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
*************************************************************************/
double APlug::gammafunction(double x)
{

    double p;
    double pp;
    double q;
    double qq;
    double z;
    long long int i;
    double sgngam;
    double result;

    sgngam = 1;
    q = fabs(x);
    if( ae_fp_greater(q,33.0) )
    {
        if( ae_fp_less(x,0.0) )
        {
            p = floor(q);
            i = floor(p+.5);
            if( i%2==0 )
            {
                sgngam = -1;
            }
            z = q-p;
            if( ae_fp_greater(z,0.5) )
            {
                p = p+1;
                z = q-p;
            }
            z = q*sin(PI*z);
            z = fabs(z);
            z = PI/(z*gammafunc_gammastirf(q));
        }
        else
        {
            z = gammafunc_gammastirf(x);
        }
        result = sgngam*z;
        return result;
    }
    z = 1;
    while(ae_fp_greater_eq(x,3))
    {
        x = x-1;
        z = z*x;
    }
    while(ae_fp_less(x,0))
    {
        if( ae_fp_greater(x,-0.000000001) )
        {
            result = z/((1+0.5772156649015329*x)*x);
            return result;
        }
        z = z/x;
        x = x+1;
    }
    while(ae_fp_less(x,2))
    {
        if( ae_fp_less(x,0.000000001) )
        {
            result = z/((1+0.5772156649015329*x)*x);
            return result;
        }
        z = z/x;
        x = x+1.0;
    }
    if( ae_fp_eq(x,2) )
    {
        result = z;
        return result;
    }
    x = x-2.0;
    pp = 1.60119522476751861407E-4;
    pp = 1.19135147006586384913E-3+x*pp;
    pp = 1.04213797561761569935E-2+x*pp;
    pp = 4.76367800457137231464E-2+x*pp;
    pp = 2.07448227648435975150E-1+x*pp;
    pp = 4.94214826801497100753E-1+x*pp;
    pp = 9.99999999999999996796E-1+x*pp;
    qq = -2.31581873324120129819E-5;
    qq = 5.39605580493303397842E-4+x*qq;
    qq = -4.45641913851797240494E-3+x*qq;
    qq = 1.18139785222060435552E-2+x*qq;
    qq = 3.58236398605498653373E-2+x*qq;
    qq = -2.34591795718243348568E-1+x*qq;
    qq = 7.14304917030273074085E-2+x*qq;
    qq = 1.00000000000000000320+x*qq;
    result = z*pp/qq;
    return result;

}


/*************************************************************************
Complemented incomplete gamma integral

The function is defined by


 igamc(a,x)   =   1 - igam(a,x)

                           inf.
                             -
                    1       | |  -t  a-1
              =   -----     |   e   t   dt.
                   -      | |
                  | (a)    -
                            x


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

Tested at random a, x.
               a         x                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double APlug::incompletegammac(double a, double x)
{
    double igammaepsilon;
    double igammabignumber;
    double igammabignumberinv;
    double ans;
    double ax;
    double c;
    double yc;
    double r;
    double t;
    double y;
    double z;
    double pk;
    double pkm1;
    double pkm2;
    double qk;
    double qkm1;
    double qkm2;
    double tmp;
    double result;


    igammaepsilon = 0.000000000000001;
    igammabignumber = 4503599627370496.0;
    igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
    if( ae_fp_less_eq(x,0)||ae_fp_less_eq(a,0) )
    {
        result = 1;
        return result;
    }
    if( ae_fp_less(x,1)||ae_fp_less(x,a) )
    {
        result = 1-incompletegamma(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, &tmp);
    if( ae_fp_less(ax,-709.78271289338399) )
    {
        result = 0;
        return result;
    }
    ax = exp(ax);
    y = 1-a;
    z = x+y+1;
    c = 0;
    pkm2 = 1;
    qkm2 = x;
    pkm1 = x+1;
    qkm1 = z*x;
    ans = pkm1/qkm1;
    do
    {
        c = c+1;
        y = y+1;
        z = z+2;
        yc = y*c;
        pk = pkm1*z-pkm2*yc;
        qk = qkm1*z-qkm2*yc;
        if( !ae_fp_eq(qk,0) )
        {
            r = pk/qk;
            t = abs((ans-r)/r);
            ans = r;
        }
        else
        {
            t = 1;
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if( ae_fp_greater(fabs(pk),igammabignumber) )
        {
            pkm2 = pkm2*igammabignumberinv;
            pkm1 = pkm1*igammabignumberinv;
            qkm2 = qkm2*igammabignumberinv;
            qkm1 = qkm1*igammabignumberinv;
        }
    }
    while(ae_fp_greater(t,igammaepsilon));
    result = ans*ax;
    return result;
}


/*************************************************************************
Incomplete gamma integral

The function is defined by

                          x
                           -
                  1       | |  -t  a-1
 igam(a,x)  =   -----     |   e   t   dt.
                 -      | |
                | (a)    -
                          0


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30       200000       3.6e-14     2.9e-15
   IEEE      0,100      300000       9.9e-14     1.5e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double APlug::incompletegamma(double a, double x)
{
    double igammaepsilon;
    double ans;
    double ax;
    double c;
    double r;
    double tmp;
    double result;


    igammaepsilon = 0.000000000000001;
    if( ae_fp_less_eq(x,0)||ae_fp_less_eq(a,0) )
    {
        result = 0;
        return result;
    }
    if( ae_fp_greater(x,1)&&ae_fp_greater(x,a) )
    {
        result = 1-incompletegammac(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, &tmp);
    if( ae_fp_less(ax,-709.78271289338399) )
    {
        result = 0;
        return result;
    }
    ax = exp(ax);
    r = a;
    c = 1;
    ans = 1;
    do
    {
        r = r+1;
        c = c*x/r;
        ans = ans+c;
    }
    while(ae_fp_greater(c/ans,igammaepsilon));
    result = ans*ax/a;
    return result;
}
