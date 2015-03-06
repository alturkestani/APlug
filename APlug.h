// Programmers: Majed Sahli, Tariq Alturkestani
// Date: January 2014 
// Adaptive Parallelism for Scientific Applications(APlug) class header 

#pragma once


#include <vector>
#include <deque>
#include <map>
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <stdio.h>

#define MAX_NUM_WORKERS 5000
#define PI 3.1415926535897932384626433832795
#define GAMMA_INTERVAL 10
using namespace std;

struct record
{
	unsigned int numOfWorkers;
	double tasksPerWorker;
	double parallelTime;
	double speedUpEfficiency;
	
	record()
	{
		numOfWorkers   = 0;
		tasksPerWorker = 0;
		parallelTime   = 0;
		speedUpEfficiency= 0;
	}
	
	record(unsigned int a,double b, 
		double c, double d)
	{
		numOfWorkers   = a;
		tasksPerWorker = b;
		parallelTime   = c;
		speedUpEfficiency= d;
	}
		 
};

class APlug
{
	
public:

	APlug(void);
	APlug(unsigned int numOfTasks);
	APlug(unsigned int numOfTasks, unsigned int sampleSize);
		
	void init(void);
	
	void setNumOfTasks(unsigned int numOfTasks);
	void setSampleSize(unsigned int sampleSize);
	void setSampleSizeAndNumOfTasks(unsigned int sampleSize, unsigned int numOfTasks);
	void setBucketsGranularity(unsigned int bucketsGranularity);
	void loadSamples(double sampleTimes[], int numOfSamples);
	double getSum(void);
	double getMean(void);
	double getStdev(void);
	double getMin(void);
	double getMax(void);
	double getEstSerialTimeSecGamma(void);
	double getEstSerialTimeSecHistogram(void);
	double getUpRate(void);
	
	unsigned int getBucketsGranularity(void);
	unsigned int getMaxWorkers(void);
	unsigned int getNumOfTasks(void);
	unsigned int getSampleSize(void);
	unsigned int getNumOfBuckets(void);
	unsigned int getRecommendedNumOfWorkers(void);
	
	int changeDecomposition(int last);
	void pushTime(double time);
	
	void getRecommendation(int minWorkers = 1, int maxWorkers = MAX_NUM_WORKERS, int step = 15); // prints the results 
	double getParallelTime(int numOfWorkers);
	double getSpeedUpEff(int numOfWorkers);
	
private:

	map<unsigned int , struct record> estimations;
	map<unsigned int , struct record> estimationsGamma;
	map<unsigned int , struct record>::iterator it;
	
	vector<double> sampleTaskTimes;
	vector<double> bucketsCount;
	vector<double> bucketsTime;
	vector<double> bucketsTimeAVG;
	vector<double>  estAllTasksTime;
	vector<double>  estAllTasksTimeGamma;
	
	unsigned int numOfTasks;
	unsigned int sampleSize;
	unsigned int numOfBuckets;
	unsigned int bucketsGranularity;
	double upRate; // numOfTasks/sampleSize
	double estSerialTimeSec;
	double estSerialTimeSecUsingGamma;
	double sum;
	double mean;
	double stdev;
	double min;
	double max;
	double recommendedNumOfWorkers;
	
	bool isSorted;
	bool isComputed;
	bool isStatSet;
	void setStats(void);
	void sortTimes(void);
	void buildHistogram(void);
	void makeTasks(void);
	void makeTasksUsingGamma(void);
	void speedupEstimationHistogram(int minWorkers, int maxWorkers, int step);	
	void speedupEstimationGamma(int minWorkers, int maxWorkers, int step);	
	void compute(int minWorkers, int maxWorkers, int step);
	double getTaskHistogram(void);
	double getTaskGamma(void);
	void recommendationByGammaToFile(void);
	void recommendationByHistogramToFile(void);
	void histogramToFile(void);
	
	// The following code was obtained from: 
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
	double incompletegamma(double a, double x);  // lower incomplete gamma distribution 
	double incompletegammac(double a, double x); // upper incomplete gamma distribution
	double gammafunction(double x);
	
	bool ae_fp_less(double v1, double v2);
	bool ae_fp_less_eq(double v1, double v2);
	bool ae_fp_greater(double v1, double v2);
	bool ae_fp_greater_eq(double v1, double v2);
	bool ae_fp_eq(double v1, double v2);
	double lngamma(double x, double* sgngam);
	double gammafunc_gammastirf(double x);
	
};
