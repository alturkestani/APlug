Programmers: Majed Sahli, Tariq Alturkestani
Date: January 2014 
Adaptive Parallelism for Scientific Applications

1 - Usage:
	A- If you already have a file of execution times of some subset, then you can follow the instruction in section 2.
	B- You can plug in the class into your code by including the header file APlug.h. See section 3.
	
2- Instruction on getting the recommendations when you have a subset of execution times
	A- 	First compile main.cpp, APlug.h and APlug.cpp using the attached makefile
		on the command prompt type: make 
	B-	Run: ./APlug  -t numOfTasks -s sampleSize -f filename -m mimNumOfWorkers -x maxNumOfWorkers -p steps
		where: 
			numOfTasks: 	 total number of tasks to execute
			sampleSize: 	 generally less that the total number of tasks. e.g. 10% of numOfTasks
			filename:   	 file that contains execution times of tasks is seconds. time per line. 
			mimNumOfWorkers: set the minimum number of works you want to get recommendation for.
			maxNumOfWorkers: set the maximum number of works you want to get recommendation for.
			steps:			 get recommendation for [mimNumOfWorkers, mimNumOfWorkers+ STEPS, mimNumOfWorkers+ 2*STEPS, ..., maxNumOfWorkers)
	C- Read output file (See section 4).
	
3- Instructions on plugging the class into your code
	A-	Copy APlug.h and APlug.cpp into the same directory of your code
	B-	Add #include <APlug.h> in you code
	C-    Use the following logic flow to get the recommendations
		int main ( ....) 
		{
			..
			..	
			APlug myPlug;
			..
			..
			myPlug.pushTime(num) // this operation will repeat sampleSize times
			..
			..
			
			myPlug.setSampleSizeAndNumOfTasks( sampleSize, numOfTasks);
			myPlug.getRecommendation(mimNumOfWorkers, maxNumOfWorkers, steps)  or myPlug.getRecommendation()
			// the default values of mimNumOfWorkers, maxNumOfWorkers, and steps are: 1, 5000, 15 
		}
	
	D- check the output files (Section 4)
		
4- Output files:
	A-	"recommendations.csv" 
		comma seperated values [NumberOfWorkers, TasksPerWorker, ParallelTime, SpeedupEfficiency]
		
	B-	"histogram.csv" 
		[BucketNumber, TasksInBucket, StartTimeOfBucket, MiddleOfBucket]
		note that you may also set the granularity of the buckets using the class member function setBucketsGranularity(num)


