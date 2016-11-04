#define _USE_MATH_DEFINES
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>

/***************************************************************************
Richard Clayton, University of SHeffield, March 2016  
Modified by 
	Eugene Chang, University of Sheffield and 
	Yen Ting Lin, The University of Manchester
Implemented in C/C++ by Yen Ting Lin

This script implements a simple Moe type cellular automaton that could be 
used for simulating aF over long time periods.
***************************************************************************/

using namespace std;

struct parameter
{
	int refractoryPeriod;
	int neighbourThreshold;	
	int numberParticles;
	int varRP;
	int numTimeSteps;
	int sinusPeriod;
	int snNumber;
	int * sinusNode;
	int nbhdMaxNumber;
	int ** nbhdMap;
	int pvActivationCandidateNumber;
	int pvActivationMaxNumber;
	int ** pvActivationMap;
	int initialDefectNumber;
	double * initialDefectProb;
	double pvBurstOnToOff;
	double pvBurstOffToOn;
	double pvBurstActivationProb;
	double restitutionA;
	double restitutionB;
	double restitutionTau;
	
	double detailedRecordStart;
	double detailedRecordEnd;
	int detailedRecordDt;

	double coarseRecordStart;
	double coarseRecordEnd;
	int coarseRecordDt;
	
};

struct state
{
    int* currentState;
    int* nextState;
    int* epsilonRefractory;
    int* refractory;
    int* refractory1;
    int* actTime;
    int* cycleLength;
    int* diastolInt;
	int* cellType;				// a label to specify the type of the cell. 0: normal, 1: inexcitable, 2: PV
	int* activeNeighbours;
	int nextBeat;
	bool pvBurst;
};

struct twin
{
	double n1, n2;
};

double rndexclusive();
twin normal();

int main(int argc, char *argv[])
{

	int ensN = 1;		// number of sample paths to be generated
	
	unsigned RS1;
	
	double newState;
	twin temp_normal;
	
	clock_t begin = clock();

	double maxPaceBeatNum = 4;    // duration of the PV burst
	int maxWaitBeatNum = 10;	  // duration of observation window
	
	parameter par;
	state var; 
	
	int temp_var, temp_var2;
	
	/**************************/
	/*initiation of parameter0s*/
	/**************************/
	
    par.refractoryPeriod = 120;   	// number of time steps for a cell to recover (default 25)
    par.neighbourThreshold = 8;     // number of active neigbours
	par.initialDefectNumber = 300;
    par.varRP = 8;               	// variability in refractory period (2)
    par.sinusPeriod = 400;       	// number of timesteps between regular beats
	par.numTimeSteps = 20*par.sinusPeriod; 		// duration of simulation
    par.pvBurstOnToOff= 0;
    par.pvBurstOffToOn= 0;;
	par.pvBurstActivationProb= 0.05;
	par.restitutionA = 121;
    par.restitutionB = 1.0;
    par.restitutionTau = 40.0;
    
	
	// Here set up the window for recording "coarse variables"
	par.coarseRecordStart = 0;
	par.coarseRecordEnd = (maxPaceBeatNum+maxWaitBeatNum)*par.sinusPeriod;
	par.coarseRecordDt = 1;
	
	// Here set up the window for detailed output (i.e., every point)
    par.detailedRecordStart = 4001;;
    par.detailedRecordEnd = 4400;;
	par.detailedRecordDt = 1;
    

	// Reading the geometry, topology, and initiation of the variables
	
	ifstream input;
	
	input.open("geometry.txt");
	input >> par.numberParticles;
	
	var.currentState = new int  [par.numberParticles];
	var.nextState = new int  [par.numberParticles];
	var.epsilonRefractory = new int  [par.numberParticles];
	var.refractory = new int  [par.numberParticles];
    var.refractory1 = new int  [par.numberParticles];
	var.actTime = new int  [par.numberParticles];
	var.cycleLength = new int  [par.numberParticles];
    var.diastolInt = new int  [par.numberParticles];
	var.cellType = new int [par.numberParticles];
	var.activeNeighbours = new int [par.numberParticles];
	par.initialDefectProb = new double [par.numberParticles];
	
	int last_0 = 0;
	int inAF = 0;
	int threshold = 10*par.sinusPeriod;
	
	int registryLength = (int)threshold+500;
	int ** tempRegistry = new int* [registryLength];
	
	for (int i=0;i<registryLength;i++)
	{
		tempRegistry[i] = new int [par.numberParticles];
	}
	
	
	
	for (int i=0;i<par.numberParticles;i++)
	{
		
		var.currentState[i] = 0;
		
		var.nextState[i] = 0;
		var.epsilonRefractory[i] = (int)((double)par.varRP * rndexclusive());
		var.refractory[i] = par.refractoryPeriod + var.epsilonRefractory[i];
        var.refractory1[i] = par.refractoryPeriod + var.epsilonRefractory[i];
		var.actTime[i] = 0;
		var.cycleLength[i] = par.sinusPeriod;
        var.diastolInt[i] = par.sinusPeriod;
		var.cellType[i] = 0;
		var.activeNeighbours[i] = 0;
		var.pvBurst = true;
		input >> par.initialDefectProb[i];
		
	}

	
	input >> par.pvActivationCandidateNumber;
	input >> par.pvActivationMaxNumber;
	
	par.pvActivationMap = new int * [par.pvActivationCandidateNumber];
	
	for (int i=0;i<par.pvActivationCandidateNumber;i++)
	{
		par.pvActivationMap[i] = new int [par.pvActivationMaxNumber];
		
		for (int j=0;j<par.pvActivationMaxNumber;j++)
		{
			input >> par.pvActivationMap[i][j];
		}
	}
	
	input >> par.nbhdMaxNumber;
	
	par.nbhdMap = new int * [par.numberParticles];
	
	for (int i=0;i<par.numberParticles;i++)
	{
		par.nbhdMap[i] = new int [par.nbhdMaxNumber];
		
		for (int j=0;j<par.nbhdMaxNumber;j++)
		{
			input >> par.nbhdMap[i][j];
		}
		
	}
	
	input >> par.snNumber;
	
	par.sinusNode = new int [par.snNumber];

	for (int i=0;i<par.snNumber;i++)
	{
		input >> par.sinusNode[i];
	}
	
	
	
	input.close();
	
	
	
	
	
	/************************/
	/* declare output files */
	/************************/
	
	int Success = 0;
	int ensIndex;
		
	ofstream output;
	stringstream filename;
	filename << "coarseEvolution.txt";
	output.open(filename.str().c_str());
				
	ofstream output2;
	stringstream filename2;
	filename2 << "detailedEvolution.txt";
	output2.open(filename2.str().c_str());
	
	for (ensIndex = 0; ensIndex < ensN ; ensIndex ++)
	{
		
        inAF = 0;
        last_0 = 0;
        
		cout << "currently running " << ensIndex << "success = " << Success << ", success rate = " << (double)Success/(ensIndex+1) << endl;

		RS1 = 1357246+ensIndex;	// this random seed does not induce reentry
		RS1 = 1357337+ensIndex; // this random seed induces reentry
		
		srand(RS1);	
				
		for (int i=0;i<par.numberParticles;i++)
		{
			var.currentState[i] = 0;
			var.nextState[i] = 0;
			var.epsilonRefractory[i] = (int)((double)par.varRP * rndexclusive());
			var.refractory[i] = par.refractoryPeriod + var.epsilonRefractory[i];
			var.refractory1[i] = var.refractory[i];
			var.actTime[i] = 0;
			var.cycleLength[i] = par.sinusPeriod;
			var.diastolInt[i] = par.sinusPeriod;
			var.cellType[i] = 0;
			var.activeNeighbours[i] = 0;
		}
		var.pvBurst = true;
			
		
		/******************************************************************************/
		/* 	Specify initial state of age-related remodelling. This sets up region of  */
		/* inexcitable cells that are assumed to result from age-related              */
		/* remodelling, and provide initial vulnerability to AF. This can be changed  */
		/* in various ways to alter vulnerability. Looping to 1000 (or reducing the   */
		/* SD of the remodelled region) will guarantee persistent AF arrives earlier. */
		/******************************************************************************/
		
		int counting = 0;
			
		while (counting < par.initialDefectNumber)
		{
			
			// randomly pick a point according to probability of having initial defect
			// using inverse sampling method
			
			double r = rndexclusive();
			int index = 0;
			double tempr = 0.0;
			
			while (tempr<r)
			{
				tempr += par.initialDefectProb[index];
				index ++;
			}
			
			index --;
			
			if (var.cellType[index] == 0)
			{
				counting ++;
				var.cellType[index] = 1;
				var.currentState[index] = -10;
				var.nextState[index]= -10;
			}
			
		}
		
		/*************************/
		/* main loop starts here */
		/*************************/
		
		
		var.nextBeat = 0;
		
		for (int t = 0;t < (int)(maxPaceBeatNum*par.sinusPeriod);t++)
		{
			
			if (t%par.sinusPeriod==0)
				cout << "Executing beat =" << t/par.sinusPeriod << "PV bursts: ON" << endl;
			
			if ((t%par.detailedRecordDt==0)&&(t>=par.detailedRecordStart)&&(t<=par.detailedRecordEnd))
				output2 << t << "\t";
			
			
			for (int i=0;i<par.numberParticles;i++)
			{
				tempRegistry[t%registryLength][i] = var.currentState[i];
			}
			
			if (t==var.nextBeat)
			{
				var.nextBeat += par.sinusPeriod;
				
				for (int i=0;i<par.snNumber;i++)
				{
					var.currentState[par.sinusNode[i]] = var.refractory[par.sinusNode[i]];
					
					if (var.actTime[par.sinusNode[i]] > 0)
					{
						var.cycleLength[par.sinusNode[i]] = t - var.actTime[par.sinusNode[i]];
						var.diastolInt[par.sinusNode[i]] = var.cycleLength[par.sinusNode[i]]-var.refractory[par.sinusNode[i]];
					}
					var.actTime[par.sinusNode[i]] = t;
				}
			
			}
			

			if ((var.pvBurst) && (rndexclusive() < par.pvBurstActivationProb))
			{
				// excite one set of cells (with uniform probability)
				
				int index = (int) (rndexclusive() * par.pvActivationCandidateNumber);
				
				for (int target = 1;target <= par.pvActivationMap[index][0];target++)
				{		
					int index2 = par.pvActivationMap[index][target];
					
					if (var.cellType[index2]==0)
					{
						if (var.actTime[index2] > 0)
						{
							var.cycleLength[index2] = t - var.actTime[index2];
							var.diastolInt[index2] = var.cycleLength[index2]-var.refractory[index2];
						}
								
						var.actTime[index2] = t;
						
						
						newState = max(floor((double)var.refractory1[index2]/2) + 5, floor(par.restitutionA*(1-par.restitutionB*exp( -(double)var.diastolInt[index2]/par.restitutionTau)) ) );
						//newState = max(floor((double)var.refractory[index2]/2) + 5, floor(sqrt(par.restitutionA*(double)var.cycleLength[index2]) + par.refractoryPeriod- sqrt(par.restitutionA*(double)par.sinusPeriod)+ (double)var.epsilonRefractory[index2]));
						//newState = floor(sqrt(8*(double)var.cycleLength[index2]) + 10 + (double)var.epsilonRefractory[index2]);
			
						var.refractory[index2] = newState;
						var.currentState[index2] = var.refractory[index2];
						
					}					
				}
				
			}
			
			/****************************************************************************/		
			/* determine next state of each "cell" based on number of active cells in   */
			/* % neighbourhood															*/
			/****************************************************************************/
			
			for (int i=0;i<par.numberParticles;i++)
			{
				
				
					if (var.cellType[i] == 0)
					{
						if (var.currentState[i] == 0)
						{				    
							if (var.activeNeighbours[i] >= (double)par.neighbourThreshold)
							{
								if (var.actTime[i] > 0)
								{
									var.cycleLength[i] = t - var.actTime[i];
									var.diastolInt[i] = var.cycleLength[i]-var.refractory[i];
								}
								
								var.actTime[i] = t;
								
								newState = max(floor((double)var.refractory1[i]/2) + 5, floor(par.restitutionA*(1-par.restitutionB*exp( -(double)var.diastolInt[i]/par.restitutionTau)) ) );
								//newState = max(floor((double)var.refractory[i]/2) + 5, floor(sqrt(par.restitutionA*(double)var.cycleLength[i]) +  par.refractoryPeriod- sqrt(par.restitutionA*(double)par.sinusPeriod) + (double)var.epsilonRefractory[i]));
								//newState = floor(sqrt(8*(double)var.cycleLength[i]) + 10.0 + (double)var.epsilonRefractory[i]);
								
								var.refractory[i] = newState;
								var.nextState[i] = var.refractory[i];
							}	   
						}
						else if (var.currentState[i] > 0)
						{
							var.nextState[i] = var.currentState[i] - 1;
						}
					}
					
					
					
					// reset the number of No. of active neighbours
					var.activeNeighbours[i] = 0;
					
			}
				
			/**********************/
			/* Update and storage */
			/**********************/
			
			// in the format: 
			// first domanSize^2 column: flattened currentState
			// second to the last column: inAF
			// the last column: eg
			
			temp_var = 0;
			temp_var2 = 0;
		
			
			for (int i=0;i<par.numberParticles;i++)
			{
				
				//cout << var.currentState[i] << " ";
				var.currentState[i] = var.nextState[i];
		
				if (var.currentState[i] >= var.refractory[i]-4)
				{
					// this cell is active, it contributes to its neighbours
					for (int j = 1;j <= par.nbhdMap[i][0]; j++)
					{
						int target = par.nbhdMap[i][j];
						
						var.activeNeighbours[target]++;						
					}
				}
				
			
				if (var.currentState[i] >= 0)
				{
					temp_var++;
					

					if (var.currentState[i] >= 1)
					{	
						
						temp_var2++;
					
					}
				}			
					
				if ((t%par.detailedRecordDt==0)&&(t>=par.detailedRecordStart)&&(t<=par.detailedRecordEnd))
					output2 << var.currentState[i] << "\t";
				
			}
		
			
			if ((t%par.detailedRecordDt==0)&&(t>=par.detailedRecordStart)&&(t<=par.detailedRecordEnd))
				output2 << endl;
			
			if ((t%par.coarseRecordDt==0)&&(t>=par.coarseRecordStart)&&(t<=par.coarseRecordEnd))
				output << t << "\t" << (double)temp_var2/(double)temp_var << endl;		
			

			
		}

		var.pvBurst = false;
		
		for (int t = (int)(maxPaceBeatNum*par.sinusPeriod);t< (int)((maxPaceBeatNum+maxWaitBeatNum)*par.sinusPeriod);t++)
		{
			
			if (t%par.sinusPeriod==0)
				cout << "Executing beat =" << t/par.sinusPeriod << "PV bursts: OFF" << endl;			
			
			if ((t%par.detailedRecordDt==0)&&(t>=par.detailedRecordStart)&&(t<=par.detailedRecordEnd))
				output2 << t << "\t";
			
			for (int i=0;i<par.numberParticles;i++)
			{
				tempRegistry[t%registryLength][i] = var.currentState[i];
			}
			
			if (t==var.nextBeat)
			{
				var.nextBeat += par.sinusPeriod;
				
				for (int i=0;i<par.snNumber;i++)
				{
					var.currentState[par.sinusNode[i]] = var.refractory[par.sinusNode[i]];
					
					if (var.actTime[par.sinusNode[i]] > 0)
					{
						var.cycleLength[par.sinusNode[i]] = t - var.actTime[par.sinusNode[i]];
						var.diastolInt[par.sinusNode[i]] = var.cycleLength[par.sinusNode[i]]-var.refractory[par.sinusNode[i]];
					}
					var.actTime[par.sinusNode[i]] = t;
				}
			
			}
			

			if ((var.pvBurst) && (rndexclusive() < par.pvBurstActivationProb))
			{
				// excite one set of cells (with uniform probability)
				
				int index = (int) (rndexclusive() * par.pvActivationCandidateNumber);
				
				for (int target = 1;target <= par.pvActivationMap[index][0];target++)
				{		
					int index2 = par.pvActivationMap[index][target];
					
					if (var.cellType[index2]==0)
					{
						if (var.actTime[index2] > 0)
						{
							var.cycleLength[index2] = t - var.actTime[index2];
							var.diastolInt[index2] = var.cycleLength[index2]-var.refractory[index2];
						}
								
						var.actTime[index2] = t;
						
						
						newState = max(floor((double)var.refractory1[index2]/2) + 5, floor(par.restitutionA*(1-par.restitutionB*exp( -(double)var.diastolInt[index2]/par.restitutionTau)) ) );
			
						var.refractory[index2] = newState;
						var.currentState[index2] = var.refractory[index2];
						
					}					
				}
				
			}
			
			/****************************************************************************/		
			/* determine next state of each "cell" based on number of active cells in   */
			/* % neighbourhood															*/
			/****************************************************************************/
			
			for (int i=0;i<par.numberParticles;i++)
			{
				
				
				
					if (var.cellType[i] == 0)
					{
						if (var.currentState[i] == 0)
						{				    
							if (var.activeNeighbours[i] >= (double)par.neighbourThreshold)
							{
								if (var.actTime[i] > 0)
								{
									var.cycleLength[i] = t - var.actTime[i];
									var.diastolInt[i] = var.cycleLength[i]-var.refractory[i];
								}
								
								var.actTime[i] = t;
								
								newState = max(floor((double)var.refractory1[i]/2) + 5, floor(par.restitutionA*(1-par.restitutionB*exp( -(double)var.diastolInt[i]/par.restitutionTau)) ) );
								//newState = max(floor((double)var.refractory[i]/2) + 5, floor(sqrt(par.restitutionA*(double)var.cycleLength[i]) +  par.refractoryPeriod- sqrt(par.restitutionA*(double)par.sinusPeriod) + (double)var.epsilonRefractory[i]));
								//newState = floor(sqrt(8*(double)var.cycleLength[i]) + 10.0 + (double)var.epsilonRefractory[i]);
								
								var.refractory[i] = newState;
								var.nextState[i] = var.refractory[i];
							}	   
						}
						else if (var.currentState[i] > 0)
						{
							var.nextState[i] = var.currentState[i] - 1;
						}
					}
								
					var.activeNeighbours[i] = 0;
					
			}
				
			/**********************/
			/* Update and storage */
			/**********************/
			
			
			temp_var = 0;
			temp_var2 = 0;
		
			
			for (int i=0;i<par.numberParticles;i++)
			{
				
				//cout << var.currentState[i] << " ";
				var.currentState[i] = var.nextState[i];
		
				if (var.currentState[i] >= var.refractory[i]-4)
				{
					// this cell is active, it contributes to its neighbours
					for (int j = 1;j <= par.nbhdMap[i][0]; j++)
					{
						int target = par.nbhdMap[i][j];
						
						var.activeNeighbours[target]++;						
					}
				}
				
			
				if (var.currentState[i] >= 0)
				{
					temp_var++;
					

					if (var.currentState[i] >= 1)
					{	
						
						temp_var2++;
					
					}
				}			
					

				if ((t%par.detailedRecordDt==0)&&(t>=par.detailedRecordStart)&&(t<=par.detailedRecordEnd))
					output2 << var.currentState[i] << "\t";
				
			}
		
			
			if ((t%par.detailedRecordDt==0)&&(t>=par.detailedRecordStart)&&(t<=par.detailedRecordEnd))
				output2 << endl;
			
			if ((t%par.coarseRecordDt==0)&&(t>=par.coarseRecordStart)&&(t<=par.coarseRecordEnd))
				output << t << "\t" << (double)temp_var2/(double)temp_var << endl;		
			
			
			if (temp_var2 == 0)
			{
				last_0 = t;
				inAF = 0;
			}
			
			if ((t - last_0 > threshold)&&(inAF==0))
			{	

				inAF = 1;
				
			}
			
			
			
		}

		
		cout << "Random seed=" << RS1 << ", inAF=" << inAF << endl;
		
		if (inAF)
			ensIndex = ensN;
		
		
	}
	
	output.close();	
	output2.close();
	
	
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "success=" << Success << ", execution time=" << elapsed_secs << " secs" << endl;
	
	return 0;
	

}



double rndexclusive()
{
	double i=0;
	
	while ((i==0)||(i==1))
	{
		i = (double)rand()/RAND_MAX;
		// or use your own favorite/efficient random number generators
	}
	
	return i;
}  

twin normal()
{
	twin out;
	
	double r1 = rndexclusive();
	double r2 = rndexclusive();
	
	//cout << M_PI << endl;
	
	out.n1 = sqrt(-2*log(r1))* cos(2*M_PI*r2);
	out.n2 = sqrt(-2*log(r1))* sin(2*M_PI*r2);
	
	return out;
}
