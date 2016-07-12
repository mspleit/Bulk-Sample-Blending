// LM_Blending.cpp : Application that uses CPLEX to optimize a bulk sample's proportions of sub-samples to ensure grades are within 
desired target bounds.
//
using namespace std;
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "time.h"

#include <ilcplex/ilocplex.h>

#include "StdAfx.h"

class Sim {
public:
	float feh, dtwr, fec, sic;
};

class Sample {

public:
	char source[25], label[25], litho[25];
	int bin, sim, weight;
	int numSims;
	vector<float> feh, dtwr, fec, sic;
};

int main(int argc, char **argv)
{
	//read raw data from file into array
		//
	
	vector<Sample> samples;
	const char* filename = "blend5b.csv";
	//const char* filename = "noPcoreLC1.csv";
	//string filename = "blend.csv";
	int numSamples=0;
	int s=0;
	/* initialize random seed: */
	srand (time(NULL));

	FILE *infile;  
	infile = fopen(filename,"r");
	//ifstream infile;
	//infile.open(filename.c_str());
	if(!infile){   // file couldn't be opened
      cout << "Error: file could not be opened" << endl;
      exit(1);
	}

	cout << "Reading sample data file..." << endl;
	int i=0;
	int ret;
	
	while ( !feof(infile)){  // keep reading until end-of-file
		Sample mySample;
		float myfeh, mydtwr,myfec,mysic;
		ret = fscanf(infile,"%25[^,],", mySample.source);
		ret = fscanf(infile,"%25[^,],", mySample.label);
		ret = fscanf(infile,"%25[^,],", mySample.litho);
		ret = fscanf(infile,"%d,", &mySample.bin);
		ret = fscanf(infile,"%d,", &mySample.sim);
		ret = fscanf(infile,"%d,", &mySample.weight);
	
		if (numSamples>0) { //check we're on at least sample #2 so we know there's a numSamples-1
			//if the previous sample had the same litho and bin, do NOT push_back the sample, but add another simulation (the grades)
			cout << samples[numSamples-1].litho << endl;
			cout << samples[numSamples-1].bin << endl;
			if (//samples[numSamples].litho.compare(samples[numSamples-1].litho)) {
				strcmp(mySample.litho,samples[numSamples-1].litho)==0 && mySample.bin==samples[numSamples-1].bin) {
				samples[numSamples-1].numSims++;
				ret = fscanf(infile, "%f,%f,%f,%f", &myfeh, &mydtwr, &myfec, &mysic);
				samples[numSamples-1].feh.push_back(myfeh);
				samples[numSamples-1].dtwr.push_back(mydtwr);
				samples[numSamples-1].fec.push_back(myfec);
				samples[numSamples-1].sic.push_back(mysic);
			}
			else {//otherwise, we DO push_back the sample as the first one for that bin, then add the grades
				samples.push_back(mySample);
				samples[numSamples].numSims=1;
				ret = fscanf(infile, "%f, %f, %f, %f", &myfeh, &mydtwr, &myfec, &mysic);
				samples[numSamples].feh.push_back(myfeh);
				samples[numSamples].dtwr.push_back(mydtwr);
				samples[numSamples].fec.push_back(myfec);
				samples[numSamples].sic.push_back(mysic);
				numSamples++;
			}
		}
		else {
			// this is the first sample, so add it to the vector and then add the grades as sim1
			samples.push_back(mySample);
			samples[numSamples].numSims=1;
			ret = fscanf(infile, "%f, %f, %f, %f", &myfeh, &mydtwr, &myfec, &mysic);
			samples[numSamples].feh.push_back(myfeh);
			samples[numSamples].dtwr.push_back(mydtwr);
			samples[numSamples].fec.push_back(myfec);
			samples[numSamples].sic.push_back(mysic);
			numSamples++;
		}
	}
	   
	fclose(infile);
	cout << "Last line of input:" << endl;
	cout << "Input file read." << endl;

	// Generate 100 simulations
	Sim mySim;
	vector<Sim> sims[100];
	int random_sim;
	int j;
	for (i=0;i<100;i++){
		for (j=0;j<numSamples;j++){
			/* generate random number between 0 and number of simulations-1 for current sample: */
			random_sim = rand() % samples[j].numSims;

			mySim.feh=samples[j].feh[random_sim];
			mySim.dtwr=samples[j].dtwr[random_sim];
			mySim.fec=samples[j].fec[random_sim];
			mySim.sic=samples[j].sic[random_sim];
			sims[i].push_back(mySim);
		}	
	}
	cout <<"Simulations generated." <<endl;
	
	ILOSTLBEGIN
	
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);
	float sic_max = 0.022;
	float sic_min = 0.022;
	
	float LC_max=0.179;
	float LC_min=0.133;
	float JUIF_max=0.207;
	float JUIF_min=0.153;
	float GC_max=0.023;
	float GC_min=0.017;
	float URC_max=0.097;
	float URC_min=0.071;
	float PGC_max=0.164;
	float PGC_min=0.122;
	float LRC_max=0.054;
	float LRC_min=0.040;
	float LRGC_max=0.424;
	float LRGC_min=0.314;

	IloExpr objective_value(env);
	IloNumVarArray q(env, numSamples,0,IloInfinity,ILOFLOAT);
	IloExprArray blend_sic_upper(env,100);
	IloExprArray blend_sic_lower(env,100);
	IloExprArray blend_sic(env,100);
	IloExpr blend_weight(env);

	IloNumVarArray dev_sic_upper, dev_sic_lower;
	
	IloExpr LC_weight(env);
	IloExpr JUIF_weight(env);
	IloExpr GC_weight(env);
	IloExpr PGC_weight(env);
	IloExpr URC_weight(env);
	IloExpr LRC_weight(env);
	IloExpr LRGC_weight(env);

	IloNumVar dev_LC_upper(env,0,IloInfinity,ILOFLOAT), dev_LC_lower(env,0,IloInfinity,ILOFLOAT);
	IloNumVar dev_JUIF_upper(env,0,IloInfinity,ILOFLOAT), dev_JUIF_lower(env,0,IloInfinity,ILOFLOAT);
	IloNumVar dev_GC_upper(env,0,IloInfinity,ILOFLOAT), dev_GC_lower(env,0,IloInfinity,ILOFLOAT);
	IloNumVar dev_URC_upper(env,0,IloInfinity,ILOFLOAT), dev_URC_lower(env,0,IloInfinity,ILOFLOAT);
	IloNumVar dev_PGC_upper(env,0,IloInfinity,ILOFLOAT), dev_PGC_lower(env,0,IloInfinity,ILOFLOAT);
	IloNumVar dev_LRC_upper(env,0,IloInfinity,ILOFLOAT), dev_LRC_lower(env,0,IloInfinity,ILOFLOAT);
	IloNumVar dev_LRGC_upper(env,0,IloInfinity,ILOFLOAT), dev_LRGC_lower(env,0,IloInfinity,ILOFLOAT);

	dev_sic_upper = IloNumVarArray(env, 100, 0, IloInfinity, ILOFLOAT);
	dev_sic_lower = IloNumVarArray(env, 100, 0, IloInfinity, ILOFLOAT);
	

	for (i=0;i<numSamples;i++) {
		model.add(q[i]<=samples[i].weight); //selected sample weight must be less than what's available
		blend_weight += q[i];

		if (strcmp(samples[i].litho,"LC")==0) {
				LC_weight += q[i];}

			if (strcmp(samples[i].litho,"JUIF")==0) {
				JUIF_weight += q[i];}
			
			if (strcmp(samples[i].litho,"GC")==0) {
				GC_weight += q[i];}

			if (strcmp(samples[i].litho,"URC")==0) { 
				URC_weight += q[i];}
			
			if (strcmp(samples[i].litho,"PGC")==0) {
				PGC_weight += q[i];}
			
			if (strcmp(samples[i].litho,"LRC")==0) {
				LRC_weight += q[i];}
			
			if (strcmp(samples[i].litho,"LRGC")==0) {
				LRGC_weight += q[i];}
			
			if (strcmp(samples[i].litho,"2014")==0) {
				LC_weight += q[i]*0.1427;
				JUIF_weight += q[i]*0.1287;
				GC_weight += q[i]*0.02722;
				URC_weight += q[i]*0.08310;
				PGC_weight += q[i]*0.2038;
				LRC_weight += q[i]*0.07063;
				LRGC_weight += q[i]*0.3440;
			}
	}

	for (int s=0;s<100;s++){
		blend_sic_upper[s] = IloExpr(env);
		blend_sic_lower[s] = IloExpr(env);
		blend_sic[s] = IloExpr(env);

		for(i=0;i<numSamples;i++){
			blend_sic_upper[s] += q[i]*(sims[s][i].sic - sic_max);
			blend_sic_lower[s] += q[i]*(sic_min - sims[s][i].sic);
			blend_sic[s] += q[i]*sims[s][i].sic;
		}
		//model.add(blend_sic_upper[s] - dev_sic_upper[s] <= 0); //deviation for each simulation from the blend target (over)
		//model.add(blend_sic_lower[s] - dev_sic_lower[s] <= 0); //deviation for each simulation from the blend target (under)

		model.add(blend_sic[s] - dev_sic_upper[s] <= sic_max*blend_weight);
		model.add(blend_sic[s] + dev_sic_lower[s] >= sic_min*blend_weight);
		
		model.add(LC_weight - dev_LC_upper <= LC_max*blend_weight);
		model.add(LC_weight + dev_LC_lower >= LC_min*blend_weight);
		model.add(JUIF_weight - dev_JUIF_upper <= JUIF_max*blend_weight);
		model.add(JUIF_weight + dev_JUIF_lower >= JUIF_min*blend_weight);
		model.add(GC_weight - dev_GC_upper <= GC_max*blend_weight);
		model.add(GC_weight + dev_GC_lower >= GC_min*blend_weight);
		model.add(URC_weight - dev_URC_upper <= URC_max*blend_weight);
		model.add(URC_weight + dev_URC_lower >= URC_min*blend_weight);
		model.add(PGC_weight - dev_PGC_upper <= PGC_max*blend_weight);
		model.add(PGC_weight + dev_PGC_lower >= PGC_min*blend_weight);
		model.add(LRC_weight - dev_LRC_upper <= LRC_max*blend_weight);
		model.add(LRC_weight + dev_LRC_lower >= LRC_min*blend_weight);
		model.add(LRGC_weight - dev_LRGC_upper <= LRGC_max*blend_weight);
		model.add(LRGC_weight + dev_LRGC_lower >= LRGC_min*blend_weight);
		/*
		model.add(blend_LC_upper + dev_LC_upper <= 0);
		model.add(blend_LC_lower + dev_LC_lower <= 0);
		model.add(blend_JUIF_upper + dev_JUIF_upper <= 0);
		model.add(blend_JUIF_lower + dev_JUIF_lower <= 0);
		model.add(blend_GC_upper + dev_GC_upper <= 0);
		model.add(blend_GC_lower + dev_GC_lower <= 0);
		model.add(blend_URC_upper + dev_URC_upper <= 0);
		model.add(blend_URC_lower + dev_URC_lower <= 0);
		model.add(blend_PGC_upper + dev_PGC_upper <= 0);
		model.add(blend_PGC_lower + dev_PGC_lower <= 0);
		model.add(blend_LRC_upper + dev_LRC_upper <= 0);
		model.add(blend_LRC_lower + dev_LRC_lower <= 0);
		model.add(blend_LRGC_upper + dev_LRGC_upper <= 0);
		model.add(blend_LRGC_lower + dev_LRGC_lower <= 0);*/
	}
	
	/*objective = max( sum over samples(qi) - penalties for SiC deviations - penalties for litho_blend deviations) 
	penalties for SiC:
	sum over simulations (cost_factor * SiC_dev(s))     i.e. the sum of the blended SiC deviation for each simulation */
	
	for (i=0;i<numSamples;i++){
		objective_value += q[i];
	}
	
	for (s=0;s<100;s++){
		objective_value -= 3.1*(dev_sic_upper[s] + dev_sic_lower[s]);
		/* Shouldn't these penalties be weighted by q[i]? 
		   Ok. They are. dev_sic_upper and _lower are weighted.*/
		/* Doesn't this just penalize deviation of the Expected Value from the target, but not control probability range? 
		   No, for each simulation s, it is penalizing the deviation of the Expected Value, so aggregate probability is being considered.*/
	}
	
	objective_value -= (0*(dev_LC_upper + dev_LC_lower) + 0*(dev_JUIF_upper + dev_JUIF_lower + dev_GC_upper + dev_GC_lower +
		dev_URC_upper + dev_URC_lower + dev_PGC_upper + dev_PGC_lower + dev_LRC_upper + dev_LRC_lower + dev_LRGC_upper + dev_LRGC_lower));
	
	model.add(IloMaximize(env, objective_value)); 

	try{
		cplex.extract(model);
	}
	catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }
   
   // Optimize the problem and obtain solution.
    if ( !cplex.solve() ) {
		env.error() << "Failed to optimize LP" << endl;
//        throw(-1);
    }
	else {
		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value  = " << cplex.getObjValue() << endl;
		env.out() << "Sample weights for solution: " << endl;
		ofstream outfile;
		outfile.open("output.csv");
		if(!outfile)
		{   // file couldn't be opened
			cerr << "Error: Output file could not be opened" << endl;
			exit(1);
		}	
			
		float solution_weight[18];
		for(i=0;i<numSamples;i++){
			solution_weight[i]=cplex.getValue(q[i]);
			env.out() << solution_weight[i] << endl;
			outfile << samples[i].litho << "," << samples[i].weight << "," << solution_weight[i];
			for (s=0;s<100;s++){
				outfile << "," << sims[s][i].sic;
			}
			outfile << endl;
		}
		float solution_sic_upper, solution_sic_lower;
		env.out() << "SiC deviations (upper, lower for each simulation)" << endl;
		int s;
		for (s=0;s<100;s++){
			env.out() << cplex.getValue(dev_sic_upper[s]) << " " << cplex.getValue(dev_sic_lower[s]) << endl;
			//getsolution_sic_upper=cplex.getValue(dev_sic_upper[s]);
			//solution_sic_lower=cplex.getValue(dev_sic_lower[s]);
			//env.out() << solution_sic_upper << " " << solution_sic_lower << endl;
	}
	
	cout << "Optimization complete." << endl;
	outfile.close();
	cout << "Output file created." << endl;
	}
	std::cin.get(); //wait for keypress
	return 0;
}
