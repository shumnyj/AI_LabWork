#pragma once
//#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <time.h>
#include <cstdlib>
#include <chrono>
#include <stdio.h>
#include <math.h>
#include <limits>
#include <random>


#define	MAXPOP 200
#define MAXITER 50000
#define DIM 10
#define ELITES 4	//even number



struct gene {
	double alleles[DIM];
	double fitness;
	double likelihood = 0;

	// Test for equality.
	bool operator==(gene gn) {
		for (int i = 0; i < DIM; i++) {
			if (gn.alleles[i] != alleles[i]) return false;
		}
		return true;
	}
};


class GenAISolver 
{
	public:
		GenAISolver(double (*fitfunc)(gene &gn), double L, double R, double mlim = 0.4);		// Constructor 

		double Solve();								// Solve the equation.
		// Returns a given gene.
		gene GetGene(int i) { return population[i]; }
		gene GetBest();
		

	protected:
		std::default_random_engine localgen;
		std::uniform_real_distribution<double> LimitedRand;
		std::uniform_real_distribution<double> MutRandRGA;

		double peak;
		gene population[MAXPOP];					// Population.
		double rlim, llim;
		double mut_lim;
		int best;

		double(*FitnessFunc)(gene &gn);
		//gene GenAISolver::Breed(int p1, int p2)
		void Breed(int p1, int p2, gene &c1, gene &c2);
		int TournamentRound();
		int CreateFitnesses();
		void CreateNewPopulation();

		//deprecated
		int FetchIndex(float val);	
		void GenerateLikelihoods();		// Generate likelihoods to choose gene with FetchIndex()
		float MultInv();
};
