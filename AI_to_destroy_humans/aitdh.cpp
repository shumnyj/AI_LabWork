#include "aitdh.h"

GenAISolver::GenAISolver(double(*fitfunc)(gene &gn), double L, double R, double mlim )
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine dreA(seed);
	std::uniform_real_distribution<double> urdA(L,R);
	std::uniform_real_distribution<double> urdB(-mlim, mlim);
	localgen = dreA;
	LimitedRand = urdA;
	MutRandRGA = urdB;
	best = 0;

	FitnessFunc = fitfunc;
	rlim = R;
	llim = L;
}

double GenAISolver::Solve() 
{
	int fitness = -1;
	int iterations = 0;
	// Generate initial population.
	srand((unsigned)time(NULL));

	for (int i = 0; i < MAXPOP; i++)				// Initial fill
	{
		for (int j = 0; j < DIM; j++)
			population[i].alleles[j] = LimitedRand(localgen);
	}
	CreateFitnesses();
	/*if (CreateFitnesses())
		return best;*/

	while (iterations < MAXITER)	// Repeat until solution found, or iteraton limit reached.
	{		
		//GenerateLikelihoods();		// Create the likelihoods.
		CreateNewPopulation();
		//CreateFitnesses();
		/*if (CreateFitnesses())		//if extr found - negligible chance
			return best;*/
		iterations++;
	}
	for (int i = 0; i < MAXPOP; i++)				// Finding best result
	{
		if (population[best].fitness < population[i].fitness)
			best = i;
	}

	return best;
}

gene GenAISolver::GetBest()
{ 
	printf("Alleles: ");
	for (int i = 0; i < DIM; i++)
		printf("%8f ", population[best].alleles[i]);
	printf("\nValue = %.8f\n", 1/population[best].fitness);
	return population[best]; 
}

int GenAISolver::CreateFitnesses()	//DOES NOT WORK FOR NEGATIVES // might define some value that 100% below minimum and compare to it in fitness
{
	//float avgfit = 0;
	int fitness = 0;
	for (int i = 0; i < MAXPOP; i++) 
	{
		population[i].fitness = FitnessFunc(population[i]);
		//avgfit += population[i].fitness;
		if (population[i].fitness == DBL_MAX)		//if extr found return
		{
			best = i;
			return 1;
		}
	}
	return 0;
}

float GenAISolver::MultInv() 
{
	float sum = 0;
	int i;
	for (i = 0; i < MAXPOP; i++)
		sum += (float)population[i].fitness;			//watch out for overflow

	return sum;
}

void GenAISolver::GenerateLikelihoods() 
{
	float multinv = MultInv();

	float last = 0;
	for (int i = 0; i < MAXPOP; i++)
		population[i].likelihood = last = last + (population[i].fitness / multinv * 100);
}

int GenAISolver::FetchIndex(float val) 
{
	float last = 0;
	for (int i = 0; i < MAXPOP; i++) 
	{
		if (last <= val && val <= population[i].likelihood) 
			return i;
		else last = population[i].likelihood;
	}

	return 0;
}

int GenAISolver::TournamentRound()		//binary tournament
{
	int a, b, i=0;
	do {
		a = rand() % MAXPOP;
		b = rand() % MAXPOP;
		i++;
	} while ((a == b || population[a] == population[b]) && i < 25);
	if (population[a].fitness >= population[b].fitness)
		return a;
	else
		return b;
}

/*gene GenAISolver::Breed(int p1, int p2) 
{
	int crossover = rand() % DIM;					// Create the crossover point (not first).
	//int first = rand() % 100;						// Which parent comes first?

	gene child = population[p1];					// Child is all first parent initially.

	int initial = 0, final = DIM-1;					// The crossover boundaries.
	if (rand() % 100 < 50) initial = crossover;		// If first parent first. start from crossover.
	else final = crossover + 1;						// Else end at crossover.

	for (int i = initial; i < final; i++)						// Crossover
	{				
		child.alleles[i] = population[p2].alleles[i];
		if (rand() % 101 < 5) child.alleles[i] = LimitedRand(localgen);	//mutation
	}

	return child;									// Return the child
}*/

void GenAISolver::Breed(int p1, int p2, gene &c1, gene &c2)
{
	gene h1, h2, h3;		// Child is all first parent initially.

	for (int i = 0; i < DIM; i++)
	{
		h1.alleles[i] = (population[p1].alleles[i] + population[p2].alleles[i]) / 2;
		h2.alleles[i] = (3*population[p1].alleles[i] - population[p2].alleles[i]) / 2;
		h3.alleles[i] = (-population[p1].alleles[i] + 3*population[p2].alleles[i]) / 2;
		if (rand() % 101 < 5) h1.alleles[i] += MutRandRGA(localgen);
		if (rand() % 101 < 5) h2.alleles[i] += MutRandRGA(localgen);
		if (rand() % 101 < 5) h3.alleles[i] += MutRandRGA(localgen);
	}
	h1.fitness = FitnessFunc(h1);
	h2.fitness = FitnessFunc(h2);
	h3.fitness = FitnessFunc(h3);

	c1 = h1;							//max h1 h2 
	c2 = h2;
	if (h3.fitness > c1.fitness)		//max h2 h3
		c1 = h3;
	else if (h3.fitness > c2.fitness)	//max h1 h3
		c2 = h3;
}

void GenAISolver::CreateNewPopulation() 
{
	gene temppop[MAXPOP];
	int E[ELITES];
	int i = 0, j = 0, k = 0, p = 0;
	bool f;
	for (i = 0; i < ELITES; i++)		//fill all elites
	{
		E[i] = 0;
		for (j = 0; j < MAXPOP; j++)	//scan population
		{
			if (population[j].fitness < population[E[i]].fitness)	//if found better/equal than current
			{
				f = true;
				for (k = 0; k < i; k++)
				{
					if (E[k] == j)			//check if not unique
						f = false;
				}
				if (f)
					E[i] = j;				//else replace current 
			}
		}
		temppop[MAXPOP - i - 1] = population[E[i]];	//saving elites into new population, any order works
	}
	for (i = 0; i < MAXPOP-ELITES; i+=2)
	{
		int parent1 = 0, parent2 = 0, iterations = 0;
		/*while (parent1 == parent2 || population[parent1] == population[parent2]) 
		{
			parent1 = FetchIndex((float)(rand() % 1001)/10);
			parent2 = FetchIndex((float)(rand() % 1001)/10);

			if (++iterations > 25) break;
		}*/
		parent1 = TournamentRound();
		parent2 = TournamentRound();
		
		Breed(parent1, parent2, temppop[i], temppop[i + 1]);		// Create a child.
	}

	for (int i = 0; i < MAXPOP - ELITES; i++)		//put results in new population
		population[i] = temppop[i];
}
