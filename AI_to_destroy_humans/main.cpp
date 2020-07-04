#include "aitdh.h"

// Troian Borys KV-62, FAM
// Variant 1.1.2
//	Genetic algorythm
//	Real-coded
//	Tournament+Elitism selection
//	Linear crossover
//	Random mutation
//	Children replace parents


double HSphereFitness(gene &gn)		//[-100; 100]
{
	double res = 0;

	for (int i = 0; i < DIM; i++)
		res += gn.alleles[i] * gn.alleles[i];
	if (res == 0)
		return DBL_MAX;
	else
		return 1 / res;
}

double AckleyFitness(gene &gn)		//[-32.768; 32.768]
{
	double a = 20, b = 0.2, c = 2 * M_PI;
	double p1=0, p2=0, A1, A2;
	for (int i = 0; i < DIM; i++)
	{
		p1 += gn.alleles[i] * gn.alleles[i];
		p2 += cos(2 * c*gn.alleles[i]);
	}
	A1 = exp((-b)*sqrt(p1 / DIM));
	A2 = exp(p2 / DIM);
	return 1/(-a * A1 - A2  + a - M_E);
}

double GriewankFitness(gene &gn)	//[-600; 600]
{
	double p1 = 0, p2 = 1;
	for (int i = 0; i < DIM; i++)
	{
		p1 += gn.alleles[i] * gn.alleles[i];
		p2 *= cos(gn.alleles[i] / sqrt(i + 1));
	}
	return 1/(p1/4000 - p2 + 1);
}

double RastriginFitness(gene &gn)	//[-5.12; 5.12]
{
	double p1 = 0;
	for (int i = 0; i < DIM; i++)
	{
		p1 += (gn.alleles[i] * gn.alleles[i]-10*cos(2*M_PI*gn.alleles[i]));
	}
	return 1 / (10*DIM + p1);
}

double RosenbrockFitness(gene &gn)		//[-5;10]
{
	double p1 = 0, x = 0;
	for (int i = 0; i < DIM-1; i++)
	{
		x = gn.alleles[i + 1] - gn.alleles[i] * gn.alleles[i];
		p1 += 100 * x*x + (gn.alleles[i] - 1)*(gn.alleles[i] - 1);
	}
	return 1 / p1;
}

int main()
{
	int ans;
	GenAISolver  A(*HSphereFitness, -100, 100);
	//GenAISolver  A(*AckleyFitness, -32.768, 32.768);
	//GenAISolver  A(*GriewankFitness, -600, 600, 20);
	//GenAISolver  A(*RastriginFitness, -5.12, 5.12, 0.05);
	//GenAISolver  A(*RosenbrockFitness, -5, 10, 0.1);
	A.Solve();
	gene aa = A.GetBest();

	//printf("%f\n", 1/A.GetBest().fitness);
	getchar();
	return 0;
}