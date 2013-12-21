// Johannes Ulén, 2013
// MATLAB wrapper for Roof duality (QPBO) by Vladimir Kolmogorov
// http://pub.ist.ac.at/~vnk/software.html
#include "mex.h"
#include "QPBO-v1.3.src/QPBO.h"
#include "utils/mexutils.h"
#include "utils/cppmatrix.h"

// Catch errors
static void erfunc(char *err) {
	mexErrMsgTxt(err);
}

void mexFunction(int nlhs,      						/* number of expected outputs */
                 mxArray        *plhs[],    /* mxArray output pointer array */
                 int            nrhs,       /* number of inputs */
                 const mxArray  *prhs[]     /* mxArray input pointer array */)
{
	// Parse input
	ASSERT(nrhs == 7 || nrhs == 8);
	ASSERT(nlhs == 4);

	int curarg = 0;
	matrix<double> U0(prhs[curarg++]);
	matrix<double> U1(prhs[curarg++]);

	matrix<double> E00(prhs[curarg++]);
	matrix<double> E01(prhs[curarg++]);
	matrix<double> E10(prhs[curarg++]);
	matrix<double> E11(prhs[curarg++]);
	matrix<unsigned int> connectivity(prhs[curarg++]);

  	MexParams params(nrhs-curarg, prhs+curarg); //Structure to hold and parse additional parameters
  	bool improve = params.get<bool>("improve", false); // RD/QPBO Improve

	ASSERT(U0.M == U1.M);
	ASSERT(U1.N == U1.N);
	ASSERT(U0.N == 1);
	ASSERT(U1.N == 1);

	ASSERT(E00.N == E01.N);
	ASSERT(E01.N == E10.N);
	ASSERT(E10.N == E11.N);
	ASSERT(E11.N == connectivity.N);

	ASSERT(E00.M == 1);
	ASSERT(E01.M == 1);
	ASSERT(E10.M == 1);
	ASSERT(connectivity.M == 2);

	int unary_terms = U0.M;
	int pairwise_terms = E00.N;

	// Create problem instance
	QPBO<double> graph(unary_terms, pairwise_terms, erfunc);
	graph.AddNode(unary_terms);

	for (int p = 0; p < pairwise_terms; p++)
		graph.AddPairwiseTerm(connectivity(0,p), connectivity(1,p), E00(p), E01(p), E10(p), E11(p));

	for (int u = 0; u <  unary_terms; u++)
		graph.AddUnaryTerm(u, U0(u), U1(u));

	// Merge edges
	graph.MergeParallelEdges();

	// Solve for optimimum and label weak persitencies
	graph.Solve();
	graph.ComputeWeakPersistencies();

	// Defining outout
	matrix<double> labelling(unary_terms);
	matrix<double> energy(1);
	matrix<double> lower_bound(1);
	matrix<double> num_unlabelled(1);
	
	plhs[0] = labelling;
	plhs[1] = energy;
	plhs[2] = lower_bound;
	plhs[3] = num_unlabelled;

	// Count unlabelled before improve
	num_unlabelled(0) = 0;
	for (int u = 0; u < unary_terms; u++)
	{
		if ((graph.GetLabel(u) < 0))
			num_unlabelled(0)++;
	}

	// if any labels are unlabelled run improve
	if (improve && (num_unlabelled(0) > 0))
		graph.Improve();

	// Set t zero for time being
	energy(0) = 0;
	for (int u = 0; u < unary_terms; u++)
		labelling(u) = graph.GetLabel(u);

	energy(0) = graph.ComputeTwiceEnergy()/2;
	lower_bound(0) = graph.ComputeTwiceLowerBound()/2;
}