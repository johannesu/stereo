//[L energy lower_bound] = vgg_trw_bp(UE, PI, PE, options);

// $Id: vgg_trw_bp.cxx,v 1.3 2009/08/31 22:05:09 ojw Exp $

typedef double REAL; 

#include <mex.h>
#include "MRFEnergy.h"

#include <string>
#include <stdexcept>
using std::string;
using std::runtime_error;
#include <sstream>
#define ASSERT(cond) if (!(cond)) { std::stringstream sout; \
									sout << "Error (file " << __FILE__ << ", line " << __LINE__ << "): " << #cond; \
                                    throw runtime_error(sout.str()); }

void flush_output()
{
	mexEvalString("drawnow");
}

#define debug() { mexPrintf("Reached line %d in file %s. \n", __LINE__, __FILE__);	flush_output(); }


// Define types
#ifdef _MSC_VER
typedef unsigned __int16 uint16_t;
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

static void erfunc(char *err) {mexErrMsgTxt(err);}
template<class TYPE> static inline void wrapper_func(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], int options[], int *nlabels, REAL lambda);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of inputs
	if (nlhs < 1 || nlhs > 3)
		mexErrMsgTxt("Unexpected number of outputs.");

	for (int i = 0; i < nrhs; i++) {
		if (mxIsComplex(prhs[i]))
			mexErrMsgTxt("Inputs must be real.");
	}

	// Read threshold
	double *data  = (double *)mxGetData(prhs[3]);
	REAL lambda = data[0];

	// Check unary terms
	if (!mxIsCell(prhs[0]))
		mexErrMsgTxt("UE must be a cell array.");

	int n_nodes = mxGetNumberOfElements(prhs[0]);
	int *nlabels = new int[n_nodes];
	int max_labels = 0;
	int min_labels = 65537;
	for (int i = 0; i < n_nodes; i++) {
		mxArray *data_array = mxGetCell(prhs[0], i);

		if (!mxIsDouble(data_array))
			mexErrMsgTxt("UE cells must be doubles.");

		if (mxIsComplex(data_array))
			mexErrMsgTxt("UE cells must be real.");
		nlabels[i] = mxGetN(data_array);
		max_labels = nlabels[i] > max_labels ? nlabels[i] : max_labels;
		min_labels = nlabels[i] < min_labels ? nlabels[i] : min_labels;
	}
	if (max_labels > 65536)
		mexErrMsgTxt("A maximum of 65536 nodes per label are supported.");

	// Check options
	int options[] = {-1, max_labels, 1, 0, 30}; // ParamsPerEdge, max_labels, UseTRW, Type, max_iters

		if (!mxIsInt32(prhs[2]))
			mexErrMsgTxt("options should be int32s");
		const int32_t *params = (const int32_t *)mxGetData(prhs[2]);
		switch (mxGetNumberOfElements(prhs[2])) {
			default:
			case 3:
				options[4] = (int)params[2];
			case 2:
				options[3] = (int)params[1];
			case 1:
				options[2] = (int)params[0] != 0;
			case 0:
				break;
		}


	// Call wrapper 
	// (I keep the strucutre from OJW's code even though it's kind of useless when we only have one type).
	wrapper_func<TypeStereo>(nlhs, plhs, nrhs, prhs, options, nlabels, lambda);

	delete nlabels;
	return;
}

// Functions for Stereo type
static inline MRFEnergy<TypeStereo>::NodeId add_node(MRFEnergy<TypeStereo> *mrf, int local_modes, TypeStereo::REAL *graph_data)
{
	return mrf->AddNode(TypeStereo::LocalSize(local_modes), TypeStereo::NodeData(graph_data));
}

static inline void add_edge( 
	MRFEnergy<TypeStereo> *mrf, MRFEnergy<TypeStereo>::NodeId node1, MRFEnergy<TypeStereo>::NodeId node2, 
	TypeStereo::REAL alpha, TypeStereo::REAL lambda,
	TypeStereo::REAL * stacked_data
	)
{
	mrf->AddEdge(node1, node2, TypeStereo::EdgeData(lambda, alpha, stacked_data) );
}

template<class TYPE> static inline void wrapper_func(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], int options[], int *nlabels, REAL lambda)
{
	ASSERT(nrhs == 13);
	int n_nodes = mxGetNumberOfElements(prhs[0]);

	// edges
	int Pindices = mxGetM(prhs[1]);
	if (Pindices != 2 && Pindices != 3)
		mexErrMsgTxt("Unexpected dimensions for PI");
	int nP = mxGetN(prhs[1]);

	// Number of Pairwise terms	
	if (!mxIsUint32(prhs[1]))
		mexErrMsgTxt("PI should be uint32s");
	const uint32_t *PI = (const uint32_t *)mxGetData(prhs[1]);

	// options[1] = max_labels
	int max_labels = options[1];

	MRFEnergy<TYPE> *mrf = new MRFEnergy<TYPE>(typename TYPE::GlobalSize(max_labels), erfunc);
	typename MRFEnergy<TYPE>::NodeId *nodes =  new typename MRFEnergy<TYPE>::NodeId[n_nodes];
	typedef typename TYPE::REAL REAL; 

	// Check data
	ASSERT(nP == mxGetN(prhs[4]))
	for (int i = 5; i < 13; i++)
	{
		ASSERT(nP == mxGetN(prhs[i]));
		ASSERT(mxIsCell(prhs[i]));
	}

 	REAL *graph_data = (REAL *)mxCalloc(max_labels*max_labels, sizeof(REAL));

	// Add unary energies
	for (int i = 0; i < n_nodes; i++) {
		mxArray *data_array = mxGetCell(prhs[0], i);
		int pitch = mxGetM(data_array);
		const double *data = mxGetPr(data_array);
		for (int j = 0; j < nlabels[i]; j++)
			graph_data[j] = (REAL)data[j*pitch];
		nodes[i] = add_node(mrf, nlabels[i], graph_data);
	}

	// All data sent in as one long vector 
	REAL *stacked_data = (REAL *)mxCalloc(2*(5*max_labels), sizeof(REAL));
	REAL alpha;

	// Add pairwise energies
	for (int i = 0; i < nP; i++, PI += Pindices) 
	{
		mxArray *data_array;
		const REAL *data;

		// Node indices
		int n1 = PI[0] - 1;
		int n2 = PI[1] - 1;
		
		data = mxGetPr(prhs[4]);
		alpha = data[i];
	
		// Stacka all data
		int j = 0; 		
		int curarg = 5;

        // Points
		for (curarg = 5; curarg < 7; curarg ++)
		{
			data_array = mxGetCell(prhs[curarg],i);	
			data = (REAL *)mxGetData(data_array);		

			for (int k = 0; k < mxGetNumberOfElements(data_array); j++, k++)
				stacked_data[j] = data[k];
		}

		// Indices
		for (curarg = 7; curarg < 13; curarg++)
		{		
			data_array = mxGetCell(prhs[curarg],i);	
			data = (REAL *)mxGetData(data_array);		

			for (int k = 0; k < mxGetNumberOfElements(data_array); j++, k++)
				stacked_data[j] = data[k]-1;
		}

	

		add_edge(mrf, nodes[n1], nodes[n2],	
			 	 alpha,    lambda,	stacked_data);
	}

	mxFree(stacked_data);

	// Function below is optional - it may help if, for example, nodes are added in a random order
	mrf->SetAutomaticOrdering();
		
	typename MRFEnergy<TYPE>::Options mrf_options;
	mrf_options.m_iterMax = options[4]; // maximum number of iterations

	REAL energy, lowerBound = 0;


	if (options[2]) {
		mexPrintf("Graph loaded. Starting optimization using TRW-S.\n");
		/////////////////////// TRW-S algorithm //////////////////////
		mrf->Minimize_TRW_S(mrf_options, lowerBound, energy);
	} else {
		mexPrintf("Graph loaded. Starting optimization using BP.\n");
		//////////////////////// BP algorithm ////////////////////////
		mrf->ZeroMessages(); // in general not necessary - it may be faster to start 
		// with messages computed in previous iterations
		mrf->Minimize_BP(mrf_options, energy);
	}

	// Read solution
	plhs[0] = mxCreateNumericMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxUINT16_CLASS, mxREAL);
	uint16_t *L = (uint16_t *)mxGetData(plhs[0]);
	for (int i = 0; i < n_nodes; i++ )
		L[i] = (uint16_t)mrf->GetSolution(nodes[i]) + 1;

	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleScalar((double)energy);
		if (nlhs > 2) {
			plhs[2] = mxCreateDoubleScalar((double)lowerBound);
		}
	}

	// Clean up
	delete nodes;
	delete mrf;
}
