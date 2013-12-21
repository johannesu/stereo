
// Johannes Ulén, 2013
// MATLAB wrapper for Convergent Tree-reweighted Message Passing for Energy Minimization (TRW-S) by Vladimir Kolmogorov
// http://pub.ist.ac.at/~vnk/papers/TRW-S.html
//
// The wrapper is inspired by Olivier Woodfoord's imrender
// and uses it's modification to TRW-S to make it compile with MATLAB.
// http://www.robots.ox.ac.uk/~ojw/software.htm
#include "mex.h"
#include "utils/mexutils.h"
#include "utils/cppmatrix.h"
#include "trw-s/MRFEnergy.h"

typedef std::pair<double, int> Pair;

struct CmpPair
{
    bool operator()(const Pair& a, const Pair& b)
    { return a.first < b.first; }
};

// Catch errors
static void erfunc(char *err) {
	mexErrMsgTxt(err);
}

template<typename TYPE>
static void solve_mrf(int nlhs, mxArray *plhs[], int nrhs,  const mxArray  *prhs[])
{
  int curarg = 1;
  matrix<double> unary(prhs[curarg++]);
  matrix<unsigned int> connectivity(prhs[curarg++]);
  matrix<double> q(prhs[curarg++]);
  matrix<double> qprim(prhs[curarg++]);
  matrix<double> alphas(prhs[curarg++]);
  matrix<double> tol(prhs[curarg++]);
  double lambda = tol(0);

  MexParams params(nrhs-curarg, prhs+curarg); //Structure to hold and parse additional parameters
  double maxiter = params.get<double>("maxiter", 1000);
  double max_relgap = params.get<double>("max_relgap", 0);

  ASSERT(connectivity.M == 2);
  ASSERT(q.N == qprim.N);
  ASSERT(q.N == connectivity.N);

  ASSERT(unary.M == q.M);
  ASSERT(unary.M == qprim.M);

  ASSERT(alphas.M == q.N);
  ASSERT(alphas.N == 1);
  ASSERT(tol.numel() == 1);

  int n_labels = unary.M; // Number of labels
  int n_nodes = unary.N; // Number of unary terms
  int n_edges = q.N; // Number of pairwise terms

  // Create MRF instance
  typedef typename TYPE::REAL REAL;
  MRFEnergy<TYPE> *mrf = new MRFEnergy<TYPE>(typename TYPE::GlobalSize(n_labels),
                                             erfunc);
  typename MRFEnergy<TYPE>::NodeId * nodes =  new typename MRFEnergy<TYPE>::NodeId[n_nodes];

  // Add unary terms
  matrix<REAL> term(n_labels);
  for (int u = 0; u < n_nodes; u++)
  {
    for (int l = 0; l < n_labels; l++)
      term(l) = unary(l,u);

    nodes[u] = mrf->AddNode(typename TYPE::LocalSize(n_labels),
                            typename TYPE::NodeData( term.data ));
  }

  // Sort the points
  std::vector<Pair> q_pair;
  std::vector<Pair> qprim_pair;

  int number_of_indices = 2*n_labels;
  matrix<REAL>  stacked_data(number_of_indices);
  matrix<int>   stacked_indices(number_of_indices);

  // Add pairwise energies
  for (int p = 0; p < n_edges; p++)
  {
      double alpha = alphas(p);

      // Pair up and get ordering
      for (int j = 0; j < n_labels; j++)
      {
          q_pair.push_back(Pair(q(j,p), j));
          qprim_pair.push_back(Pair(qprim(j,p), j));

          // q_pair.second() contains the order
          std::sort(q_pair.begin(), q_pair.end(), CmpPair());
          std::sort(qprim_pair.begin(), qprim_pair.end(), CmpPair());
      }

      // Stack all data
      int stack_id = 0;
      for (int j = 0; j < q_pair.size(); j++, stack_id++)
      {
        stacked_data(stack_id) = q(j,p);
        stacked_indices(stack_id) = q_pair[j].second;
      }

      for (int j = 0; j < qprim_pair.size(); j++, stack_id++)
      {
        stacked_data(stack_id) = qprim(j,p);
        stacked_indices(stack_id) = qprim_pair[j].second;
      }

      mrf->AddEdge(nodes[connectivity(0,p)],
                   nodes[connectivity(1,p)],
                   typename TYPE::EdgeData(lambda, alpha, stacked_data.data, stacked_indices.data) );

      q_pair.clear();
      qprim_pair.clear();
  }

  mrf->SetAutomaticOrdering();

  // Options struct
  typename MRFEnergy<TYPE>::Options mrf_options;
  mrf_options.m_iterMax = (int)maxiter;
  mrf_options.m_relgapMax = (double) max_relgap; 
  mrf_options.m_printMinIter = false;

  // Solve
  REAL energy, lowerBound = 0;
  int iterations = mrf->Minimize_TRW_S(mrf_options, lowerBound, energy);

  // OUTPUT
  matrix<double> labelling(n_nodes);

  // Read solution
  for (int u = 0; u < n_nodes; u++ )
    labelling(u) = mrf->GetSolution(nodes[u])+1;

  plhs[0] = labelling;
  plhs[1] = mxCreateDoubleScalar((double)energy);
  plhs[2] = mxCreateDoubleScalar((double)lowerBound);
  plhs[3] = mxCreateDoubleScalar((double)iterations);

  delete nodes;
  delete mrf;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,  const mxArray  *prhs[])
{
  // Parse input
  ASSERT(nrhs == 8);
  ASSERT(nlhs == 4);
  matrix<int> kernel(prhs[0]);

  switch (kernel(0))
  {
    case 1: solve_mrf<TypeStereoLinear>(nlhs, plhs, nrhs, prhs);
        break;
    case 2 :solve_mrf<TypeStereoQuadratic>(nlhs, plhs, nrhs, prhs);
        break;
    default: mexErrMsgTxt("Unsupported kernel");
  }
}