/******************************************************************
typeStereo.h
// TODO: In Woodford et al. we acctually deal with rounded energies
// TODO: More effective EDGE::initilize instead of shuffling memory around
- Johannes


Energy function with truncated linear interactions:
   E(x)   =   \sum_i D_i(x_i)   +   \sum_ij V_ij(x_i,x_j)
   where x_i \in {0, 1, ..., K-1},
   V_ij(ki, kj) = min { alpha_ij*|ki-kj|, lambda_ij }.
   alpha_ij and lambda_ij must be non-negative.

Example usage:

Minimize function E(x,y) = Dx(x) + Dy(y) + min { alpha*|x - y| , lambda } where 
  x,y \in {0,1,2},
  Dx(0) = 0, Dx(1) = 1, Dx(2) = 2,
  Dy(0) = 3, Dy(1) = 4, Dy(2) = 5,
  alpha = 6,
  lambda = 7




#include <stdio.h>
#include "MRFEnergy.h"

void testTruncatedLinear()
{
	MRFEnergy<TypeStereo>* mrf;
	MRFEnergy<TypeStereo>::NodeId* nodes;
	MRFEnergy<TypeStereo>::Options options;
	TypeStereo::REAL energy, lowerBound;

	const int nodeNum = 2; // number of nodes
	const int K = 3; // number of labels
	TypeStereo::REAL D[K];
	int x, y;

	mrf = new MRFEnergy<TypeStereo>(TypeStereo::GlobalSize(K));
	nodes = new MRFEnergy<TypeStereo>::NodeId[nodeNum];

	// construct energy
	D[0] = 0; D[1] = 1; D[2] = 2;
	nodes[0] = mrf->AddNode(TypeStereo::LocalSize(), TypeStereo::NodeData(D));
	D[0] = 3; D[1] = 4; D[2] = 5;
	nodes[1] = mrf->AddNode(TypeStereo::LocalSize(), TypeStereo::NodeData(D));
	mrf->AddEdge(nodes[0], nodes[1], TypeStereo::EdgeData(6, 7));

	// Function below is optional - it may help if, for example, nodes are added in a random order
	// mrf->SetAutomaticOrdering();

	/////////////////////// TRW-S algorithm //////////////////////
	options.m_iterMax = 30; // maximum number of iterations
	mrf->Minimize_TRW_S(options, lowerBound, energy);

	// read solution
	x = mrf->GetSolution(nodes[0]);
	y = mrf->GetSolution(nodes[1]);

	printf("Solution: %d %d\n", x, y);

	//////////////////////// BP algorithm ////////////////////////
	mrf->ZeroMessages(); // in general not necessary - it may be faster to start 
	                     // with messages computed in previous iterations

	options.m_iterMax = 30; // maximum number of iterations
	mrf->Minimize_BP(options, energy);

	// read solution
	x = mrf->GetSolution(nodes[0]);
	y = mrf->GetSolution(nodes[1]);

	printf("Solution: %d %d\n", x, y);

	// done
	delete nodes;
	delete mrf;
}

*******************************************************************/

#ifndef __TYPESTEREO_H__
#define __TYPESTEREO_H__

#include <string.h>
#include <assert.h>

void flush_output();
#define debug() { mexPrintf("Reached line %d in file %s. \n", __LINE__, __FILE__);	flush_output(); }

template <class T> class MRFEnergy;
class TypeStereo
{
public:
	// types declarations
	struct Edge; // stores edge information and either forward or backward message
	struct Vector; // node parameters and messages
	typedef int Label;
	typedef double REAL;
	struct GlobalSize; // global information about number of labels
	struct LocalSize; // local information about number of labels (stored at each node)
	struct NodeData; // argument to MRFEnergy::AddNode()
	struct EdgeData; // argument to MRFEnergy::AddEdge()


	struct GlobalSize
	{
		GlobalSize(int K);

	private:
	friend struct Vector;
	friend struct Edge;
		int		m_K; // number of labels
	};

	struct LocalSize // number of labels is stored at MRFEnergy::m_Kglobal
	{
	};

	struct NodeData
	{
		NodeData(REAL* data); // data = pointer to array of size MRFEnergy::m_Kglobal

	private:
		friend struct Vector;
		friend struct Edge;
		REAL*		m_data;
	};

	struct EdgeData
	{
		EdgeData(REAL alpha, REAL lambda, REAL * q_distance);
/*
EdgeData(REAL alpha, REAL lambda,
			 int * q_order, REAL * q_distance, 
			 int * qprim_prev_id, int * qprim_next_id,
			 REAL * qprim_prev_dist, REAL * qprim_next_dist);

*/
	private:
		friend struct Vector;
		friend struct Edge;
		REAL		m_alpha;
		REAL		m_lambda;
		
		REAL * 		m_q_distance;
		
		//int * 		m_q_order;
/*
		int * 		m_qprim_prev_id;
		int * 		m_qprim_next_id;

		
		REAL * 		m_qprim_prev_dist;
		REAL * 		m_qprim_next_dist;
*/
	};




	//////////////////////////////////////////////////////////////////////////////////
	////////////////////////// Visible only to MRFEnergy /////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////


friend class MRFEnergy<TypeStereo>;

	struct Vector
	{
		static int GetSizeInBytes(GlobalSize Kglobal, LocalSize K); // returns -1 if invalid K's
		void Initialize(GlobalSize Kglobal, LocalSize K, NodeData data);  // called once when user adds a node
		void Add(GlobalSize Kglobal, LocalSize K, NodeData data); // called once when user calls MRFEnergy::AddNodeData()

		void SetZero(GlobalSize Kglobal, LocalSize K);                            // set this[k] = 0
		void Copy(GlobalSize Kglobal, LocalSize K, Vector* V);                    // set this[k] = V[k]
		void Add(GlobalSize Kglobal, LocalSize K, Vector* V);                     // set this[k] = this[k] + V[k]
		REAL GetValue(GlobalSize Kglobal, LocalSize K, Label k);                  // return this[k]
		REAL ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin);            // return vMin = min_k { this[k] }, set kMin
		REAL ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K);              // same as previous, but additionally set this[k] -= vMin (and kMin is not returned)

	// Additional data is stored here
	private:
	friend struct Edge;
	//	int 		m_q_order[1]; 
		REAL		m_data[1]; // actual size is MRFEnergy::m_Kglobal
	};

	struct Edge
	{
		static int GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data); // returns -1 if invalid data
		static int GetBufSizeInBytes(int vectorMaxSizeInBytes); // returns size of buffer need for UpdateMessage()
		void Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj); // called once when user adds an edge
		Vector* GetMessagePtr();
		void Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj); // if the client calls this function, then the meaning of 'dir'
								                                               // in distance transform functions is swapped

		// When UpdateMessage() is called, edge contains message from dest to source.
		// The function must replace it with the message from source to dest.
		// The update rule is given below assuming that source corresponds to tail (i) and dest corresponds
		// to head (j) (which is the case if dir==0).
		//
		// 1. Compute Di[ki] = gamma*source[ki] - message[ki].  (Note: message = message from j to i).
		// 2. Compute distance transform: set
		//       message[kj] = min_{ki} (Di[ki] + V(ki,kj)). (Note: message = message from i to j).
		// 3. Compute vMin = min_{kj} m_message[kj].
		// 4. Set m_message[kj] -= vMin.
		// 5. Return vMin.
		//
		// If dir==1 then source corresponds to j, sink corresponds to i. Then the update rule is
		//
		// 1. Compute Dj[kj] = gamma*source[kj] - message[kj].  (Note: message = message from i to j).
		// 2. Compute distance transform: set
		//       message[ki] = min_{kj} (Dj[kj] + V(ki,kj)). (Note: message = message from j to i).
		// 3. Compute vMin = min_{ki} m_message[ki].
		// 4. Set m_message[ki] -= vMin.
		// 5. Return vMin.
		//
		// If Edge::Swap has been called odd number of times, then the meaning of dir is swapped.
		//
		// Vector 'source' must not be modified. Function may use 'buf' as a temporary storage.
		REAL UpdateMessage(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Vector* source, REAL gamma, int dir, void* buf);

		// If dir==0, then sets dest[kj] += V(ksource,kj).
		// If dir==1, then sets dest[ki] += V(ki,ksource).
		// If Swap() has been called odd number of times, then the meaning of dir is swapped.
		void AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir);

	private:
		// edge information
		REAL		m_alpha;
		REAL		m_lambda;

		REAL 		m_data[1];
		// Added for distance transform
		// actual size is MRFEnergy::m_Kglobal
		// THESE SHOULD PROBABLY BE MOVED!
/*		
		REAL * m_q_distance;

		int  * m_qprim_prev_id;
		int  * m_qprim_next_id;
		REAL * m_qprim_prev_dist;
		REAL * m_qprim_next_dist;
*/

		// message
		Vector		m_message;
	};
};

//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Implementation ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

inline TypeStereo::GlobalSize::GlobalSize(int K)
{
	m_K = K;
}

///////////////////// NodeData and EdgeData ///////////////////////

inline TypeStereo::NodeData::NodeData(REAL* data)
{
	m_data = data;
}

/*
inline TypeStereo::EdgeData::EdgeData(
	REAL alpha, REAL lambda, 
	int * q_order, 
	REAL * q_distance, 
	int * qprim_prev_id, 
	int * qprim_next_id,
	REAL * qprim_prev_dist, 
	REAL * qprim_next_dist)
*/

inline TypeStereo::EdgeData::EdgeData(REAL alpha, REAL lambda, REAL * q_distance)
{
	m_alpha = alpha;
	m_lambda = lambda;

	m_q_distance = q_distance;
	
//	m_q_order = q_order;

	/*	
	

	m_qprim_prev_id = qprim_prev_id;
	m_qprim_next_id  = qprim_next_id;
	m_qprim_prev_dist = qprim_prev_dist;
	m_qprim_next_dist = qprim_next_dist;
	*/
}

///////////////////// Vector ///////////////////////

inline int TypeStereo::Vector::GetSizeInBytes(GlobalSize Kglobal, LocalSize K)
{
	if (Kglobal.m_K < 1)
	{
		mexPrintf("Kglobal.m_K < 1 line %d of %s \n", __LINE__, __FILE__);
		return -1;
	}
	return Kglobal.m_K*sizeof(REAL);
}
inline void TypeStereo::Vector::Initialize(GlobalSize Kglobal, LocalSize K, NodeData data)
{
	memcpy(m_data, data.m_data, Kglobal.m_K*sizeof(REAL));
}

inline void TypeStereo::Vector::Add(GlobalSize Kglobal, LocalSize K, NodeData data)
{
	for (int k=0; k<Kglobal.m_K; k++)
	{
		m_data[k] += data.m_data[k];
	}
}

inline void TypeStereo::Vector::SetZero(GlobalSize Kglobal, LocalSize K)
{
	memset(m_data, 0, Kglobal.m_K*sizeof(REAL));
}

inline void TypeStereo::Vector::Copy(GlobalSize Kglobal, LocalSize K, Vector* V)
{
	memcpy(m_data, V->m_data, Kglobal.m_K*sizeof(REAL));
}

inline void TypeStereo::Vector::Add(GlobalSize Kglobal, LocalSize K, Vector* V)
{
	for (int k=0; k<Kglobal.m_K; k++)
	{
		m_data[k] += V->m_data[k];
	}
}

inline TypeStereo::REAL TypeStereo::Vector::GetValue(GlobalSize Kglobal, LocalSize K, Label k)
{
	assert(k>=0 && k<Kglobal.m_K);
	return m_data[k];
}

inline TypeStereo::REAL TypeStereo::Vector::ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin)
{
	REAL vMin = m_data[0];
	kMin = 0;
	for (int k=1; k<Kglobal.m_K; k++)
	{
		if (vMin > m_data[k])
		{
			vMin = m_data[k];
			kMin = k;
		}
	}

	return vMin;
}

inline TypeStereo::REAL TypeStereo::Vector::ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K)
{
	REAL vMin = m_data[0];
	for (int k=1; k<Kglobal.m_K; k++)
	{
		if (vMin > m_data[k])
		{
			vMin = m_data[k];
		}
	}
	for (int k=0; k<Kglobal.m_K; k++)
	{
		m_data[k] -= vMin;
	}

	return vMin;
}

///////////////////// EdgeDataAndMessage implementation /////////////////////////

inline int TypeStereo::Edge::GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data)
{
	int messageSizeInBytes = (Kglobal.m_K)*sizeof(REAL);

	return (sizeof(Edge) - sizeof(REAL) + (Kglobal.m_K)*sizeof(REAL) + messageSizeInBytes);

	//return sizeof(Edge);
	//return 3*sizeof(int)*Kglobal.m_K + 4*sizeof(REAL)*Kglobal.m_K+2*sizeof(REAL);
}

inline int TypeStereo::Edge::GetBufSizeInBytes(int vectorMaxSizeInBytes)
{
	// Need as much space as number of labels*sizeof(REAL).
	return vectorMaxSizeInBytes;
}

// TODO: Change to meaningfull names!
// !!!! Right now it's copying to m_data !!!!

inline void TypeStereo::Edge::Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj)
{
	// Normal copy
	m_alpha = data.m_alpha;
	m_lambda = data.m_lambda;
	

	// Be very careful here...
	int int_data_size = Kglobal.m_K*sizeof(int);
	int real_data_size = Kglobal.m_K*sizeof(REAL);

	// Copying vectors
	memcpy( ((Edge*)this)->m_data, data.m_q_distance, real_data_size);

	// Move the m_message pointer to correct location.	
	// Maye after removing dir sizeof(REAL) might not need to be subtracted
	m_message = (Vector*)((char*)this + sizeof(Edge) - sizeof(REAL) + real_data_size);
/*
	m_q_order = data.m_q_order;
	m_q_distance = data.m_q_distance;

	m_qprim_prev_id =  data.m_qprim_prev_id;
	m_qprim_next_id  = data.m_qprim_next_id;
	m_qprim_prev_dist =  data.m_qprim_prev_dist;
	m_qprim_next_dist =  data.m_qprim_next_dist;
*/
}

inline TypeStereo::Vector* TypeStereo::Edge::GetMessagePtr()
{
	return &m_message;
}

inline void TypeStereo::Edge::Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj)
{
}

inline TypeStereo::REAL TypeStereo::Edge::UpdateMessage(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Vector* source, REAL gamma, int dir, void* _buf)
{
	mexPrintf("Inside update message \n");
	flush_output();
	int k;
	REAL vMin;
	REAL* buf = (REAL*) _buf; // Store intial pass over qs.
	
	debug();
	mexPrintf("lambda = %g \n", m_lambda);
	mexPrintf("alpha = %g \n", m_alpha);
	for (int k = 0; k < 3; k++)
	{
		mexPrintf("id = %d\n",k);
		//mexPrintf("m_q_order[%d] = %d\n", k, m_q_order[k]);
		
				
		//mexPrintf("m_qprim_prev_id[%d] = %d\n", k, m_qprim_prev_id[k]);
		//mexPrintf("m_qprim_next_id[%d] = %d\n", k, m_qprim_next_id[k]);
		mexPrintf("m_q_distance[%d] = %g \n", k, m_q_distance[k]);
		//mexPrintf("m_qprim_prev_dist[%d] = %g \n", k, m_qprim_prev_dist[k]);
		//mexPrintf("m_qprim_next_dist[%d] = %g \n", k, m_qprim_next_dist[k]);
		//mexPrintf("\n\n");
	}
	debug();
/*
	mexPrintf("gamma*source->m_data[m_q_order[0]] = %g \n", gamma*source->m_data[m_q_order[0]]);
	debug();
	
	mexPrintf("m_message.m_data[m_q_order[0]] = %g \n", m_message.m_data[m_q_order[0]]);
	debug();
		
	mexPrintf("buf[0] = %g \n", buf[0]);
	debug();


	buf[0] = gamma*source->m_data[m_q_order[0]] - m_message.m_data[m_q_order[0]];
    	debug();
   	 // Two pass lower envelope
	 // This updates the cost at the qs
	for (k=1; k<Kglobal.m_K; k++)
	{
		buf[k] = gamma*source->m_data[m_q_order[k]] - m_message.m_data[m_q_order[k]];
		if (buf[k] > buf[k-1] + m_alpha*m_q_distance[k-1])
		{
			buf[k] = buf[k-1] + m_alpha*m_q_distance[k-1];
		} 
	}

	debug();

	// Passing backwards
	for (k--; k>=0; k--)
	{
		if (buf[k] > buf[k+1] + m_alpha*m_q_distance[k])
		{
			buf[k] = buf[k+1] + m_alpha*m_q_distance[k]; // First element in distance is distance between 0 and 1.
		}
	}
	
	debug();

	// Now it's time to
	// 1) Update the message this is a single pass
	// 2) Truncate
	// 3) vMin
	vMin = m_lambda;
	for (k=0; k<Kglobal.m_K; k++)
	{
		// 1)
		m_message.m_data[k] = buf[m_qprim_prev_id[k]] + m_alpha*m_qprim_prev_dist[k];
		
		if (m_message.m_data[k] > buf[m_qprim_next_id[k]] + m_alpha*m_qprim_next_dist[k])
				m_message.m_data[k] = buf[m_qprim_next_id[k]] + m_alpha*m_qprim_next_dist[k];

		// 2)
		if (m_message.m_data[k] > m_lambda)
			m_message.m_data[k] = m_lambda;
		
		// 3)
		if (vMin > m_message.m_data[k])
			vMin = m_message.m_data[k];		
	}

	debug();

	// Subtract vMin from each message
	for (k=0; k<Kglobal.m_K; k++)
	{
		m_message.m_data[k] -= vMin;
	}

	return vMin;
*/
	return 0;
}

inline void TypeStereo::Edge::AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir)
{
	assert(ksource>=0 && ksource<Kglobal.m_K);

	int k;

	for (k=0; k<ksource; k++)
	{
		dest->m_data[k] += (ksource-k)*m_alpha < m_lambda ? (ksource-k)*m_alpha : m_lambda;
	}
	for (k++; k<Kglobal.m_K; k++)
	{
		dest->m_data[k] += m_alpha*(k-ksource) < m_lambda ? m_alpha*(k-ksource) : m_lambda;
	}
}

//////////////////////////////////////////////////////////////////////////////////

#endif
