/******************************************************************
 * typeStereoQuadratic.h
 *- Johannes Ulén 2013
 *******************************************************************/

#include <string>
#include <stdexcept>
using std::string;
using std::runtime_error;
#include <sstream>

#ifndef __TypeStereoQuadratic_H__
#define __TypeStereoQuadratic_H__

#include <string.h>
#include <assert.h>
#include <algorithm>
#include <limits>
using std::min;

template <class T> class MRFEnergy;

class TypeStereoQuadratic
{
public:
    struct Edge; // stores edge information and either forward or backward message
    struct Vector; // node parameters and messages

    // types declarations
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
        int     m_K; // number of labels
    };

    struct LocalSize // number of labels is stored at MRFEnergy::m_Kglobal
    {
        LocalSize(int K);

    private:
        friend struct Vector;
        friend struct Edge;
        int     m_K; // number of labels
    };

    struct NodeData
    {
        NodeData(REAL* data); // data = pointer to array of size MRFEnergy::m_Kglobal

    private:
        friend struct Vector;
        friend struct Edge;
        REAL*       m_data;
    };

    struct EdgeData
    {
        EdgeData(REAL lambda, REAL alpha, REAL* data, int* indices);

    private:
        friend struct Vector;
        friend struct Edge;
        REAL    m_lambda;
        REAL    m_alpha;
        REAL*   m_stacked_data;
        int*    m_stacked_indices;
    };

    //////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// Visible only to MRFEnergy /////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    friend class MRFEnergy<TypeStereoQuadratic>;

    struct Vector
    {
        static int GetSizeInBytes(GlobalSize Kglobal, LocalSize K); // returns -1 if invalid K's
        void Initialize(GlobalSize Kglobal, LocalSize K, NodeData data);  // called once when user adds a node
        void Add(GlobalSize Kglobal, LocalSize K, NodeData data); // called once when user calls MRFEnergy::AddNodeData()

        void SetZero(GlobalSize Kglobal, LocalSize K);                            // set this[k] = 0
        void Copy(GlobalSize Kglobal, LocalSize K, Vector* V);                    // set this[k] = V[k]
        void Add(GlobalSize Kglobal, LocalSize K, Vector* V);                     // set this[k] = this[k] + V[k]
        REAL GetValue(GlobalSize Kglobal, LocalSize K, Label k);                  // return this[k]
        REAL ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin);            // return min_k { this[k] }, set kMin
        REAL ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K);              // same as previous, but additionally set this[k] -= vMin (and kMin is not returned)

    private:
        friend struct Edge;
        REAL        m_data[1]; // actual size is MRFEnergy::m_Kglobal
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

        // Truncate smooth cost
        REAL Smooth(REAL val);

    protected:

        REAL        m_lambda;
        REAL        m_alpha;

        // message
        Vector*     m_message;
    private:
        int     m_dir; // 0 if Swap() was called even number of times, 1 otherwise
        int    * m_inds;
        REAL    m_data[1]; // Will be made larger inside initilize
    };
};




//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Implementation ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

inline TypeStereoQuadratic::GlobalSize::GlobalSize(int K)
{
    m_K = K;
}
inline TypeStereoQuadratic::LocalSize::LocalSize(int K)
{
    m_K = K;
}

///////////////////// NodeData and EdgeData ///////////////////////

inline TypeStereoQuadratic::NodeData::NodeData(REAL* data)
{
    m_data = data;
}

inline TypeStereoQuadratic::EdgeData::EdgeData(REAL lambda, REAL alpha, REAL* stacked_data, int* stacked_indices)
{
    m_alpha  = alpha;
    m_lambda = lambda;
    m_stacked_data = stacked_data;
    m_stacked_indices = stacked_indices;
}

///////////////////// Vector ///////////////////////

inline int TypeStereoQuadratic::Vector::GetSizeInBytes(GlobalSize Kglobal, LocalSize K)
{
    if (K.m_K < 1)
    {
        return -1;
    }
    return K.m_K*sizeof(REAL);
}
inline void TypeStereoQuadratic::Vector::Initialize(GlobalSize Kglobal, LocalSize K, NodeData data)
{
    memcpy(m_data, data.m_data, K.m_K*sizeof(REAL));
}

inline void TypeStereoQuadratic::Vector::Add(GlobalSize Kglobal, LocalSize K, NodeData data)
{
    for (int k=0; k<K.m_K; k++)
    {
        m_data[k] += data.m_data[k];
    }
}

inline void TypeStereoQuadratic::Vector::SetZero(GlobalSize Kglobal, LocalSize K)
{
    memset(m_data, 0, K.m_K*sizeof(REAL));
}

inline void TypeStereoQuadratic::Vector::Copy(GlobalSize Kglobal, LocalSize K, Vector* V)
{
    memcpy(m_data, V->m_data, K.m_K*sizeof(REAL));
}

inline void TypeStereoQuadratic::Vector::Add(GlobalSize Kglobal, LocalSize K, Vector* V)
{
    for (int k=0; k<K.m_K; k++)
    {
        m_data[k] += V->m_data[k];
    }
}

inline TypeStereoQuadratic::REAL TypeStereoQuadratic::Vector::GetValue(GlobalSize Kglobal, LocalSize K, Label k)
{
    assert(k>=0 && k<K.m_K);
    return m_data[k];
}

inline TypeStereoQuadratic::REAL TypeStereoQuadratic::Vector::ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin)
{
    REAL vMin = m_data[0];
    kMin = 0;
    for (int k=1; k<K.m_K; k++)
    {
        if (vMin > m_data[k])
        {
            vMin = m_data[k];
            kMin = k;
        }
    }

    return vMin;
}

inline TypeStereoQuadratic::REAL TypeStereoQuadratic::Vector::ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K)
{
    REAL vMin = m_data[0];
    for (int k=1; k<K.m_K; k++)
    {
        if (vMin > m_data[k])
        {
            vMin = m_data[k];
        }
    }
    for (int k=0; k<K.m_K; k++)
    {
        m_data[k] -= vMin;
    }

    return vMin;
}

///////////////////// EdgeDataAndMessage implementation /////////////////////////

inline int TypeStereoQuadratic::Edge::GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data)
{
    int messageSizeInBytes = ((Ki.m_K > Kj.m_K) ? Ki.m_K : Kj.m_K)*sizeof(REAL);

    int size_of_data = 2*Kglobal.m_K*sizeof(REAL);
    int size_of_indices = 2*Kglobal.m_K*sizeof(int);
    return sizeof(Edge) - sizeof(REAL) + size_of_data + size_of_indices + messageSizeInBytes;
}

inline int TypeStereoQuadratic::Edge::GetBufSizeInBytes(int vectorMaxSizeInBytes)
{
    int K = vectorMaxSizeInBytes / sizeof(REAL);
    return (2*K+1)*sizeof(REAL) + K*sizeof(int);

}

inline void TypeStereoQuadratic::Edge::Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj)
{

    m_lambda    = data.m_lambda;
    m_alpha     = data.m_alpha;
    int size_of_data = 2*Kglobal.m_K*sizeof(REAL);
    int size_of_indices = 2*Kglobal.m_K*sizeof(int);

    // Default m_dir
    ((Edge*)this)->m_dir = 0;

    // Take pointer to data move it 2 doubles and change to int
    // it now points at correct location<
    double * m_first = ((Edge*)this)->m_data;
    m_first += 2*Kglobal.m_K;
    ((Edge*)this)->m_inds = (int *) m_first;

    memcpy(((Edge*)this)->m_data, data.m_stacked_data, size_of_data);
    memcpy(((Edge*)this)->m_inds, data.m_stacked_indices, size_of_indices);

    m_message = (Vector*)((char*)this + sizeof(Edge) - sizeof(REAL) + size_of_data + size_of_indices);
}

inline TypeStereoQuadratic::Vector* TypeStereoQuadratic::Edge::GetMessagePtr()
{
    return m_message;
}

inline void TypeStereoQuadratic::Edge::Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj)
{
    ((Edge*)this)->m_dir = 1 - ((Edge*)this)->m_dir;
}


inline TypeStereoQuadratic::REAL TypeStereoQuadratic::Edge::Smooth(REAL val)
{
        return ( m_alpha* min( val * val , m_lambda) );
}

inline TypeStereoQuadratic::REAL TypeStereoQuadratic::Edge::UpdateMessage(GlobalSize Kglobal,
        LocalSize Ksource,
        LocalSize Kdest,
        Vector* source,
        REAL gamma,
        int dir,
        void* _buf)
{
    REAL *q;
    REAL *qprim;
    int *q_order;
    int *qprim_order;
    int *qprim_prev_id;
    int *qprim_next_id;

    // All indices
    if (dir == m_dir)
    {
        q               = m_data + 1*Kglobal.m_K;
        qprim           = m_data + 0*Kglobal.m_K;

        q_order         = m_inds + 1*Kglobal.m_K;
        qprim_order     = m_inds + 0*Kglobal.m_K;
    } else
    {
        q               = m_data + 0*Kglobal.m_K;
        qprim           = m_data + 1*Kglobal.m_K;

        q_order         = m_inds + 0*Kglobal.m_K;
        qprim_order     = m_inds + 1*Kglobal.m_K;
    }

    // buffer variables
    REAL* H = (REAL*) _buf;
    REAL* z = (REAL*) _buf + Kglobal.m_K;
    int * v = (int*) ( (REAL*)_buf +  2*Kglobal.m_K +1);

    // Distance transform very similar to
    // Efficent Belief Propgation for early Vision
    // Felzenszwalb and Huttenlocher.
    int j,k,l;
    REAL s, val;

    // Temporary variables optimization will take care of this
    REAL hk, hj, qk, qj, qprimk;

    // Keeping track on largest cost and smallest cost
    // (both initialized at +infty is correct).
    REAL vTrunc = std::numeric_limits<REAL>::infinity();
    REAL vMin = std::numeric_limits<REAL>::infinity();

    // Read position of "roots" of every parabola
    for (k = 0; k < Kglobal.m_K; k++)
    {
        H[k] = gamma*source->m_data[k] - m_message->m_data[k];
        vTrunc = min(vTrunc, H[k]);
    }

	// No regularization
	if (m_alpha == 0)
	{
		for (k = 0; k < Kglobal.m_K; k++)
				m_message->m_data[k] = vTrunc;

        vMin = vTrunc;

    // Regularization
	} else
	{
	    // Add truncated regularization cost
	    vTrunc += m_alpha*m_lambda;

        // Init
        j = 0;
        v[0] = q_order[0];
        z[0] = -std::numeric_limits<REAL>::infinity();
        z[1] = std::numeric_limits<REAL>::infinity();

	    for (k = 1; k <  Kglobal.m_K; k++)
	    {
		   // Current Height
	        hk = H[ q_order[k] ];
	        qk = q[ q_order[k] ];

	        // "goto row 6"
	        for (l = k; l >= 0; l--)
	        {
	            // Compare height
	            hj = H[ v[j] ];
	            qj = q[ v[j] ];
   
	            // If the point on the lower envelope is equal we choose the lowest one
	            if ( (qk - qj) < 1e-8)
	            {
	                // If hk > hj then qk will never be useful
	                // and identically hj > hk then qj never useful.
	                if (hj > hk)
	                {
	                    // Special case:
	                    // First parabola is replaced
	                    if (j == 0)
	                    {
	                        v[0] = q_order[k];
	                        z[0] = -std::numeric_limits<REAL>::infinity();
	                        z[1] = std::numeric_limits<REAL>::infinity();
	                        break;
	                    }
	                    else
	                    {
	                        // Search backwards and intersect with previous parabola
	                        j--;
	                    }
	                }
	                // This parabola is useless.
	                else
	                {
	                    break;
	                }
	            } 
	            else 
	            {
	                // Intersection point (c = alpha in the Efficent paper)
	                s = (+ (hk + m_alpha * qk*qk)
	                     - (hj + m_alpha * qj*qj)
	                    ) / ( 2*m_alpha*( qk - qj));

	                // parabolas intersect prior to current lower envelope boundary
	                // need to adjust z[j].
	                if (s <= z[j])
	                {
	                    j--;
	                }
	                // Adding new intersection
	                else
	                {
                        j++;
                        v[j] = q_order[k];
                        z[j] = s;
                        z[j+1] = std::numeric_limits<REAL>::infinity();
	                    break;
	                }
	            }
	        }
	    }

	    j=0;
	    for (k = 0; k <  Kglobal.m_K; k++)
	    {
	        // q' point
	        qprimk = qprim[ qprim_order[k] ];

	        // Find correct parabola
	        while( z[j+1] < qprimk ) { j++; };

	        // The base point of the parabola
	        qj = q[ v[j] ];
	        hj = H[ v[j] ];

            val = qprimk - qj;

	        m_message->m_data[ qprim_order[k] ] =  min( vTrunc, m_alpha*val*val + hj);
	    }

        for (k = 0; k < Kglobal.m_K; k++)
            vMin = min( vMin, m_message->m_data[k] );
	}

    // Subtract vMin
    for (k = 0; k < Kglobal.m_K; k++)
        m_message->m_data[k] -= vMin;

    return vMin;
}

// If dir==0, then sets dest[kj] += V(ksource,kj).
// If dir==1, then sets dest[ki] += V(ki,ksource).
inline void TypeStereoQuadratic::Edge::AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir)
{
    REAL* m_data = ((Edge*)this)->m_data;

    REAL * q;
    REAL * qprim;

    q = m_data;
    qprim = m_data + Kglobal.m_K;

    int k;

    if (dir == ((Edge*)this)->m_dir)
    {
        for (k=0; k<Kdest.m_K; k++)
        {
            dest->m_data[k] +=  Smooth(qprim[ksource] - q[k]);
        }
    }
    else
    {
        for (k=0; k<Kdest.m_K; k++)
        {
            dest->m_data[k] += Smooth(qprim[k] - q[ksource]);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////

#endif
