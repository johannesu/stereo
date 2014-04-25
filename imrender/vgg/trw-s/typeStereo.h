/******************************************************************
 * typeStereo.h
 *- Johannes
 *******************************************************************/

void flush_output();
#define debug() { mexPrintf("Reached line %d in file %s. \n", __LINE__, __FILE__);	flush_output(); }

#include <string>
#include <stdexcept>
using std::string;
using std::runtime_error;
#include <sstream>
#define ASSERT(cond) if (!(cond)) { std::stringstream sout; \
									sout << "Error (file " << __FILE__ << ", line " << __LINE__ << "): " << #cond; \
                                    throw runtime_error(sout.str()); }


#ifndef __TypeStereo_H__
#define __TypeStereo_H__

#include <string.h>
#include <assert.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <math.h>

// Define types
#ifdef _MSC_VER
typedef unsigned __int16 uint16_t;
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

template <class T> class MRFEnergy;

class TypeStereo
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
        int		m_K; // number of labels
    };
    
    struct LocalSize // number of labels is stored at MRFEnergy::m_Kglobal
    {
        LocalSize(int K);
        
    private:
        friend struct Vector;
        friend struct Edge;
        int		m_K; // number of labels
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
        EdgeData(REAL lambda, REAL alpha, REAL* data);
        
    private:
        friend struct Vector;
        friend struct Edge;
        REAL	       m_lambda;
        REAL            m_alpha;
        REAL*	 m_stacked_data;
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
        REAL ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin);            // return min_k { this[k] }, set kMin
        REAL ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K);              // same as previous, but additionally set this[k] -= vMin (and kMin is not returned)
        
    private:
        friend struct Edge;
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
        
        // Truncate smooth cost
        REAL Smooth(REAL val);
        
    protected:
        
        REAL		m_lambda;
        REAL        m_alpha;
        
        // message
        Vector*		m_message;
    private:
        int     m_dir; // 0 if Swap() was called even number of times, 1 otherwise
        REAL	m_data[1]; // Will be made larger inside initilize
    };
};




//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Implementation ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

inline TypeStereo::GlobalSize::GlobalSize(int K)
{
    m_K = K;
}
inline TypeStereo::LocalSize::LocalSize(int K)
{
    m_K = K;
}

///////////////////// NodeData and EdgeData ///////////////////////

inline TypeStereo::NodeData::NodeData(REAL* data)
{
    m_data = data;
}

inline TypeStereo::EdgeData::EdgeData(REAL lambda, REAL alpha, REAL* stacked_data)
{
    m_alpha  = alpha;
    m_lambda = lambda;
    m_stacked_data = stacked_data;
}

///////////////////// Vector ///////////////////////

inline int TypeStereo::Vector::GetSizeInBytes(GlobalSize Kglobal, LocalSize K)
{
    if (K.m_K < 1)
    {
        return -1;
    }
    return K.m_K*sizeof(REAL);
}
inline void TypeStereo::Vector::Initialize(GlobalSize Kglobal, LocalSize K, NodeData data)
{
    memcpy(m_data, data.m_data, K.m_K*sizeof(REAL));
}

inline void TypeStereo::Vector::Add(GlobalSize Kglobal, LocalSize K, NodeData data)
{
    for (int k=0; k<K.m_K; k++)
    {
        m_data[k] += data.m_data[k];
    }
}

inline void TypeStereo::Vector::SetZero(GlobalSize Kglobal, LocalSize K)
{
    memset(m_data, 0, K.m_K*sizeof(REAL));
}

inline void TypeStereo::Vector::Copy(GlobalSize Kglobal, LocalSize K, Vector* V)
{
    memcpy(m_data, V->m_data, K.m_K*sizeof(REAL));
}

inline void TypeStereo::Vector::Add(GlobalSize Kglobal, LocalSize K, Vector* V)
{
    for (int k=0; k<K.m_K; k++)
    {
        m_data[k] += V->m_data[k];
    }
}

inline TypeStereo::REAL TypeStereo::Vector::GetValue(GlobalSize Kglobal, LocalSize K, Label k)
{
    assert(k>=0 && k<K.m_K);
    return m_data[k];
}

inline TypeStereo::REAL TypeStereo::Vector::ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin)
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

inline TypeStereo::REAL TypeStereo::Vector::ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K)
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

inline int TypeStereo::Edge::GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data)
{
    int messageSizeInBytes = ((Ki.m_K > Kj.m_K) ? Ki.m_K : Kj.m_K)*sizeof(REAL);
    
    int size_of_data = 8*Kglobal.m_K*sizeof(REAL);
    return sizeof(Edge) - sizeof(REAL) + size_of_data + messageSizeInBytes;
    
    
}

inline int TypeStereo::Edge::GetBufSizeInBytes(int vectorMaxSizeInBytes)
{
    // Need to keep track of "data" and "reg" cost
    return 2*vectorMaxSizeInBytes;
}

inline void TypeStereo::Edge::Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj)
{
    m_lambda    = data.m_lambda;
    m_alpha     = data.m_alpha;
    int size_of_data = 8*Kglobal.m_K*sizeof(REAL);
    
    ((Edge*)this)->m_dir = 0;
    memcpy(((Edge*)this)->m_data, data.m_stacked_data, size_of_data);
    m_message = (Vector*)((char*)this + sizeof(Edge) - sizeof(REAL) + size_of_data);
}

inline TypeStereo::Vector* TypeStereo::Edge::GetMessagePtr()
{
    return m_message;
}

inline void TypeStereo::Edge::Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj)
{
    ((Edge*)this)->m_dir = 1 - ((Edge*)this)->m_dir;
}


inline TypeStereo::REAL TypeStereo::Edge::Smooth(REAL val)
{
        return ( m_alpha*std::min( std::abs(val) , m_lambda) );
}

inline TypeStereo::REAL TypeStereo::Edge::UpdateMessage(GlobalSize Kglobal,
        LocalSize Ksource,
        LocalSize Kdest,
        Vector* source,
        REAL gamma,
        int dir,
        void* _buf)
{                  
    REAL* buf = (REAL*) _buf;
    REAL* my_m = buf + Kglobal.m_K;     
      
    REAL * q;
    REAL * qprim;
    REAL * q_order;
    REAL * qprim_prev_id;
    REAL * qprim_next_id;
     
    
    if (dir == m_dir)
    {
        q               = m_data + 1*Kglobal.m_K;
        qprim           = m_data + 0*Kglobal.m_K;
        
        // This should be stored as int in the future
        q_order         = m_data + 3*Kglobal.m_K;
        qprim_prev_id   = m_data + 6*Kglobal.m_K;
        qprim_next_id   = m_data + 7*Kglobal.m_K;  
    } else
    {       
        q               = m_data + 0*Kglobal.m_K;
        qprim           = m_data + 1*Kglobal.m_K;
        
        // This should be stored as int in the future
        q_order         = m_data + 2*Kglobal.m_K;
        qprim_prev_id   = m_data + 4*Kglobal.m_K;
        qprim_next_id   = m_data + 5*Kglobal.m_K;
    }
  
    // Two pass into buffer
    int k;
    buf[(int)q_order[0]] = gamma*source->m_data[(int)q_order[0]] - m_message->m_data[(int)q_order[0]];
       
    REAL vMax =  buf[(int)q_order[0]];  // yes max.
    REAL vMin = std::numeric_limits<REAL>::max();
            
    REAL dist;  
    for (k=1;  k < Kglobal.m_K; k++)
    {       
        // Keep my cost or choose label of closest neighbor and pay
        // the regularization cost
        dist = q[(int)q_order[k]] - q[(int)q_order[k-1]]; 
        buf[(int)q_order[k]] = gamma*source->m_data[(int)q_order[k]] - m_message->m_data[(int)q_order[k]];
     
        // Find lowest D_i
        vMax = std::min(vMax,  buf[(int)q_order[k]]);
        
        // Choose smallest and update d and reg buffers (Omitted r_buf from LHS since == 0)
        buf[(int)q_order[k]] = std::min( buf[(int)q_order[k]], 
                                    buf[(int)q_order[k-1]] + m_alpha*dist);
    }
    
    // Add truncation cost
    vMax += m_alpha*m_lambda;
    
    k--;
    for (k--; k>=0; k--)
    {
        dist = q[(int)q_order[k+1]] - q[(int)q_order[k]];
        buf[(int)q_order[k]] = std::min( buf[(int)q_order[k]],
                                    buf[(int)q_order[k+1]] + m_alpha*dist
                                  ); 
    }
  
    // In this last single pass almost everything is calculated:
    // 1) Find costs
    // 2) Truncate
    // 3) Calculate vMin
    REAL prev_dist;
    REAL next_dist;
    
    for (k=0; k < Kglobal.m_K; k++)
    {
        // Choose smallest of two neighbors
        // Note qprim_prev/next points to
        // to the two q values which are closest.      
        prev_dist   = std::abs(qprim[k] - q[(int)qprim_prev_id[k]]);
        next_dist   = std::abs(qprim[k] - q[(int)qprim_next_id[k]]);
        
        // my_m[k] = m_message->m_data[k]
        // Choosing smallest cost               
        m_message->m_data[k] = std::min ( 
                buf[ (int)qprim_prev_id[k] ] + m_alpha*prev_dist,
                buf[ (int)qprim_next_id[k] ] + m_alpha*next_dist
                );
       
        // Truncate, highest possible cost for a message.
        if (m_message->m_data[k] > vMax)
            m_message->m_data[k] = vMax;
      
                        
        // vMin
        if (vMin > m_message->m_data[k])
            vMin = m_message->m_data[k];       
    }
    
    // Subtract vMin
    for (k = 0; k < Kglobal.m_K; k++)
    {
        m_message->m_data[k] -= vMin;
    }
    
    
//     // RESTORE points to be removed
//     for (int i = 0; i < Kglobal.m_K; i++)
//         q[i]  =     (REAL)m_data[0*Kglobal.m_K+ i]; 
//   
//     for (int i = 0; i < Kglobal.m_K; i++)
//         qprim[i]  = (REAL)m_data[1*Kglobal.m_K+ i]; 
//     
//     // Orders
//     for (int i = 0; i < Kglobal.m_K; i++)
//         q_order[i]          = (int)m_data[2*Kglobal.m_K + i];
//     
//     for (int i = 0; i < Kglobal.m_K; i++)
//         qprim_prev_id[i]    = (int)m_data[4*Kglobal.m_K + i];
//     
//     for (int i = 0; i < Kglobal.m_K; i++)
//         qprim_next_id[i]    = (int)m_data[5*Kglobal.m_K + i];
//     
//     // Straight copy of type General
//     Vector* v_buf = (Vector*) _buf;
//        
//     int ksource, kdest;   
//     for (ksource=0; ksource<Ksource.m_K; ksource++)
//     {
//         v_buf->m_data[ksource] = gamma*source->m_data[ksource] - m_message->m_data[ksource];
//     }
//     
//     // The key is here.
//     if (dir == ((Edge*)this)->m_dir)
//     {
//         for (kdest=0; kdest<Kdest.m_K; kdest++)
//         {
//             // Varying of qprim
//             vMin = v_buf->m_data[0] +  Smooth(qprim[0]  -  q[kdest]);
//             
//             for (ksource=1; ksource<Ksource.m_K; ksource++)
//             {
//                 vMin = std::min ( vMin,
//                                   v_buf->m_data[ksource] + Smooth(qprim[ksource] - q[kdest])
//                                  );
//             }
//             
//             m_message->m_data[kdest] = vMin;
//         }
//     }
//     else
//     {
//         for (kdest=0; kdest<Kdest.m_K; kdest++)
//         {
//             // Varying of q
//             vMin = v_buf->m_data[0] + Smooth(qprim[kdest] - q[0]);
//             
//             for (ksource=1; ksource<Ksource.m_K; ksource++)
//             {
//                 vMin = std::min ( vMin,
//                                   v_buf->m_data[ksource] + Smooth(qprim[kdest] - q[ksource])
//                                  );
//             }
//             
//             m_message->m_data[kdest] = vMin;
//         }
//     }
//     
//     // Find vMin
//     vMin = m_message->m_data[0];
//     for (kdest=1; kdest<Kdest.m_K; kdest++)
//     {
//         if (vMin > m_message->m_data[kdest])
//         {
//             vMin = m_message->m_data[kdest];
//         }
//     }
//     
//     // Subtract
//     for (k=0; k < Kglobal.m_K; k++)
//     {
//         m_message->m_data[k] -= vMin;
//     }
//     
//   
//     int VERBOSE = 0;
//     
//     for (int i = 0; i <3; i++)
//     {
//         REAL dif = abs( (m_message->m_data[i] + vMin) - my_m[i]);
//         
//         if (dif > 10)
//         {    
//             VERBOSE = 1;
//             mexPrintf("k =%d, diff = %g \n",i, dif  );
//         }
//     }
//         
//     if (VERBOSE == 1)
//     {
//         
//          
//          if (dir == m_dir)
//          {
//             mexPrintf("dir == m_dir\n");   
//          } else
//          {
//            mexPrintf("dir != m_dir\n");
//          }
// 
//          mexPrintf("smooth = @(F) %g *min(abs(F),%g) \n",m_alpha,m_lambda);
//          mexPrintf("vMax = %g \n", vMax);
//          
//              
//                   
//          for (int i = 0; i <3; i++)
//              mexPrintf("q[%d] = %g \n", i, q[i]);
//          
//          for (int i = 0; i <3; i++)
//              mexPrintf("qprim[%d] = %g \n", i, qprim[i]);
//          
//                  
//         
//          for (int i = 0; i <3; i++)
//              mexPrintf("m_message->m_data[%d] + vMin = %g\n", i , m_message->m_data[i] + vMin);
//          
//          for (int i = 0; i <3; i++)
//              mexPrintf("my_m[%d] = %g \n", i , my_m[i]);
//          
//              mexErrMsgTxt("stop \n");
// 
//     }
     
  
    return vMin;
}

// If dir==0, then sets dest[kj] += V(ksource,kj).
// If dir==1, then sets dest[ki] += V(ki,ksource).
inline void TypeStereo::Edge::AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir)
{
    REAL* m_data = ((Edge*)this)->m_data;
    
    REAL * q;
    REAL * qprim;
        
	q = m_data;
	qprim = m_data + Kglobal.m_K;
 
    
    // Vary over qprim
    int k;
    
    if (dir == ((Edge*)this)->m_dir)
    {
        for (k=0; k<Kdest.m_K; k++)
        {
            dest->m_data[k] +=  Smooth(qprim[ksource] - q[k]);
        }
    }
    // Vary of q
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
