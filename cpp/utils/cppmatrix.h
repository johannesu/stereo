//
// CPPMATRIX
//
// C++ wrapper for Matlab matrices
//
// Petter Strandmark 2010
//


#ifndef CPPMATRIX_MATLAB_HEADER
#define CPPMATRIX_MATLAB_HEADER

#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "mex.h"

#include <sstream>
#define ASSERT(cond) if (!(cond)) { std::stringstream sout; \
                                    sout << "Error (file " << __FILE__ \
                                    	 <<", line " << __LINE__ << "): " \
                                         << #cond; \
                                    throw runtime_error(sout.str()); }

// Extremely annoying macros need to be undefined in
// order to make std::min and std::max work.
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif


//
// typeToID<T>() returns the correct Matlab identifier
// for the type T. It throws an exception for unknown
// types.
//
namespace {
	template<typename T>
	mxClassID typeToID()
	{
		mexErrMsgTxt("Unknown type!");
		return mxUINT8_CLASS;
	}

	template<>
	mxClassID typeToID<double>()
	{
		return mxDOUBLE_CLASS;
	}
	template<>
	mxClassID typeToID<float>()
	{
		return mxSINGLE_CLASS;
	}
	template<>
	mxClassID typeToID<unsigned char>()
	{
		return mxUINT8_CLASS;
	}
	template<>
	mxClassID typeToID<signed char>()
	{
		return mxINT8_CLASS;
	}
	template<>
	mxClassID typeToID<unsigned int>()
	{
		ASSERT(sizeof(unsigned int)==4);
		return mxUINT32_CLASS;
	}
	template<>
	mxClassID typeToID<signed int>()
	{
		ASSERT(sizeof(signed int)==4);
		return mxINT32_CLASS;
	}
	template<>
	mxClassID typeToID<unsigned short>()
	{
		ASSERT(sizeof(unsigned short)==2);
		return mxUINT16_CLASS;
	}
	template<>
	mxClassID typeToID<signed short>()
	{
		ASSERT(sizeof(signed short)==2);
		return mxINT16_CLASS;
	}
	template<>
	mxClassID typeToID<unsigned long long>()
	{
		ASSERT(sizeof(unsigned long long)==8);
		return mxUINT64_CLASS;
	}
	template<>
	mxClassID typeToID<signed long long>()
	{
		ASSERT(sizeof(signed long long)==8);
		return mxINT64_CLASS;
	}
	template<>
	mxClassID typeToID<bool>()
	{
		return mxLOGICAL_CLASS;
	}
}


template<typename T>
class matrix
{
public:

	T* data;
    mxArray* array;
	mwSize M,N,O,P;

	matrix(const mxArray* array)
	{
		this->array = (mxArray*)array; //Hack to allow const mxArrays
		ASSERT(!mxIsSparse(array));
		ASSERT(mxGetClassID(array) == typeToID<T>());

		int ndim = mxGetNumberOfDimensions(array);
		ASSERT(ndim<=4);

		if (ndim <= 2) {
			M = mxGetM(array);
			N = mxGetN(array);
			O = 1;
			P = 1;
		}
		else if (ndim==3) {
			const mwSize* dims = mxGetDimensions(array);
			M = dims[0];
			N = dims[1];
			O = dims[2];
			P = 1;
		}
		else {
			const mwSize* dims = mxGetDimensions(array);
			M = dims[0];
			N = dims[1];
			O = dims[2];
			P = dims[3];
		}
		data = (T*)mxGetPr(array);

		shouldDestroy = false;
	}

	matrix(int M, int N=1, int O=1, int P=1)
	{
		this->M = M;
		this->N = N;
		this->O = O;
		this->P = P;
		if (O==1 && P==1) {
			mwSize size[] = {M,N};
			this->array = mxCreateNumericArray(2,size,typeToID<T>(),mxREAL);
		}
		else if (P==1) {
			mwSize size[] = {M,N,O};
			this->array = mxCreateNumericArray(3,size,typeToID<T>(),mxREAL);
		}
		else {
			mwSize size[] = {M,N,O,P};
			this->array = mxCreateNumericArray(4,size,typeToID<T>(),mxREAL);
		}
		data = (T*)mxGetPr(array);

		shouldDestroy = true;
	}

	matrix(const matrix& org)
	{
		*this = org;
	}

	matrix()
	{
		M = N = O = P = 0;
		data = 0;
		array = 0;
		shouldDestroy = false;
	}

	~matrix()
	{
		if (shouldDestroy) {
			mxDestroyArray(array);
		}
	}

	mwSize numel() const
	{
		return M*N*O*P;
	}

	void operator=(const matrix& org)
	{
		if (org.shouldDestroy) {
			throw std::runtime_error("matrix() : Cannot copy a managed matrix");
		}
		shouldDestroy = false;
		M = org.M;
		N = org.N;
		O = org.O;
		P = org.P;
		data = org.data;
		array = org.array;
	}

	T& operator[](mwSize i)
	{
		ASSERT(i>=0 && i<numel());
		return data[i];
	}
	T& operator()(mwSize i)
	{
		return operator[](i);
	}
	T& operator()(mwSize i, mwSize j)
	{
		return operator[](i + M*j);
	}
	T& operator()(mwSize i, mwSize j, mwSize k)
	{
		return operator[](i + M*j + M*N*k);
	}
	T& operator()(mwSize i, mwSize j, mwSize k, mwSize l)
	{
		return operator[](i + M*j + M*N*k + M*N*O*l);
	}

	T operator[](mwSize i) const
	{
		ASSERT(i>=0 && i<numel());
		return data[i];
	}
	T operator()(mwSize i) const
	{
		return operator[](i);
	}
	T operator()(mwSize i, mwSize j) const
	{
		return operator[](i + M*j);
	}
	T operator()(mwSize i, mwSize j, mwSize k) const
	{
		return operator[](i + M*j + M*N*k);
	}
	T operator()(mwSize i, mwSize j, mwSize k, mwSize l) const
	{
		return operator[](i + M*j + M*N*k + M*N*O*l);
	}

	operator mxArray*()
	{
		shouldDestroy = false;
		return array;
	}

	int ndim() const
	{
		return  mxGetNumberOfDimensions(array);
	}

    T min() const
    {
        return *std::min_element(data, data+numel());
    }

    T max() const
    {
        return *std::max_element(data, data+numel());
    }

private:

	bool shouldDestroy;

};

#endif
