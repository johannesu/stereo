//
// Written by Petter Strandmark 2010.
//
// Contains a function for flushing the output buffer to
// MATLAB
//
#ifndef MEX_UTILS_PETTER
#define MEX_UTILS_PETTER

#include <map>
#include <string>
#include <stdexcept>
#include <vector>
using std::map;
using std::string;
using std::vector;
using std::runtime_error;

#include "mex.h"
#include <sstream>
#define ASSERT(cond) if (!(cond)) { std::stringstream sout; \
                                    sout << "Error (file " << __FILE__ \
                                    	 <<", line " << __LINE__ << "): " \
                                         << #cond; \
                                    throw runtime_error(sout.str()); }

#include "cppmatrix.h"

#if defined (_WIN32)
    #include <windows.h>
#elif defined (__linux__)
    #include <unistd.h>
#endif
//Detect Ctrl+C
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
#else
    extern bool utIsInterruptPending();
#endif


namespace {

	void flush_output()
	{
		//Hack
		mexEvalString("drawnow");
	}



class MexParams
{
public:

	MexParams(int nrhs, const mxArray	*prhs[])
	{
		if (nrhs < 1) {
			return;
		}

		int i = 0;
		if (mxGetClassID(prhs[i])==mxSTRUCT_CLASS) {
			ASSERT(mxGetNumberOfElements(prhs[i]) == 1);
			int nfield = mxGetNumberOfFields(prhs[i]);
			for (int field=0;field<nfield;++field) {
				const char*    fieldname = mxGetFieldNameByNumber(prhs[i],field);
				params[fieldname] = mxGetFieldByNumber(prhs[i],0,field);
			}
			i++;
		}

		for (;i<nrhs-1; i+=2) {
			const mxArray* mxKey = prhs[i];
			const mxArray* mxVal = prhs[i+1];
			char buffer[1024];
			if (mxGetString(mxKey, buffer, 1024)) {
				throw runtime_error("Expected string parameter name");
			}
			params[buffer] = mxVal;
		}
	}
	~MexParams()
	{
	}

	template<typename T> T getdef(const string& key, T def );
	template<typename T> T get(const string& key, T def=T() )
	{
		return getdef<T>(key,def);
	}

private:
	map<string,const mxArray*> params;
};

//
// Default: return a scalar
//
template<typename T> T MexParams::getdef(const string& key, T def)
{
	const mxArray* mxVal = params[key];
	if (!mxVal) {
		//throw runtime_error("No such key: " + key);
		return def;
	}

	matrix<T> my_val(mxVal);

	ASSERT( my_val.numel() == 1 );
	return my_val[0];
}

//
// const mxArray*
//
template<> const mxArray* MexParams::getdef(const string& key, const mxArray* def)
{
	const mxArray* mxVal = params[key];
	if (!mxVal) {
		return def;
	}
	return mxVal;
}

//
// Vector of doubles
//
template<> vector<double> MexParams::getdef(const string& key, vector<double> def)
{
	const mxArray* mxVal = params[key];
	if (!mxVal) {
		//throw runtime_error("No such key: " + key);
		return def;
	}

	ASSERT( mxIsDouble(mxVal) );
	ASSERT( !mxIsComplex(mxVal) );
	ASSERT( !mxIsSparse(mxVal) );
	size_t n = mxGetNumberOfElements(mxVal);

	vector<double> vec(n);
	double* ptr = (double*)mxGetPr(mxVal);
	for (size_t i=0;i<n;++i) {
		vec[i] = ptr[i];
	}

	return vec;
}

//
// String
//
template<> string MexParams::getdef(const string& key, string def)
{
	const mxArray* mxVal = params[key];
	if (!mxVal) {
		//throw runtime_error("No such key: " + key);
		return def;
	}

	char buffer[1024];
	if (mxGetString(mxVal, buffer, 1024)) {
		throw runtime_error("Expected string parameter value");
	}

	return string(buffer);
}

}

#endif
