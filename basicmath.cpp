#include <immintrin.h>

#include "basicmath.h"

namespace delsquared {

template void cXc(dsComplex<double>& cout, const dsComplex<double> & cin1, const dsComplex<double> & cin2, const unsigned int vlength);
template void cXc(dsComplex<double>& cout, const dsComplex<float>  & cin1, const dsComplex<double> & cin2, const unsigned int vlength);
template void cXc(dsComplex<double>& cout, const dsComplex<double> & cin1, const dsComplex<float>  & cin2, const unsigned int vlength);
template void cXc(dsComplex<double>& cout, const dsComplex<float>  & cin1, const dsComplex<float>  & cin2, const unsigned int vlength);

template void cXc(dsComplex<float>&  cout, const dsComplex<double> & cin1, const dsComplex<double> & cin2, const unsigned int vlength);
template void cXc(dsComplex<float>&  cout, const dsComplex<float>  & cin1, const dsComplex<double> & cin2, const unsigned int vlength);
template void cXc(dsComplex<float>&  cout, const dsComplex<double> & cin1, const dsComplex<float>  & cin2, const unsigned int vlength);
template void cXc(dsComplex<float>&  cout, const dsComplex<float>  & cin1, const dsComplex<float>  & cin2, const unsigned int vlength);

//template void rXc(dsComplex<double>& cout, const dsReal<double> & rin1, const dsComplex<double> & cin2, const unsigned int vlength);
//template void rXc(dsComplex<double>& cout, const dsReal<float>  & rin1, const dsComplex<double> & cin2, const unsigned int vlength);
//template void rXc(dsComplex<double>& cout, const dsReal<double> & rin1, const dsComplex<float>  & cin2, const unsigned int vlength);
//template void rXc(dsComplex<double>& cout, const dsReal<float>  & rin1, const dsComplex<float>  & cin2, const unsigned int vlength);
//
//template void rXc(dsComplex<float>&  cout, const dsReal<float>  & rin1, const dsComplex<float>  & cin2, const unsigned int vlength);
//template void rXc(dsComplex<float>&  cout, const dsReal<double> & rin1, const dsComplex<float>  & cin2, const unsigned int vlength);
//template void rXc(dsComplex<float>&  cout, const dsReal<float>  & rin1, const dsComplex<double> & cin2, const unsigned int vlength);
//template void rXc(dsComplex<float>&  cout, const dsReal<double> & rin1, const dsComplex<double> & cin2, const unsigned int vlength);

template void su3mXsu3v( su3v<double>& vout, const su3m<double>& min, const su3v<double>& vin, const unsigned int vlength, su3v<double>* const t);
template void su3mXsu3v( su3v<double>& vout, const su3m<float>&  min, const su3v<double>& vin, const unsigned int vlength, su3v<double>* const t);
template void su3mXsu3v( su3v<double>& vout, const su3m<double>& min, const su3v<float>&  vin, const unsigned int vlength, su3v<double>* const t);
template void su3mXsu3v( su3v<double>& vout, const su3m<float>&  min, const su3v<float>&  vin, const unsigned int vlength, su3v<double>* const t);

template void su3mXsu3v( su3v<float>& vout, const su3m<float>&  min, const su3v<float>&  vin, const unsigned int vlength, su3v<float>* const t);
template void su3mXsu3v( su3v<float>& vout, const su3m<double>& min, const su3v<float>&  vin, const unsigned int vlength, su3v<float>* const t);
template void su3mXsu3v( su3v<float>& vout, const su3m<float>&  min, const su3v<double>& vin, const unsigned int vlength, su3v<float>* const t);
template void su3mXsu3v( su3v<float>& vout, const su3m<double>& min, const su3v<double>& vin, const unsigned int vlength, su3v<float>* const t);

template void su3mXsu3v_direct( su3v<double>& vout, const su3m<double>& min, const su3v<double>& vin, const unsigned int vlength, su3v<double>* const t);
template void su3mXsu3v_direct( su3v<double>& vout, const su3m<float>&  min, const su3v<double>& vin, const unsigned int vlength, su3v<double>* const t);
template void su3mXsu3v_direct( su3v<double>& vout, const su3m<double>& min, const su3v<float>&  vin, const unsigned int vlength, su3v<double>* const t);
template void su3mXsu3v_direct( su3v<double>& vout, const su3m<float>&  min, const su3v<float>&  vin, const unsigned int vlength, su3v<double>* const t);

template void su3mXsu3v_direct( su3v<float>& vout, const su3m<float>&  min, const su3v<float>&  vin, const unsigned int vlength, su3v<float>* const t);
template void su3mXsu3v_direct( su3v<float>& vout, const su3m<double>& min, const su3v<float>&  vin, const unsigned int vlength, su3v<float>* const t);
template void su3mXsu3v_direct( su3v<float>& vout, const su3m<float>&  min, const su3v<double>& vin, const unsigned int vlength, su3v<float>* const t);
template void su3mXsu3v_direct( su3v<float>& vout, const su3m<double>& min, const su3v<double>& vin, const unsigned int vlength, su3v<float>* const t);

//template void rXc(dsComplex<__m128> &cout, dsReal<__m128> & rin1, dsComplex<__m128> & cin2);

}

