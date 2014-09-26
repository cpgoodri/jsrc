#ifndef STD_INCLUDE_H
#define STD_INCLUDE_H

#include <Eigen/Core>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <string>
#include <sstream>
#include <stdlib.h>

#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include <iterator>

#include <math.h>
#include <cmath>

#include <sys/stat.h>
#include <stdexcept>
#include <assert.h>


#define SQRT2 1.4142135623730950488016887242096980785696718753769480
#define SQRT3 1.7320508075688772935274463415058723669428052538103806
#define SQRT6 2.4494897427831780981972840747058913919659474806566701

#define PI M_PI

#define DB_STRING_SIZE 128


//#define DONT_USE_SUITESPARSE
//#define DONT_USE_ARPACK
//#define DONT_USE_NETCDF


//using namespace std;
using std::vector;
using std::list;
using std::string;
using std::ostream;
using std::istream;


typedef double dbl;
typedef std::complex<dbl> cdbl;



#define LOGM_0()        {printf("DEBUGGING:         FILE %s,  LINE %d\n", __FILE__, __LINE__); fflush(stdout);}
#define LOGM_1(index)   {printf("DEBUGGING (%03i):   FILE %s,  LINE %d\n", index, __FILE__, __LINE__); fflush(stdout);}
#define LOGM_X(x,A,FUNC,...) FUNC
#define LOGM(...) LOGM_X(,##__VA_ARGS__,LOGM_1(__VA_ARGS__),LOGM_0(__VA_ARGS__))






/*
namespace Constants
{
	const dbl sqrt2 = 1.4142135623730950488016887242096980785696718753769480;
	const dbl sqrt3 = 1.7320508075688772935274463415058723669428052538103806;

	const dbl FCC_DENSITY          = 0.740480489693061041169313498343;  //pi*sqrt(2)/6
	const dbl BCC_DENSITY          = 0.680174761587831693972779346616;  //pi*sqrt(3)/8
	const dbl SIMPLE_CUBIC_DENSITY = 0.523598775598298873077107230547;  //pi/6
};


}
*/
#endif //STD_INCLUDE_H




















