#ifndef SOFTSPOTCOMPUTER
#define SOFTSPOTCOMPUTER

#include "../Resources/std_include.h"
#include "../Resources/Clusterizer.h"

namespace LiuJamming
{

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int Dim>
class CSoftSpotCalculator
{
	typedef Eigen::Matrix<dbl,Dim,1>   dvec;
	typedef Eigen::Matrix<dbl,Dim,Dim> dmat;

	static bool comparator(std::pair<int,dbl> const &i, std::pair<int,dbl> const &j){return (i.second > j.second);};
	static bool comparator2(std::vector<int> const &i, std::vector<int> const &j){return (i.size() > j.size());};

public:

	static void Binarize(int N, dbl *vecs, int Nm, int Np, vector<bool> &softparticles);
	static void Clusterize(vector< vector<int> > const &nbrs, vector<bool> &softparticles, vector< vector<int> > &softspots, int ClusterSizeCutoff=0);

	static void ConvertSSListToArray(int N, vector< vector<int> > const &softspotsList, vector<int> &softspotsArray);
	static void ConvertSSArrayToList(vector<int> const &softspotsArray, vector< vector<int> > &softspotsList);
	static bool CheckSSConsistency(vector<int> const &softspotsArray, vector< vector<int> > const &softspotsList);
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION    //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


/**
 *	For each of the first Nm eigenvectors in vecs, find the Np particles
 *	with the highest participation. softparticles[i] should be true if 
 *	particle i is one of these highly participating particles in any
 *	of the first Nm modes.
 *
 *	@param[in]		N				Number of particles
 *	@param[in]		vecs			Array of eigenvectors, each of length Dim*N. Must be at least Nm eigenvectors.
 *	@param[in]		Nm				Number of lowest modes to consider.
 *	@param[in]		Np				Number of highly participating particles in each mode to consider.
 *	@param[out]		softparticles	True if a particle is "soft".
 */
template <int Dim>
void CSoftSpotCalculator<Dim>::Binarize(int N, dbl *vecs, int Nm, int Np, vector<bool> &softparticles)
{
	//initialize softparticles to false
	softparticles.assign(N,false);

	//Other initializations
	int Nvar = Dim*N;
	vector< std::pair<int,dbl> > pmags(N, std::pair<int,dbl>());
	dvec ptemp;

	//Loop over the first Nm modes
	for(int m=0; m<Nm; ++m)
	{
		//Loop over all the particles
		for(int p=0; p<N; ++p)
		{
			//Extract the polarization vector
			for(int dd=0; dd<Dim; ++dd)	ptemp(dd) = vecs[Nvar*m+Dim*p+dd];
			//calculate the magnitude of the polarization vector (ptemp.norm()).
			pmags[p] = make_pair(p, ptemp.norm());
		}
		//Note, pmags is a vector of pairs. Each pair has a particle number and the polarization of that particle. This way, we can 
		//sort them and know which particles participate the most.

		//Sort pmags using the "comparator" function
		std::sort(pmags.begin(), pmags.end(), comparator); //Sorts from biggest to smallest

		//Loop over the first Np particles and designate them as "soft"
		for(int p=0; p<Np; ++p)
			softparticles[pmags[p].first] = true;
	}
}


/**
 *	Take a neighbor list and decompose a list of soft particles into clusters, or soft spots.
 *	If any cluster has less than ClusterSizeCutoff, it is discarded and the list of softparticles is changed to reflect this.
 *
 *	@param[in]		nbrs				List of neighbors. nbrs.size() should be the total nuber of particles, and equal to softparticles.size()
 *	@param[in,out]	softparticles		List of soft particles (for example, the output of CSoftSpotCalculator<Dim>::Binarize). On output, only particles in soft clusters of size greater than or equal to ClusterSizeCutoff are still considered soft.
 *	@param[out]		softspots			List of soft spots
 *	@param[in]		ClusterSizeCutoff	Minimum size of a soft spot. Default = 0.
 */
template <int Dim>
void CSoftSpotCalculator<Dim>::Clusterize(vector< vector<int> > const &nbrs, vector<bool> &softparticles, vector< vector<int> > &softspots, int ClusterSizeCutoff)
{
	//check
	assert(nbrs.size() == softparticles.size());

	int N = (int)nbrs.size();
	CClusterizer cizer(N, nbrs, softparticles);	//This performs the decomposition.

	//Get the clusters. Output is softspots
	//softspots.size() is number of soft spots (number of clusters)
	//softspots[i] is a vector of the particles in the ith soft spot.
	cizer.GetClusters(softspots);

	//Sort the clusters by size
	std::sort(softspots.begin(), softspots.end(), comparator2);

	//remove entreis where softspots[i].size()<ClusterSizeCutoff
	for(vector< vector<int> >::iterator it=softspots.begin(); it!=softspots.end(); ++it)
	{
		if( (int)it->size() < ClusterSizeCutoff )
		{
			softspots.erase(it, softspots.end());
			break;
		}
	}
	
	//find new soft particles
	if(ClusterSizeCutoff > 0)
	{
		vector<bool> temp(softparticles.size(), false);
		for(vector< vector<int> >::iterator it=softspots.begin(); it!=softspots.end(); ++it)
			for(vector<int>::iterator it2=it->begin(); it2!=it->end(); ++it2)
				temp[(*it2)] = true;

		//check that no new soft particles were created.
		vector<bool>::iterator it1=softparticles.begin(), ittemp=temp.begin();
		for( ; ittemp!=temp.end(); ++it1, ++ittemp) if((*ittemp)) assert( (*it1) );

		//swap "temp" and "softparticles"
		temp.swap(softparticles);
	}
}

template <int Dim>
void CSoftSpotCalculator<Dim>::ConvertSSListToArray(int N, vector< vector<int> > const &softspotsList, vector<int> &softspotsArray)
{
	softspotsArray.assign(N, -1);
	int k;
	for(int ssi=0; ssi<(int)softspotsList.size(); ++ssi)
		for(int j=0; j<(int)softspotsList[ssi].size(); ++j)
		{
			k = softspotsList[ssi][j];
			assert(softspotsArray[k] == -1);
			softspotsArray[k] = ssi;
		}

	if(!CheckSSConsistency(softspotsArray, softspotsList))	assert(false);
}

template <int Dim>
void CSoftSpotCalculator<Dim>::ConvertSSArrayToList(vector<int> const &softspotsArray, vector< vector<int> > &softspotsList)
{
	int ssmax = -1;
	for(vector<int>::const_iterator it=softspotsArray.begin(); it!=softspotsArray.end(); ++it)
	{
		assert((*it) >= -1);
		if((*it)>ssmax)
			ssmax = (*it);
	}
	if(ssmax==-1) return;


	softspotsList.assign(ssmax+1, vector<int>());
	for(int i=0; i<(int)softspotsArray.size(); ++i)
		if(softspotsArray[i] != -1)
			softspotsList[softspotsArray[i]].push_back(i);
	
	if(!CheckSSConsistency(softspotsArray, softspotsList))	assert(false);
}

template <int Dim>
bool CSoftSpotCalculator<Dim>::CheckSSConsistency(vector<int> const &softspotsArray, vector< vector<int> > const &softspotsList)
{
	int N = (int)softspotsArray.size();
	int Nsoftspots = (int)softspotsList.size();

	int temp;
	int Nsoftparticles1 = 0;
	for(int i=0; i<(int)softspotsArray.size(); ++i)
	{
		temp = softspotsArray[i];
		if(temp >= Nsoftspots)	return false;
		if(temp <  -1)			return false;		
		if(temp != -1) ++Nsoftparticles1;
	}

	int Nsoftparticles2 = 0;
	for(int ssi=0; ssi<Nsoftspots; ++ssi)
	{
		for(int j=0; j<(int)softspotsList[ssi].size(); ++j)
		{
			temp = softspotsList[ssi][j];
			if(temp >= N)	return false;
			if(temp <  0)	return false;
			if(softspotsArray[temp] != ssi) return false;
		}
		Nsoftparticles2 += (int)softspotsList[ssi].size();
	}

	if(Nsoftparticles1 != Nsoftparticles2)
		return false;

	return true;
}


};






















#endif

