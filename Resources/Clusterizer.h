#ifndef CLUSTERIZER_H
#define CLUSTERIZER_H

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <assert.h>
#include <math.h>

/**
 *
 *	This class is based on the paper "Fast Monte Carlo algorithm for site or bond percolation" by Newman1 and Ziff [PRE 64 016706 (2001)].
 *	Given some connectivity map (node i is connected to node j, etc.) and a subset of the nodes that are "occupied," this will decompose
 *	the occupied nodes into clusters of connected occupied nodes. 
 *
 *	Basic Implementation:
 *	Clusters are stored as a tree. The basic tree structure is stored as an integer assigned to each node (PTR[i].ptr). This value has three meanings:
 *	1. If the node is not occumpid, then ptr=EMPTY.
 *	2. If the node is a "root" site, then ptr is negative and equal to the negative size of the cluster (i.e. number of occupied nodes currently in the cluster).
 *	3. If the node is not a root site, then ptr is the index of the node's parent in the tree.
 *
 *	The recursive function findroot(int i) returns the index of the root of the cluster that includes node i. It also performs path compression
 *	for highly efficient code.
 *
 */

class CClusterizer
{
	//! This class isn't really necessary since it just stores an int. However, future implementations may include a Dim dimensional vector pointing from the node to the parent. This also makes sure ptr is initialized to -1.
	class pointer
	{
		public:
		int ptr;
		pointer(){	ptr = -1;	};
	};
	const int EMPTY;

public:
	const int N;	//!<Total number of nodes.
	
	//Data
	std::vector< std::vector<int> > nn;	//!<Neighbor list
	pointer *PTR;						//!<Array of pointers to store the cluster information

	//Constructor/destructor
	CClusterizer(int n);																					//!<Simple constructor
	CClusterizer(int n, std::vector< std::vector<int> > const &nbrs);										//!<Construct and set the neighbor list
	CClusterizer(int n, std::vector< std::vector<int> > const &nbrs, std::vector<bool> const &include);		//!<Construct, set the neighbor list, and run the decomposition
	~CClusterizer();

	void set_nn(std::vector< std::vector<int> > const &nbrs);

private:
	inline int findroot(int i);
	void amalgamate_clusters(int rbig, int rsmall, int sbig, int ssmall); //Combine two clusters.

public:

	void DecomposeClusters(std::vector<bool> const &include);		//!<Decompose the system into distinct clusters.
	void DecomposeClusters();										//!<Decompose the system into distinct clusters, where all nodes are "included." 

	void GetClusters(std::vector< std::vector<int> > &clusters);	//!<Return the clusters in a vector< vector<int> >. DecomposeCluster() must be called prior to this.
	void GetLargestCluster(std::vector<int> &cluster);				//!<Return just the largest cluster in a vector<int>. DecomposeCluster() must be called prior to this. 
	int  GetNumClusters();											//!<Return the number of distinct clusters. DecomposeCluster() must be called prior to this.
};

CClusterizer::CClusterizer(int n)
	:N(n), EMPTY(-n-1)
{
	PTR = new pointer[N];
};

CClusterizer::CClusterizer(int n, std::vector< std::vector<int> > const &nbrs)
	:N(n), EMPTY(-n-1)
{
	PTR = new pointer[N];
	set_nn(nbrs);
};

CClusterizer::CClusterizer(int n, std::vector< std::vector<int> > const &nbrs, std::vector<bool> const &include)
	:N(n), EMPTY(-n-1)
{
	PTR = new pointer[N];
	set_nn(nbrs);
	DecomposeClusters(include);
};

CClusterizer::~CClusterizer()
{
	delete[] PTR;
}

void CClusterizer::set_nn(std::vector< std::vector<int> > const &nbrs)
{
	nn = nbrs;
};

inline int CClusterizer::findroot(int i)
{
	if (PTR[i].ptr<0) return i;
	int parent = PTR[i].ptr;
	PTR[i].ptr = findroot(parent); //findroot(parent) sets PTR[parent] to be
	return PTR[i].ptr;
};

//Combine two clusters
void CClusterizer::amalgamate_clusters(int rbig, int rsmall, int sbig, int ssmall)
{
	//assertions
	assert(PTR[rbig].ptr < 0);		//rbig is a root
	assert(PTR[rsmall].ptr < 0);	//rsmall is a root

	if(sbig != rbig)		assert(PTR[sbig].ptr == rbig);
	if(ssmall != rsmall)	assert(PTR[ssmall].ptr == rsmall);

	//Add the sizes of the clusters
	PTR[rbig].ptr += PTR[rsmall].ptr;

	//Set rsmall to point to rbig
	PTR[rsmall].ptr = rbig;
}

void CClusterizer::DecomposeClusters()
{
	std::vector<bool> include(N, true);
	DecomposeClusters(include);
}

void CClusterizer::DecomposeClusters(std::vector<bool> const &include)
{
	assert((int)include.size() == N);

	int i,j;
	int s1,s2; //site 1 and 2
	int r1,r2; //root 1 and 2

	for (i=0; i<N; ++i) PTR[i].ptr = EMPTY;
	for (i=0; i<N; ++i) 
	{
		if(!include[i]) continue;	//This can be easily changed to make the input be a particular ordering of the particles, see Newman and Ziff, PRE (2001)

		r1 = s1 = i;
		PTR[s1].ptr = -1;
		for (j=0; j<(int)nn[s1].size(); ++j) 
		{
			int rtemp = findroot(s1);
			assert(rtemp == r1);
			s2 = nn[s1][j];
			if (PTR[s2].ptr!=EMPTY) //If site 2 is empty (or not placed yet), do nothing.
			{
				r2 = findroot(s2); //Find the root of site 2
				if (r2!=r1) 
				{
					//PTR[r1].ptr = -size of the cluster
					if (PTR[r1].ptr>PTR[r2].ptr) //if cluster 2 is BIGGER than cluster 1
					{
						amalgamate_clusters(r2, r1, s2, s1);
						r1 = r2;
					} else //if cluster 2 is SMALLER than cluster 1
					{
						amalgamate_clusters(r1, r2, s1, s2);
					}
				}
			}
		}
	}
}

void CClusterizer::GetLargestCluster(std::vector<int> &cluster)
{
	int size = 0;
	int lc_root = -1;
	for(int i=0; i<N; ++i)
		if(PTR[i].ptr != EMPTY)
			if(PTR[i].ptr < size)
			{
				lc_root = i;
				size = PTR[i].ptr;
			}
	size *= -1;
	assert(size>0);

	cluster.clear();
	cluster.reserve(size);
	for(int i=0; i<N; ++i)
		if(findroot(i) == lc_root)
			cluster.push_back(i);
	assert((int)cluster.size() == size);
}
	
void CClusterizer::GetClusters(std::vector< std::vector<int> > &clusters)
{
	clusters.clear();

	std::vector<int> roots(N,-1);
	int root_number = 0;
	int r;
	for(int i=0; i<N; ++i)
		if(PTR[i].ptr != EMPTY)
		{
			if(i==findroot(i))
			{
				roots[i] = root_number;
				++root_number;
			}
		}

	clusters.assign(root_number, std::vector<int>());
	for(int i=0; i<N; ++i)
		if(PTR[i].ptr != EMPTY)
		{
			r = roots[ findroot(i) ];
			assert(r>=0);
			assert(r<root_number);
			clusters[r].push_back(i);
		}
}

int CClusterizer::GetNumClusters()
{
	int num_roots = 0;
	for(int i=0; i<N; ++i)
	{
		if(PTR[i].ptr != EMPTY)
		{
			if(i==findroot(i))
			{
				assert(PTR[i].ptr < 0);
				++num_roots;
			}
		}
	}
	return num_roots;
}


#endif //CLUSTERIZER_H



