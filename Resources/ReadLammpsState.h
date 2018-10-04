#ifndef READ_LAMMPS_STATE
#define READ_LAMMPS_STATE

#include "Resources/Resources.h"

class CLammpsStateData
{
public:
	int NAtoms;
	int NAtomTypes;
	dbl xlo, xhi, ylo, yhi, zlo, zhi;
	dbl xy, xz, yz;
	Eigen::VectorXd masses; //Length NAtomTypes
	Eigen::VectorXd pos;	//Length 3*NAtoms
	Eigen::VectorXi type;	//Length NAtoms
	Eigen::VectorXi wrap;	//Length 3*NAtoms
	Eigen::VectorXd vel;	//Length 3*NAtoms

	CLammpsStateData(){};
	CLammpsStateData(char *fn)
	{
		ReadLammpsState(fn);
	};
	
	template<int DIM> void SetTransformationMatrix(Eigen::Matrix<dbl,DIM,DIM> &T);
	int ReadLammpsState(char *fn);
};

template<int DIM>
void CLammpsStateData::SetTransformationMatrix(Eigen::Matrix<dbl,DIM,DIM> &T)
{
	T = Eigen::Matrix<dbl,DIM,DIM>::Zero();
	T(0,0) = xhi-xlo;
	T(1,1) = yhi-ylo;
	T(0,1) = xy;
	if(DIM==3)
	{
		T(DIM-1,DIM-1) = zhi-zlo;
		T(0,DIM-1) = xz;
		T(1,DIM-1) = yz;   //Check that these off-diagonal terms are correct!!!
	}
}

int CLammpsStateData::ReadLammpsState(char *fn)
{
//	AssertThatFileExists(fn);
	if(!FileExists(fn))
	{
		printf("Warning: file %s does not exist\n", fn);
		return 0;
	}
	string line;
	ifstream file(fn);

	assert(file.is_open());

	std::getline(file,line); //throw away the first line
	if(line.substr(0,32) != string("LAMMPS data file via write_data,"))
	{
		printf("Warning: line 1 = %s\n", line.c_str());
		return 0;
	}
	std::getline(file,line); //throw away the second line

	//process the third line
	file >> NAtoms;
	std::getline(file,line); //read the rest of the line

	//process the fourth line
	file >> NAtomTypes;
	std::getline(file,line); //read the rest of the line

	std::getline(file,line); //throw away the fifth line

	//process the 3 box size lines
	string tmp_str;
	file >> xlo >> xhi >> tmp_str >> tmp_str;
	file >> ylo >> yhi >> tmp_str >> tmp_str;
	file >> zlo >> zhi >> tmp_str >> tmp_str;
	xy = xz = yz = 0.;
	
	//Throw away the rest of the current line
	std::getline(file,line); 


	//initialize vectors
	masses	= Eigen::VectorXd::Zero(NAtomTypes);
	pos		= Eigen::VectorXd::Zero(3*NAtoms);
	type	= Eigen::VectorXi::Zero(NAtoms);
	wrap	= Eigen::VectorXi::Zero(3*NAtoms);
	vel		= Eigen::VectorXd::Zero(3*NAtoms);

	//Throw away the next 3 lines
	std::getline(file,line); 
	std::getline(file,line); 
	if(line != string("Masses"))
	{
		printf("Warning: line 10 = %s\n", line.c_str());
		return 0;
	}
	std::getline(file,line); 

	//read the masses
	int i_temp;
	for(int i=0; i<NAtomTypes; ++i)
	{
		file >> i_temp >> masses[i];
		assert(i_temp == i+1);
	}
	//Throw away the rest of the current line
	std::getline(file,line); 

	//Throw away the next 3 lines
	std::getline(file,line); 
	std::getline(file,line); 
	if(line != string("Atoms"))
	{
		printf("Warning: line 15 = %s\n", line.c_str());
		return 0;
	}
	std::getline(file,line); 

	//read the positions, types and wraps
	int index, ptype, wx, wy, wz;
	dbl x, y, z;
	for(int i=0; i<NAtoms; ++i)
	{
		if(!(file >> index >> ptype >> x >> y >> z >> wx >> wy >> wz))
		{
			printf("Warning, reading data failed\n");
			return 0;
		}
		index--;
		assert(index>=0);
		assert(index<NAtoms);
		assert(type[index] == 0);

		type[index] = ptype;
		pos[3*index+0] = x;
		pos[3*index+1] = y;
		pos[3*index+2] = z;
		wrap[3*index+0] = wx;
		wrap[3*index+1] = wy;
		wrap[3*index+2] = wz;
	}

	//Throw away the next 3 lines
	std::getline(file,line); 
	std::getline(file,line); 
	std::getline(file,line); 

	//read the velocities
	dbl vx, vy, vz;
	for(int i=0; i<NAtoms; ++i)
	{
		if(!(file >> index >> vx >> vy >> vz))
		{
			printf("Warning, reading data failed\n");
			return 0;
		}
		index--;
		assert(index>=0);
		assert(index<NAtoms);

		vel[3*index+0] = vx;
		vel[3*index+1] = vy;
		vel[3*index+2] = vz;
	}

	return 1;
}










#endif //READ_LAMMPS_STATE


