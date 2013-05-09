#ifndef POTENTIALS_H
#define POTENTIALS_H

#include "../Resources/std_include.h"
#include "BasePotential.h"
#include "HarmonicPotential.h"

namespace LiuJamming
{

static std::map<string,CPotential*> CreatePotentialMap()
{
	std::map<string,CPotential*> m;
	CPotential *p;
	
	//For each potential...
	p = new CHarmonicPotential(); m[CHarmonicPotential::GetName()] = p;

	return m;
}

std::map<string,CPotential*> CPotential::PotentialTypes = CreatePotentialMap();


};

#endif //POTENTIALS_H




