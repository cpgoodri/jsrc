#ifndef DATABASE_H
#define DATABASE_H

#include "../Resources/std_include.h"
#include "../Resources/Exception.h"
#include "../Resources/Resources.h"
#include "netcdfcpp.h"

namespace LiuJamming
{

class CDatabase
{
public:
	string filename;
	NcFile File;

	CDatabase(string fn);

	void OpenReadOnly();
	void OpenReplace();
	void OpenWrite();

};

CDatabase::CDatabase(string fn)
	: filename(fn)
{};


}

#endif //DATABASE_H
