// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include "Node.h"
#include "elements/element.h"
#include "FEStorage.h"
#include "FESolver.h"

using namespace nla3d;

bool readCdbFile(const char *filename, FEStorage *storage, FESolver* solver, ElementType elType);
