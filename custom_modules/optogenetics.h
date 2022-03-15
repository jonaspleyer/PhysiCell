/*
 * optogenetics.h
 *
 *  Created on: 28.10.2019
 *  Author: Christian Fleck
 *  Institution: ETH Zürich
 *
 *  Project: CyGenTiG
 */

// ToDo: split the files optogenetics.h and optogenetics.cpp

#ifndef CUSTOM_MODULES_OPTOGENETICS_H_
#define CUSTOM_MODULES_OPTOGENETICS_H_

#include "../core/PhysiCell.h"
#include "domain.h"
#include "light.h"
#include "control.h"

using namespace PhysiCell;

// For debugging
#ifndef _HI_
#define _HI_ std::cout << "HERE I AM" << std::endl;
#endif

Control_Lattice create_control_lattice( Control_Lattice control_lattice );
void setup_controll_system( void );

#endif /* CUSTOM_MODULES_OPTOGENETICS_H_ */
