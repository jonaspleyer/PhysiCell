/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
/*
 * spectrum.h
 *
 *  Created on: 17.01.2020
 *      Author: cfleck
 */

#ifndef CUSTOM_MODULES_SPECTRUM_H_
#define CUSTOM_MODULES_SPECTRUM_H_

#include "domain.h"

// For debugging
#ifndef _HI_
#define _HI_ std::cout << "HERE I AM" << std::endl;
#endif


// TODO: provide a general spectrum class in which the spectrum is discretely or by a function defined
// TODO: find a way to combine the spectra of light sources.

class Spectral_Density
{
protected:
	std::string type;
	std::vector<double> param;
public:
	Spectral_Density() : param(2) {};
	virtual ~Spectral_Density() {};
	virtual void set_param( std::vector<double> _param ) = 0;
	virtual double operator()( double w ) = 0;
	std::string get_type( void ) { return type; }
	std::vector<double> get_param( void ) { return param; }
};

class Delta : public Spectral_Density
{
private:
	double width = 0.05; // band width of the approximate delta distribution in nm
	double norm = 2.0 * width;
public:
	Delta() { type = "Delta"; }
	Delta( std::vector<double> _param ) { type = "Delta"; param = _param; }
	void set_param( std::vector<double> _param ) { param = _param; }
	inline double operator()( double w ) { return ((w >= param[0] - width) && (w <= param[0] + width)) ? 1.0/norm : 0.0; }
};

class Uniform : public Spectral_Density
{
private:
	double norm;
public:
	Uniform() { type = "Uniform"; norm = 0;};
	Uniform( std::vector<double> _param ) { type = "Uniform"; param = _param; norm = param[0] - param[1]; }
	void set_param( std::vector<double> _param ) { param = _param; norm = param[0] - param[1]; }
	inline double operator()( double w ) { return ((w >= param[0]) && (w <= param[1])) ? 1.0/norm : 0.0; }
};

class Gauss : public Spectral_Density
{
private:
	double norm;
public:
	Gauss() { type = "Gauss"; norm = 0; }
	Gauss( std::vector<double> _param ) { type = "Gauss"; param = _param; norm = sqrt(2.0*PhysiCell::PhysiCell_constants::pi)*param[1]; }
	void set_param( std::vector<double> _param ) { param = _param; norm = sqrt(2.0*PhysiCell::PhysiCell_constants::pi)*param[1]; }
	inline double operator()( double w ) { return exp( -pow(w - param[0],2.0)/(2.0*param[1]*param[1]) )/norm; }
};

#endif /* CUSTOM_MODULES_SPECTRUM_H_ */
