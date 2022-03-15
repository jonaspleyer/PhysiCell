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

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include "optogenetics.h"

#ifndef _HI_
#define _HI_ std::cout << "HERE I AM" << std::endl;
#endif

using namespace BioFVM; 
using namespace PhysiCell;

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt ); 

extern Light_Control opto;
extern int index_red;
extern int index_blue;
extern int index_green;

extern std::vector<std::ofstream> cell_number_file;
extern std::vector<std::ofstream> density_file;
extern std::vector<std::ofstream> variance_file;


extern Circle domain1;
extern Rectangle domain2;
extern Rectangle domain3;

// any additional cell types (beyond cell_defaults)
extern Cell_Definition hek293_cell_A;
extern Cell_Definition hek293_cell_B;
extern Cell_Definition hek293_cell_C;
extern Cell_Definition hek293_cell_antithetic;

// custom cell phenotype functions could go here 

// setup functions to help us along 
void create_cell_types( void );
void create_HEK293_cell_A( void );
void create_HEK293_cell_B( void );
void create_HEK293_cell_C( void );
void create_HEK293_cell_antithetic( void );
void setup_tissue( void ); 

// set up the BioFVM microenvironment 
void setup_microenvironment( void ); 

// the phenotype update function for HEK293 cells
void update_cell_HEK293_A( Cell* pCell, Phenotype& phenotype, double dt );
void update_cell_HEK293_B( Cell* pCell, Phenotype& phenotype, double dt );
void update_cell_HEK293_C( Cell* pCell, Phenotype& phenotype, double dt );

// the custom cell rule function for HEK293 cells
void cell_rule_HEK293_A( Cell* pCell, Phenotype& phenotype, double dt );
void cell_rule_HEK293_B( Cell* pCell, Phenotype& phenotype, double dt );
void cell_rule_HEK293_C( Cell* pCell, Phenotype& phenotype, double dt );
void cell_rule_HEK293_antithetic( Cell* pCell, Phenotype& phenotype, double dt );

// the arrest function for HEK293 cells
bool arrest_function_HEK293( Cell* pCell, Phenotype& phenotype, double dt );

// custom pathology coloring function 
std::vector<std::string> my_coloring_function_0( Cell* );
std::vector<std::string> my_coloring_function_1( Cell* );

//void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt ){;}
//void custom_function( Cell* pCell, Phenotype& phenotype , double dt ){;}
//void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt ){;}

