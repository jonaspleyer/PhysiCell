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

#include "./custom.h"
#include "../BioFVM/BioFVM.h"  
#include "../core/PhysiCell_cell_container.h"
#include <math.h>
#include "./tinycolormap/include/tinycolormap.hpp"
using namespace BioFVM;


//#include "rrc_api.h"
//#include "rrc_types.h"
// #include "rrc_utilities.h"
//extern "C" rrc::RRHandle createRRInstance();

void chemotaxis_bias_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	// quickly find O2
	static int activator_index = microenvironment.find_density_index( "activator" );
	// sample O2
	double density = pCell->nearest_density_vector()[activator_index];

	// set direction along activator gradients
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient(activator_index);
	normalize( &( phenotype.motility.migration_bias_direction ) );
	// set speed proportional to activator, scaled by normoxic O2 ( 160 mmHg);
	// with a maximum of 1.2 micron per minute
	const double max_val = parameters.doubles("max_activator_value");
	const double min_val = parameters.doubles("min_activator_value");

	double theta = pow((max_val - density) / (max_val-min_val),2.0);
	phenotype.motility.migration_speed = 0.1*theta;

	if( phenotype.motility.migration_speed > 0.1 ) {
		phenotype.motility.migration_speed = 0.1;
	}
	// the greater the oxygen, the more biased the motion
	//phenotype.motility.migration_bias = theta;

	return;
}


void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

	/*
	   This parses the cell definitions in the XML config file. 
	*/

	initialize_cell_definitions_from_pugixml(); 

	Cell_Definition* motil_cell = cell_definitions_by_name["motil_cell"];

	motil_cell->functions.volume_update_function = standard_volume_update_function;
	motil_cell->functions.update_velocity = standard_update_cell_velocity;

	motil_cell->functions.update_migration_bias = chemotaxis_bias_function;
	motil_cell->functions.update_phenotype = NULL;
	motil_cell->functions.custom_cell_rule = NULL;//slowly_kill_cell_if_further_away;

	motil_cell->functions.add_cell_basement_membrane_interactions = NULL;
	motil_cell->functions.calculate_distance_to_membrane = NULL;


	build_cell_definitions_maps();
	display_cell_definitions( std::cout );

	// Set the update RHS rule for the default cell definition
	// This NEEDS to be done after Intracellular was initialized by the xml file!

	define_cell_parameters();
//	cell_defaults.phenotype.intracellular->setRHS(RHS_definition);

	return; 
}

void define_cell_parameters( void )
{
	cell_defaults.phenotype.intracellular->set_parameter_value(1, parameters.doubles("k1"));
	cell_defaults.phenotype.intracellular->set_parameter_value(2, parameters.doubles("k2"));
	cell_defaults.phenotype.intracellular->set_parameter_value(3, parameters.doubles("k3"));
	cell_defaults.phenotype.intracellular->set_parameter_value(4, parameters.doubles("k4"));
	cell_defaults.phenotype.intracellular->set_parameter_value(5, parameters.doubles("k5"));
	cell_defaults.phenotype.intracellular->set_parameter_value(6, parameters.doubles("k6"));
	cell_defaults.phenotype.intracellular->set_parameter_value(7, parameters.doubles("k7"));

	cell_defaults.phenotype.intracellular->set_parameter_value(8, parameters.doubles("secretion_activator"));
	cell_defaults.phenotype.intracellular->set_parameter_value(9, parameters.doubles("secretion_inhibitor"));

	return;
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = microenvironment.mesh.bounding_box[3];
	double Ymax = microenvironment.mesh.bounding_box[4];
	double Zmax = microenvironment.mesh.bounding_box[5];

	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0;
		Zmax = 0.0;
	}
	
	double Xrange = Xmax - Xmin;
	double Yrange = Ymax - Ymin;
	double Zrange = Zmax - Zmin;

	double X_offset = Xrange/2.0*(parameters.doubles("space_seperation"));
	double Y_offset = Yrange/2.0*(parameters.doubles("space_seperation"));

	double Xrange_reduced = Xrange-2*X_offset;
	double Yrange_reduced = Yrange-2*Y_offset;

	double XYRatio = Xrange_reduced/Yrange_reduced;
	int N_cells = parameters.ints("number_of_cells");
	int N_X = round(sqrt(N_cells*XYRatio));
	int N_Y = round(sqrt(N_cells/XYRatio));

	// create some of each type of cell
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k];
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
		std::cout << "Placing " << N_X << " x " << N_Y << " cells" << std::endl;

		for( int n = 1 ; n <= N_X ; n++ )
			for ( int m = 1 ; m <= N_Y ; m++ ) {
			{
				std::vector<double> position = {0,0,0};
				position[0] = Xmin + X_offset + Xrange_reduced * n / (N_X+1);
				position[1] = Ymin + Y_offset + Yrange_reduced * m / (N_Y+1);
				position[2] = 0;

				pC = create_cell( *pCD );
				pC->assign_position( position );
				pC->phenotype.intracellular->start();
			}
		}
	}
	std::cout << std::endl;

	// load cells from your CSV file (if enabled)
	// load_cells_from_pugixml();

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
    // start with flow cytometry coloring 
	static std::vector< std::string > output( 4 , "rgb(0,0,0)" );

	if ( pCell->phenotype.death.dead ) {
		return output;
	}

	static int motil_cell_id = find_cell_definition( "motil_cell" )->type;
	if ( pCell->type == motil_cell_id ) {
		// Outer padding
		output[1] = "rgb(0,0,0)";
		// Middle padding
		output[0] = "rgb(150,50,5)";
		// Most inner padding
		output[3] = "rgb(0,0,0)";
		// Color of core
		output[2] = "rgb(230,96,55)";
		return output;
	}

	static int default_cell_id = find_cell_definition( "default" )->type;
	int substrate=0;
//	std::vector<double> boundaries = get_plotting_color_boundary_vals(PhysiCell::all_cells, substrate);
	const double max_val = parameters.doubles("max_activator_value");
	const double min_val = parameters.doubles("min_activator_value");
	double dens = pCell->phenotype.intracellular->get_parameter_value(".0");
	double value = 1-(max_val - dens - min_val)/(max_val-min_val);
	const tinycolormap::Color color = tinycolormap::GetColor(value, tinycolormap::ColormapType::Viridis);

	if ( pCell->type == default_cell_id ) {
		output[1] = "rgb(0,0,0)";
		output[0] = "rgb(" + std::to_string(256*color.r()) + "," + std::to_string(256*color.g()) + "," + std::to_string(256*color.b()) + ")";
		output[3] = "rgb(0,0,0)";
		output[2] = "rgb(" + std::to_string(256*color.r()) + "," + std::to_string(256*color.g()) + "," + std::to_string(256*color.b()) + ")";
	}
	return output; 
}
