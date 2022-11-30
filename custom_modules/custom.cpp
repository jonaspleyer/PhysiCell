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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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
#include "./optogenetics/OptoGen.h"
#include <math.h>

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
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
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	static std::vector< std::string > colors(0); 
	
	// start all black 
	
	std::vector<std::string> output = { "black", "black", "black", "black" }; 
	
	// paint by number -- by cell type 
	std::string interior_color = "rgb(34,139,34)"; 
	// the final red color should be rgb(220,20,60) before it turns black when dying
	
	output[0] = interior_color; // set cytoplasm color
	
	if( pCell->phenotype.death.dead == false ) // if live, color nucleus same color 
	{
		output[2] = interior_color; 
		output[3] = interior_color; 
	}
	// apoptotic cells will retain a black nucleus 
	// if necrotic, color the nucleus brown 
	else if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[2] = "rgb(139,69,19)";
		output[3] = "rgb(139,69,19)";
	}

	// Color the inside of the cell with the current concentration of ions inside it.
	double thresh_start = PhysiCell::parameters.doubles("light_ion_thresh_death_start");
	double thresh_end = PhysiCell::parameters.doubles("light_ion_thresh_death_end");
	double conc = pCell->custom_data["light_ion_concentration"];

	// Creates a values between 0 and 1
	if (conc <= thresh_start) {
		double value = std::max(std::min(thresh_start, conc), 0.0) / thresh_start;

		// If value is at thresh, we want to have no more green 
		double r_value = (1 - value) *  34.0 + value * 220.0;
		double g_value = (1 - value) * 139.0 + value *  20.0;
		double b_value = (1 - value) *  34.0 + value *  60.0;

		output[2] = "rgb(" + std::to_string(r_value) + "," + std::to_string(g_value) +"," + std::to_string(b_value) + ")";
	} else if (conc > thresh_start && conc < thresh_end) {
		double value = std::max(std::min(thresh_end, conc), 0.0) / thresh_end;

		double r_value = (1 - value) * 220.0;
		double g_value = (1 - value) *  20.0;
		double b_value = (1 - value) *  60.0;

		output[2] = "rgb(" + std::to_string(r_value) + "," + std::to_string(g_value) +"," + std::to_string(b_value) + ")";
	} else {
		output[2] = "rgb(0,0,0)";
	}

	return output;
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Increase chance of death with more neighbors
	// Count the number of cells in the current voxel
	// int index = pCell->get_current_voxel_index();
	// int n_cells = std::count_if(PhysiCell::all_cells->begin(), PhysiCell::all_cells->end(), [index](Cell* cell){return cell->get_current_voxel_index()==index;});
	std::vector<Cell*> nearby_cells = pCell->nearby_cells();
	int n_cells = std::count_if(nearby_cells.begin(), nearby_cells.end(), [](Cell* cell){return !(cell->phenotype.death.dead);});

	// The maximum number of cells per voxel is given in the settings xml file
	int max_cells_nearby = PhysiCell::parameters.ints("max_cells_nearby");
	double factor = (1.0 * n_cells) / (1.0 * max_cells_nearby);

	// Set the death rate accordingly
	pCell->phenotype.death.rates[0] = PhysiCell::parameters.doubles("unmodified_death_rate") * factor;

	// std::cout << "Cells: " << n_cells << " ==> Death Rate: " << PhysiCell::parameters.doubles("unmodified_death_rate") << " with factor " << factor << "\n";

	return;
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	double C = pCell->custom_data["light_ion_concentration"];
	double T1 = PhysiCell::parameters.doubles("light_ion_thresh_death_start");
	double T2 = PhysiCell::parameters.doubles("light_ion_thresh_death_end");
	double time_scale = PhysiCell::parameters.doubles("light_ion_death_rate") * dt;
	// If the total ion concentration is too high, initiate the death process
	if (C >= T1 && C <= T2) {
		// Draw a random number between 0 and 1 and if the random number is lower than the calculated value then kill the cell
		double rand = PhysiCell::UniformRandom();
		if ((C - T1) / (T2 - T1) * time_scale > rand) {
			pCell->start_death(pCell->phenotype.death.find_death_model_index(100));
		}
	}

	// If the ion concentration exceeds the upper limit, kill the cell immediately
	if (C > T2) {
		pCell->start_death(pCell->phenotype.death.find_death_model_index(100));
	}
	return;
}

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }
