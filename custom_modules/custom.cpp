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
#include <math.h>


void Diff_RHS::operator() ( const state_type& X , state_type& dX , const double t )
{
	// EXAMPLE FOR 2 EXTERNAL AND 2 INTERNAL SUBSTRATES
	// This changes the external values
	// dX[0] = -P(00)*(X[3]-X[0]);
	// dX[1] = -P(10)*(X[4]-X[1]);
	// dX[2] = -P(20)*(X[5]-X[2]);

	// This changes internal and pure internal values
	// dX[3] = -dX[0];
	// dX[4] = -dX[1];
	// dX[5] = -dX[2];

	// NOTE: The index can easily lead to segmentation faults when going above the implemented substrate limit
	std::fill(dX.begin(), dX.end(), 0);
	std::cout << "Test\n";
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

	
	cell_defaults.functions.update_phenotype = NULL; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	// Set the update RHS rule for the default cell definition
	// This NEEDS to be done after Intracellular was initialized by the xml file!

	define_cell_parameters();
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void define_cell_parameters( void )
{
	// cell_defaults.phenotype.intracellular = new Diff_Intracellular();
	cell_defaults.phenotype.intracellular->set_parameter_value(00, parameters.doubles("substrate_1_target_intracellular"));
	cell_defaults.phenotype.intracellular->set_parameter_value(01, 0.0);
	cell_defaults.phenotype.intracellular->set_parameter_value(02, parameters.doubles("substrate_1_multiplier_intracellular"));
	cell_defaults.phenotype.intracellular->set_parameter_value(10, parameters.doubles("substrate_2_target_intracellular"));
	cell_defaults.phenotype.intracellular->set_parameter_value(11, 0.0);
	cell_defaults.phenotype.intracellular->set_parameter_value(12, parameters.doubles("substrate_2_multiplier_intracellular"));
	cell_defaults.phenotype.intracellular->set_parameter_value(20, parameters.doubles("killer_production"));
	cell_defaults.phenotype.intracellular->set_parameter_value(21, 0.0);
	cell_defaults.phenotype.intracellular->set_parameter_value(22, 0.0);

	Cell_Definition* differentiation_cell = cell_definitions_by_name["differentiation_cell"];

	differentiation_cell->phenotype.intracellular->set_parameter_value(00, 0.0);
	differentiation_cell->phenotype.intracellular->set_parameter_value(01, 0.0);
	differentiation_cell->phenotype.intracellular->set_parameter_value(02, 0.0);
	differentiation_cell->phenotype.intracellular->set_parameter_value(10, 0.0);
	differentiation_cell->phenotype.intracellular->set_parameter_value(11, 0.0);
	differentiation_cell->phenotype.intracellular->set_parameter_value(12, 0.0);
	differentiation_cell->phenotype.intracellular->set_parameter_value(20, 0.0);
	differentiation_cell->phenotype.intracellular->set_parameter_value(21, 0.0);
	differentiation_cell->phenotype.intracellular->set_parameter_value(22, 0.0);

	differentiation_cell->functions.update_phenotype = diff_phenotype_function;

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

	double frac = PhysiCell::parameters.doubles("fraction_box_height");
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		int N_cells = 0;
		if (k==0) {
		 	N_cells = ceil(parameters.ints("number_of_cells")*frac);
		} else if (k==1) {
			N_cells = ceil(parameters.ints("number_of_cells")*(1-frac));
		}
		for( int n = 0 ; n < N_cells ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			if (k==0) {
				position[1] = Ymin + UniformRandom()*Yrange*frac; 
			} else if (k==1) {
				position[1] = Ymin + Yrange*frac + UniformRandom()*Yrange*(1-frac); 
			}
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
			// pC->phenotype.intracellular = new Diff_Intracellular;
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
	std::string interior_color = "grey"; 
	std::string interior_color_diff_1 = "red";
	std::string interior_color_diff_2 = "blue";
	std::string interior_color_diff_3 = "yellow";
	
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
	
	if (pCell->type_name == "differentiation_cell") {
		output[2] = "white";
		if (fabs(pCell->custom_data["diff"] - PhysiCell::parameters.doubles("diff_enable_1")) < 0.1)
		{
			output[0] = interior_color_diff_1;
			output[2] = interior_color_diff_1;
		}
		if (fabs(pCell->custom_data["diff"] - PhysiCell::parameters.doubles("diff_enable_2")) < 0.1)
		{
			output[0] = interior_color_diff_2;
			output[2] = interior_color_diff_2;
		}
		if (fabs(pCell->custom_data["diff"] - PhysiCell::parameters.doubles("diff_enable_3")) < 0.1)
		{
			output[0] = interior_color_diff_3;
			output[2] = interior_color_diff_3;
		}
	}

	return output;

}


// We want the cell to differentiate with a certain probability when enough substrate is around
void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	return;
}


void diff_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	int voxel_index = microenvironment.nearest_voxel_index(pCell->position);
	std::vector<double> densities = microenvironment.nearest_density_vector(voxel_index);
	
	// only for debugging
	/* bool test = densities[0] > densities[1];
	std::cout << densities << test << "\n";*/


	// If substrate_2 is higher than substrate_1
	// ==> Differentiate in state 2
	if (densities[0] > std::max(densities[1], densities[2])) {
		pCell->custom_data["diff"] = parameters.doubles("diff_enable_1");
	}

	// If substrate_1 is higher than substrate_2
	// ==> Differentiate in state 1
	if (densities[1] > std::max(densities[0], densities[2])) {
		pCell->custom_data["diff"] = parameters.doubles("diff_enable_2");
	}

	// If substrate_1 is higher than substrate_2
	// ==> Differentiate in state 1
	if (densities[2] > std::max(densities[0], densities[1])) {
		pCell->custom_data["diff"] = parameters.doubles("diff_enable_3");
	}

	// Set the death rate according to the density of the killer substrate
	int death_index = phenotype.death.find_death_model_index(100);
	phenotype.death.rates[death_index] = std::min((densities[3] - parameters.doubles("killer_threshold"))/parameters.doubles("killer_modulation"), 0.0);

	return;
}


void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; }


void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 