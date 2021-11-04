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
using namespace BioFVM;


//#include "rrc_api.h"
//#include "rrc_types.h"
// #include "rrc_utilities.h"
//extern "C" rrc::RRHandle createRRInstance();

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

	double XYRatio = Xrange/Yrange;
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
				position[0] = Xmin + Xrange * n / (N_X+1);
				position[1] = Ymin + Yrange * m / (N_Y+1);
				position[2] = 0;

				pC = create_cell( *pCD );
				pC->assign_position( position );
				pC->phenotype.intracellular->start();
				pC->phenotype.intracellular->setUpdateFunction(update_RHS_custom);
			}
		}
	}
	std::cout << std::endl;

	// load cells from your CSV file (if enabled)
	// load_cells_from_pugixml();

	return; 
}

void update_intracellular()
{
    // BioFVM Indices
    static int N_species = (*all_cells)[0]->phenotype.molecular.internalized_total_substrates.size();

    #pragma omp parallel for
    for( int i=0; i < (*all_cells).size(); i++ )
    {
    	std::ofstream myfile;
		myfile.open ("logs.txt", std::ios_base::app);
    	for ( int j=0; j<N_species; j++ )
    	{
    		std::string substrate_name = microenvironment.density_names[j];
    		int substrate_index = microenvironment.find_density_index(substrate_name);
			if( (*all_cells)[i]->is_out_of_domain == false  )
			{
				double cell_volume = (*all_cells)[i]->phenotype.volume.total; // Cell Volume
				double substrate_density_internal = (*all_cells)[i]->phenotype.molecular.internalized_total_substrates[substrate_index]/cell_volume; // Intracellular Concentrations

				myfile << substrate_density_internal << ",";
				// Update with OdeSolver
				if ( PhysiCell_globals.current_time > 0 )
				{
					(*all_cells)[i]->phenotype.intracellular->set_parameter_value(substrate_name,substrate_density_internal);
				}

				(*all_cells)[i]->phenotype.intracellular->update(); // SBML Simulation
				(*all_cells)[i]->phenotype.intracellular->update_phenotype_parameters((*all_cells)[i]->phenotype); // Phenotype Simulation

				// Internalized Chemical Update After SBML Simulation
				(*all_cells)[i]->phenotype.molecular.internalized_total_substrates[substrate_index] = (*all_cells)[i]->phenotype.intracellular->get_parameter_value(substrate_name) * cell_volume;
			}
        }
    	myfile << std::endl;
    	myfile.close();
    }
}

void update_RHS_custom(const std::vector<double> &X, std::vector<double> &dX, const double t)
{
	std::vector<double> result;
	for ( int i=0; i<X.size(); i++ ) {
		result.push_back(100.0 + i*20 - X[i]);
	}
	dX = result;
	return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
    // start with flow cytometry coloring 
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 

	// color
    // proliferative cell
	output[0] = "rgb(20,20,20)";
	output[2] = "rgb(10,10,10)";

	return output; 
}
