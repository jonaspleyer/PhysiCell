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
#include "../core/PhysiCell.h"

// declare cell definitions here 
/*
 * TODO: movel all cell definitions to the XML file
 */

Cell_Definition hek293_cell_A;
Cell_Definition hek293_cell_B;
Cell_Definition hek293_cell_C;
Cell_Definition hek293_cell_antithetic;
//Cell_Definition motile_cell;

// needed for cell internal states
// maybe put this somewhere else?
int ri, bi, gi, opt;

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
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
//
//	cell_defaults.functions.update_phenotype = phenotype_function;
//	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = standard_elastic_contact_function;
	
	// add custom data here, if any 
	create_HEK293_cell_A();
	create_HEK293_cell_B();
	create_HEK293_cell_C();
	create_HEK293_cell_antithetic();

	build_cell_definitions_maps();
	display_cell_definitions( std::cout );

	return; 
}

void create_HEK293_cell_A( void )
{
/***************************************************************
 *
 * Define the HEK293 cells. Base it on the cell_defaults.
 * I put the several rates based on discussions with Jamie Davies; they are NOT based on real data!
 * This could be moved to the XML file under <cell_definitions>
 *
 ***************************************************************/

	hek293_cell_A = cell_defaults;
	hek293_cell_A.type = 1;
	hek293_cell_A.name = "hek293_cell_A";

// volume
	hek293_cell_A.phenotype.volume.total = parameters.doubles("HEK293_cell_total_volume");

// define them as motile
	hek293_cell_A.phenotype.motility.is_motile = true;
	hek293_cell_A.phenotype.motility.migration_speed = parameters.doubles("HEK293_cell_migration_speed");
	hek293_cell_A.phenotype.motility.migration_bias = 0.0;// completely random

// modify the adhesion
	hek293_cell_A.phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("HEK293_cell_cell_adhesion_strength");
	hek293_cell_A.phenotype.mechanics.relative_maximum_adhesion_distance = parameters.doubles("HEK293_cell_relative_maximum_adhesion_distance");
	hek293_cell_A.phenotype.mechanics.cell_cell_repulsion_strength = parameters.doubles("HEK293_cell_cell_repulsion_strength");


//	std::cout << hek293_cell_A.phenotype.mechanics.cell_cell_adhesion_strength << std::endl;
//	hek293_cell_A.phenotype.mechanics.set_relative_equilibrium_distance( 1.0 );
//	std::cout << hek293_cell_A.phenotype.mechanics.cell_cell_adhesion_strength << std::endl;

// Set the cell cycle model
// Flow Cytometry separated cycle model

	hek293_cell_A.functions.cycle_model = flow_cytometry_separated_cycle_model;
	hek293_cell_A.phenotype.sync_to_functions( hek293_cell_A.functions );

	int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );
	int G2_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G2_phase );
	int M_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::M_phase );

///*
//   std::cout << "G0G1_index = " << G0G1_index << std::endl;
//	 std::cout << "S_index = " << S_index << std::endl;
//	 std::cout << "G2_index = " << G2_index << std::endl;
//   std::cout << "M_index = " << M_index << std::endl;
//*/
// Adjust the apoptosis rate
// Attention: there are two rates: one rate (phenotype.death.rates[Index]) which gives the rate by which the cells
// are to be triggered dead, and one rate (phenotype.death.models[Index]->transition_rate(0,1))
// which gives the rate the individual cell goes through the cycle from alife to dead
	int apoptosis_index = hek293_cell_A.phenotype.death.find_death_model_index("apoptosis");

// This slows the shrinking rate down after death induction
//	// There is probably not much known about these detailed rates -> check with Jamie
	hek293_cell_A.phenotype.death.parameters[apoptosis_index].unlysed_fluid_change_rate *= 0.3;
	hek293_cell_A.phenotype.death.parameters[apoptosis_index].lysed_fluid_change_rate *= 0.3;
	hek293_cell_A.phenotype.death.parameters[apoptosis_index].nuclear_biomass_change_rate *= 0.3;
	hek293_cell_A.phenotype.death.parameters[apoptosis_index].cytoplasmic_biomass_change_rate *= 0.3;

// Set apoptosis to zero
//	hek293_cell_A.phenotype.death.rates[apoptosis_model_index] = 0.000002;

// with this rate the cells are removed from the system
	hek293_cell_A.phenotype.death.models[apoptosis_index]->transition_rate(0,1) = 1.0/(18.0*60.0);
//	std::cout << hek293_cell_A.phenotype.death.rates[apoptosis_index] << std::endl;

//  Attention: this changes the overall apoptosis rate (no induced death)
//	hek293_cell_A.phenotype.death.rates[apoptosis_index] = 1.0/(12.0*60.0);
//	std::cout << "Apopotosis transition rate = " << hek293_cell_A.phenotype.death.models[apoptosis_index]->transition_rate(0,1) << std::endl;
//	std::cout << "Current death model = " << hek293_cell_A.phenotype.death.current_model().name << std::endl;

//	 Alter the transition rate from G0G1 state to S state
	hek293_cell_A.phenotype.cycle.data.transition_rate(G0G1_index,S_index) *= 1.2;
	hek293_cell_A.phenotype.cycle.data.transition_rate(S_index, G2_index) *= 2;
	hek293_cell_A.phenotype.cycle.data.transition_rate(G2_index, M_index) *= 2.0;
	hek293_cell_A.phenotype.cycle.data.transition_rate(M_index, G0G1_index) *= 1;
/*
//	Choose a simpler cell cycle model
hek293_cell_A.functions.cycle_model = flow_cytometry_cycle_model;
hek293_cell_A.phenotype.sync_to_functions( hek293_cell_A.functions );

int G0G1_index = flow_cytometry_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );
int G2M_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G2M_phase );

// Alter the transition rate from G0G1 state to S state
hek293_cell_A.phenotype.cycle.data.transition_rate(G0G1_index,S_index) *= 2;
hek293_cell_A.phenotype.cycle.data.transition_rate(S_index, G2M_index) *= 3;
hek293_cell_A.phenotype.cycle.data.transition_rate(G2M_index, G0G1_index) *= 2;
*/
	// provide a phase arrest function
	hek293_cell_A.phenotype.cycle.pCycle_Model->
		phase_link( G0G1_index, S_index ).arrest_function = arrest_function_HEK293;
	hek293_cell_A.phenotype.cycle.pCycle_Model->
			phase_link( G2_index, M_index ).arrest_function = arrest_function_HEK293;

// add some variables for optogenetics
// these just mimic internal state due to optogenetic treatments
	ri = hek293_cell_A.custom_data.add_variable("red_light",0.0);
	bi = hek293_cell_A.custom_data.add_variable("blue_light",0.0);
	gi = hek293_cell_A.custom_data.add_variable("green_light",0.0);
	opt = hek293_cell_A.custom_data.add_variable("opto",1.0);

// provide cell custom function
	hek293_cell_A.functions.custom_cell_rule = cell_rule_HEK293_A;

// provide the phenotype update function
	hek293_cell_A.functions.update_phenotype = update_cell_HEK293_A;

/*******************************************************************************

Not needed any longer, left here for inspiration (CF)

// Now, let's define another cell type.
// It's best to just copy the default and modify it.

// make this cell type randomly motile, less adhesive, greater survival,
// and less proliferative

motile_cell = cell_defaults;
motile_cell.type = 1;
motile_cell.name = "motile tumor cell";

// make sure the new cell type has its own reference phenotype

motile_cell.parameters.pReference_live_phenotype = &( motile_cell.phenotype );

// enable random motility
motile_cell.phenotype.motility.is_motile = true;
motile_cell.phenotype.motility.persistence_time = parameters.doubles( "motile_cell_persistence_time" ); // 15.0;
motile_cell.phenotype.motility.migration_speed = parameters.doubles( "motile_cell_migration_speed" ); // 0.25 micron/minute
motile_cell.phenotype.motility.migration_bias = 0.0;// completely random

// Set cell-cell adhesion to 5% of other cells
motile_cell.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "motile_cell_relative_adhesion" ); // 0.05;

// Set apoptosis to zero
motile_cell.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "motile_cell_apoptosis_rate" ); // 0.0;

// Set proliferation to 10% of other cells.
// Alter the transition rate from G0G1 state to S state
motile_cell.phenotype.cycle.data.transition_rate(G0G1_index,S_index) *=
	parameters.doubles( "motile_cell_relative_cycle_entry_rate" ); // 0.1;

********************************************************************************/

	return;
}

void create_HEK293_cell_B( void )
{
/***************************************************************
 *
 * This defines HEK293 cells that do not cluster
 *
 ***************************************************************/
	hek293_cell_B = hek293_cell_A;
	hek293_cell_B.type = 2;
	hek293_cell_B.name = "HEK293_cell_B";

//	hek293_cell_B.parameters.o2_proliferation_threshold = 8;
//	hek293_cell_B.parameters.o2_necrosis_threshold = 8;

	// modify the adhesion
//	hek293_cell_B.phenotype.mechanics.cell_cell_adhesion_strength = 0.2*parameters.doubles("HEK293_cell_cell_adhesion_strength");
//	hek293_cell_B.phenotype.mechanics.relative_maximum_adhesion_distance = parameters.doubles("HEK293_cell_relative_maximum_adhesion_distance");

	// provide cell custom function
	hek293_cell_B.functions.custom_cell_rule = cell_rule_HEK293_B;

	// provide the phenotype update function
	hek293_cell_B.functions.update_phenotype = update_cell_HEK293_B;
	return;
}

void create_HEK293_cell_C( void )
{
/***************************************************************
 *
 * This defines HEK293 cells that do not cluster
 *
 ***************************************************************/
	hek293_cell_C = hek293_cell_A;
	hek293_cell_C.type = 3;
	hek293_cell_C.name = "HEK293_cell_C";

//	hek293_cell_C.parameters.o2_proliferation_threshold = 10;
//	hek293_cell_C.parameters.o2_necrosis_threshold = 10;

	// modify the adhesion
//	hek293_cell_C.phenotype.mechanics.cell_cell_adhesion_strength = 0.2*parameters.doubles("HEK293_cell_cell_adhesion_strength");
//	hek293_cell_C.phenotype.mechanics.relative_maximum_adhesion_distance = parameters.doubles("HEK293_cell_relative_maximum_adhesion_distance");

	// provide cell custom function
	hek293_cell_C.functions.custom_cell_rule = cell_rule_HEK293_C;

	// provide the phenotype update function
	hek293_cell_C.functions.update_phenotype = update_cell_HEK293_C	;
	return;
}

void create_HEK293_cell_antithetic( void )
{
/***************************************************************
 *
 * This defines HEK293 cells including the production of a Z1 and Z2 molecule for Mustafa's antithetic controller.
 * The phenotype update function is altered accordingly and some of the mechanical properties are also adjusted (no clustering).
 *
 ***************************************************************/
	hek293_cell_antithetic = hek293_cell_A;
	hek293_cell_antithetic.type = 3;
	hek293_cell_antithetic.name = "HEK293_cell_antithetic";

	// add the internal variables Z1 and Z2
	// these just mimic internal state due to optogenetic treatments
	hek293_cell_antithetic.custom_data.add_variable("Z1",500.0);
	hek293_cell_antithetic.custom_data.add_variable("Z2",0.0);
	hek293_cell_antithetic.custom_data.add_variable("mu_Z1",parameters.doubles("HEK293_cell_antithetic_mu_Z1"));
	hek293_cell_antithetic.custom_data.add_variable("theta_Z2",parameters.doubles("HEK293_cell_antithetic_theta_Z2"));
	hek293_cell_antithetic.custom_data.add_variable("eta_annihilation",parameters.doubles("HEK293_cell_antithetic_eta_annihilation"));

	// provide cell custom function
	// this appears not to be the right place for it, but for now I put the intracellular functions there
	hek293_cell_antithetic.functions.custom_cell_rule = cell_rule_HEK293_antithetic;

	return;
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
/* now this is in XML 	
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;

//	pC = create_cell();
//	pC->assign_position( 0.0, 0.0, 0.0 );
//
//	pC = create_cell();
//	pC->assign_position( -100, 0, 0.0 );
//
//	pC = create_cell();
//	pC->assign_position( 0, 100, 0.0 );
	
// now create a motile cell
	
//	pC = create_cell( motile_cell );
//	pC->assign_position( 15.0, -18.0, 0.0 );
    
    double sc = 1.02; // reduces the initial box
    double x_min = default_microenvironment_options.X_range[0]/sc;
    double x_max = default_microenvironment_options.X_range[1]/sc;
    double y_min = default_microenvironment_options.X_range[0]/sc;
    double y_max = default_microenvironment_options.Y_range[1]/sc;
    int nc = parameters.ints("number_of_cells");

    std::cout << std::endl;
    std::cout << "Creating N = " << nc  << " " << hek293_cell_A.name << " cells" << std::endl;
    
    for (int n = 0; n < nc; n++)
    {
        pC = create_cell( hek293_cell_A );
//    	pC = create_cell();
//        double pos_x = x_min + UniformRandom()*(x_max - x_min);
//        double pos_y = y_min + UniformRandom()*(y_max - y_min);
        double r =500*sqrt(UniformRandom());
        double ang = UniformRandom()*2*3.1415;
        double pos_x = r*cos(ang);
        double pos_y = r*sin(ang);
////        std::cout << std::endl;
        std::cout << "Creating " << pC->type_name << " cell at x=" << pos_x << ", y=" << pos_y << std::endl;
        pC->assign_position( pos_x, pos_y, 0.0 );
        // make the cell sensitive to light
        pC->custom_data.variables[opt].value = 1;
//       if( UniformRandom() <= 0.99 )
//       {
//      	pC->custom_data.variables[opt].value = 1;
//       }
    }

    std::cout << std::endl;
    std::cout << "Creating N = " << nc  << " " << hek293_cell_B.name << " cells" << std::endl;

    for (int n = 0; n < nc; n++)
    {
        pC = create_cell( hek293_cell_B );
 //    	pC = create_cell();
//      double pos_x = x_min + UniformRandom()*(x_max - x_min);
//      double pos_y = y_min + UniformRandom()*(y_max - y_min);
        double r = 500*sqrt(UniformRandom());
        double ang = UniformRandom()*2*3.1415;
        double pos_x = r*cos(ang);
        double pos_y = r*sin(ang);
//      std::cout << std::endl;
        std::cout << "Creating " << pC->type_name << " cell at x=" << pos_x << ", y=" << pos_y << std::endl;
        pC->assign_position( pos_x, pos_y, 0.0 );
        // make the cell sensitive to light
        pC->custom_data.variables[opt].value = 1;
//      if( UniformRandom() <= 0.99 )
//      {
//   	  	pC->custom_data.variables[opt].value = 1;
//      }
    }

    std::cout << std::endl;
    std::cout << "Creating N = " << nc  << " " << hek293_cell_C.name << " cells" << std::endl;

    for (int n = 0; n < 0*nc; n++)
    {
    	pC = create_cell( hek293_cell_C );
 //    	pC = create_cell();
 //     double pos_x = x_min + UniformRandom()*(x_max - x_min);
 //     double pos_y = y_min + UniformRandom()*(y_max - y_min);
    	double r = 650*sqrt(UniformRandom());
    	double ang = UniformRandom()*2*3.1415;
    	double pos_x = r*cos(ang);
    	double pos_y = r*sin(ang);
//      std::cout << std::endl;
    	std::cout << "Creating " << pC->type_name << " cell at x=" << pos_x << ", y=" << pos_y << std::endl;
    	pC->assign_position( pos_x, pos_y, 0.0 );
    	// make the cell sensitive to light
    	pC->custom_data.variables[opt].value = 1;
//     if( UniformRandom() <= 0.99 )
//     {
//     		pC->custom_data.variables[opt].value = 1;
//     }
    }
    std::cout << std::endl;
	return; 
}


std::vector<std::string> my_coloring_function_0( Cell* pCell )
{
	// output: vector of four strings:
	// the cytoplasm fill color, cytoplasm outline color,
	// nuclear color, and nuclear outline color

	// start with flow cytometry coloring

    std::vector<std::string> output = false_cell_coloring_cytometry(pCell);
//	std::vector<std::string> output (4,"black");

//	if( pCell->phenotype.death.dead == false && pCell->type == 1 )
//	{
//		 output[0] = "black";
//		 output[2] = "black";
//	}

//    if( pCell->phenotype.death.dead == false && pCell->phenotype.motility.is_motile == true )
//    {
//        output[0] = "black";
//        output[1] = "black";
//        output[2] = "yellow";
//        output[3] = "yellow";
//    }
//    else if( pCell->phenotype.death.dead == false && pCell->phenotype.motility.is_motile == false )
//    {
//    	output[0] = "gray";
//    	output[1] = "gray";
//    	output[2] = "gray";
//    	output[3] = "gray";
//    }
//    else if ( pCell->phenotype.death.dead == true )
//    {
//    	output[0] = "red";
//    	output[1] = "red";
//    	output[2] = "red";
//    	output[3] = "red";
//    }

    return output;
}

//Slightly changed paint_by_number_cell_coloring function
std::vector<std::string> my_coloring_function_1( Cell* pCell )
{
	static std::vector< std::string > colors(0);
	static bool setup_done = false;
	if( setup_done == false )
	{
		colors.push_back( "grey" ); // default color will be grey

		colors.push_back( "blue" );
		colors.push_back( "red" );
		colors.push_back( "yellow" );
		colors.push_back( "green" );

		colors.push_back( "magenta" );
		colors.push_back( "orange" );
		colors.push_back( "lime" );
		colors.push_back( "cyan" );

		colors.push_back( "hotpink" );
		colors.push_back( "peachpuff" );
		colors.push_back( "darkseagreen" );
		colors.push_back( "lightskyblue" );

		setup_done = true;
	}

	// start all black

	std::vector<std::string> output = { "black", "black", "black", "black" };

	// paint by number -- by cell type

	std::string interior_color = "white";
	if( pCell->type < 13 )
	{ interior_color = colors[ pCell->type ]; }

	output[0] = interior_color; // set cytoplasm color

	if( pCell->phenotype.death.dead == false ) // if live, color nucleus same color
	{
		output[2] = interior_color;
		output[3] = interior_color;
	}
	else
	{
		// apoptotic cells will retain a black nucleus
		// if necrotic, color the nucleus brown
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling ||
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed ||
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
		{
			output[2] = "rgb(139,69,19)";
			output[3] = "rgb(139,69,19)";
		}
	}
	return output;
}

/*
 * Call hierarchy:
 * Cell_Container::update_all_cells -> Cell::advance_bundled_phenotype_functions -> Cell::Functions::update_cell_HEK293
 * Cell_Container::update_all_cells -> Cell::Functions::cell_rule_HEK293
 *
 */
/*
 * This function is called from the function "Cell::advance_bundled_phenotype_functions", where all the other
 * phenotype functions are called.
 *
 */

inline void update_cell_HEK293_A( Cell* pCell, Phenotype& phenotype, double dt )
{
	// first call the standard update function
	update_cell_and_death_parameters_O2_based( pCell, phenotype, dt );

	// get the current position
	Point P{pCell->position[0], pCell->position[1]};
//	std::cout << P.tuple() << std::endl;

	// current time
	double time = PhysiCell_globals.current_time;

	// Neighbors?
//	std::vector<Cell*> nearby = pCell->cells_in_my_container();

//	then do HEK293 cell specific things

/* **************
 * Optogenetics *
 ****************/
	// the phenotype governs all other properties like mechanics or cell cycle
	// set the internal state according to the light intensities
	if( pCell->custom_data.variables[opt].value > 0 )
	{
		pCell->custom_data.variables[bi].value = opto.get_intensity(P, time, index_blue);
		pCell->custom_data.variables[ri].value = opto.get_intensity(P, time, index_red);
		pCell->custom_data.variables[gi].value = opto.get_intensity(P, time, index_green);
	}
	if( pCell->custom_data.variables[opt].value > 0 && pCell->custom_data.variables[bi].value > 0
			&& !pCell->phenotype.death.dead )
	{
		if( UniformRandom() <= 0.9*pCell->custom_data.variables[bi].value )
		{
//			std::cout << "Intensity: " << pCell->custom_data.variables[bi].value << std::endl;
//			std::cout << "Removing cell" << std::endl;
			pCell->start_death( phenotype.death.find_death_model_index( "apoptosis" ) );
		}
	}
	if( pCell->custom_data.variables[opt].value > 0 && pCell->custom_data.variables[gi].value > 0
					&& !pCell->phenotype.death.dead )
	{
		if( UniformRandom() <= 0.9*pCell->custom_data.variables[gi].value )
		{
//			std::cout << "Intensity: " << pCell->custom_data.variables[gi].value << std::endl;
//			std::cout << "Converting cell from " << pCell->type << " to ";
			pCell->convert_to_cell_definition( hek293_cell_B );
//			std::cout << pCell->type << std::endl;
		}
	}
/* ****************
 * END Optogenetics *
 ******************/

	// set the motility depending on the cell cycle phase
	if( pCell->phenotype.death.dead == true || phenotype.cycle.current_phase_index() == 3 )
	{
		phenotype.motility.is_motile = false;
	}
	else
	{
		phenotype.motility.is_motile = true;
	}
	return;
}

inline void update_cell_HEK293_B( Cell* pCell, Phenotype& phenotype, double dt )
{
	// first call the standard update function
	update_cell_and_death_parameters_O2_based( pCell, phenotype, dt );

	// get the current position
	Point P{pCell->position[0], pCell->position[1]};
//	std::cout << P.tuple() << std::endl;

	// current time
	double time = PhysiCell_globals.current_time;

	// Neighbors?
//	std::vector<Cell*> nearby = pCell->cells_in_my_container();

//	then do HEK293 cell specific things

/* **************
 * Optogenetics *
 ****************/
	// the phenotype governs all other properties like mechanics or cell cycle
	// set the internal state according to the light intensities
	if( pCell->custom_data.variables[opt].value > 0 )
	{
		pCell->custom_data.variables[bi].value = opto.get_intensity(P, time, index_blue);
		pCell->custom_data.variables[ri].value = opto.get_intensity(P, time, index_red);
		pCell->custom_data.variables[gi].value = opto.get_intensity(P, time, index_green);
	}
	if( pCell->custom_data.variables[opt].value > 0 && pCell->custom_data.variables[bi].value > 0
			&& !pCell->phenotype.death.dead )
	{
		if( UniformRandom() <= 0.9*pCell->custom_data.variables[bi].value )
		{
//			std::cout << "Intensity: " << pCell->custom_data.variables[bi].value << std::endl;
//			std::cout << "Removing cell" << std::endl;
			pCell->start_death( phenotype.death.find_death_model_index( "apoptosis" ) );
		}
	}
	if( pCell->custom_data.variables[opt].value > 0 && pCell->custom_data.variables[gi].value > 0
					&& !pCell->phenotype.death.dead )
	{
		if( UniformRandom() <= 0.9*pCell->custom_data.variables[gi].value )
		{
//			std::cout << "Intensity: " << pCell->custom_data.variables[gi].value << std::endl;
//			std::cout << "Converting cell from " << pCell->type << " to ";
			pCell->convert_to_cell_definition( hek293_cell_A );
//			std::cout << pCell->type << std::endl;
		}
	}
/* ****************
 * END Optogenetics *
 ******************/

	// set the motility depending on the cell cycle phase
	if( pCell->phenotype.death.dead == true || phenotype.cycle.current_phase_index() == 3 )
	{
		phenotype.motility.is_motile = false;
	}
	else
	{
		phenotype.motility.is_motile = true;
	}
	return;
}


inline void update_cell_HEK293_C( Cell* pCell, Phenotype& phenotype, double dt )
{
	// first call the standard update function
	update_cell_and_death_parameters_O2_based( pCell, phenotype, dt );

	// get the current position
	Point P{pCell->position[0], pCell->position[1]};
//	std::cout << P.tuple() << std::endl;

	// current time
	double time = PhysiCell_globals.current_time;

	// Neighbors?
//	std::vector<Cell*> nearby = pCell->cells_in_my_container();

//	then do HEK293 cell specific things

/* **************
 * Optogenetics *
 ****************/
	// the phenotype governs all other properties like mechanics or cell cycle
	// set the internal state according to the light intensities
	if( pCell->custom_data.variables[opt].value > 0 )
	{
		pCell->custom_data.variables[bi].value = opto.get_intensity(P, time, index_blue);
		pCell->custom_data.variables[ri].value = opto.get_intensity(P, time, index_red);
		pCell->custom_data.variables[gi].value = opto.get_intensity(P, time, index_green);
	}
	if( pCell->custom_data.variables[opt].value > 0 && pCell->custom_data.variables[bi].value > 0
			&& !pCell->phenotype.death.dead )
	{
		if( UniformRandom() <= 0.9*pCell->custom_data.variables[bi].value )
		{
//			std::cout << "Intensity: " << pCell->custom_data.variables[bi].value << std::endl;
//			std::cout << "Removing cell" << std::endl;
			pCell->start_death( phenotype.death.find_death_model_index( "apoptosis" ) );
		}
	}
	if( pCell->custom_data.variables[opt].value > 0 && pCell->custom_data.variables[gi].value > 0
					&& !pCell->phenotype.death.dead )
	{
		if( UniformRandom() <= 0.9*pCell->custom_data.variables[gi].value )
		{
//			std::cout << "Intensity: " << pCell->custom_data.variables[gi].value << std::endl;
//			std::cout << "Converting cell from " << pCell->type << " to ";
			pCell->convert_to_cell_definition( hek293_cell_A );
//			std::cout << pCell->type << std::endl;
		}
	}
/* ****************
 * END Optogenetics *
 ******************/

	// set the motility depending on the cell cycle phase
	if( pCell->phenotype.death.dead == true || phenotype.cycle.current_phase_index() == 3 )
	{
		phenotype.motility.is_motile = false;
	}
	else
	{
		phenotype.motility.is_motile = true;
	}
	return;
}

/*
 * This function is called from the function Cell_Container::update_all_cells.
 *
 */
inline void cell_rule_HEK293_A( Cell* pCell, Phenotype& phenotype, double dt )
{

	// Get first some information
	// Neighbors?
	std::vector<Cell*> nearby = pCell->cells_in_my_container();

	// count only the cells not already marked for death
	int cell_number = 0;
	for ( auto iter = nearby.begin(); iter < nearby.end(); ++iter )
	{
		if( !(*iter)->phenotype.death.dead && (*iter)->type == 1 ){ ++cell_number; }
	}

	// get the current position
//	Point P{pCell->position[0], pCell->position[1]};

	// current time
//	double time = PhysiCell_globals.current_time;

	// Default mechanics
	double default_adhesion = hek293_cell_A.phenotype.mechanics.cell_cell_adhesion_strength;
	double default_speed = hek293_cell_A.phenotype.motility.migration_speed;

	// Change the mechanical properties based on the cellular environment.
	// This is either to mimic mechanical contact dependent behavior or
	// secretion of some chemical and concentration dependence. In the latter case
	// should be replaced by a microenvironment sampling function

	/* **************
	 * Optogenetics *
	 ****************/
	/*
	 * Currently the optogenetics is in three places: update_cell, arrest_function and cell_rule.
	 * In phenotype all cell internal states like cell cycle, apoptosis, etc. and in
	 * cell_rule for the mechanical properties, like adhesion.
	 */

	// if there is red light, no  colonial growth
//	if( pCell->custom_data.variables[ri].value > 0.0 )
//	{
//		pCell->phenotype.motility.migration_speed = default_speed;
//		pCell->phenotype.mechanics.cell_cell_adhesion_strength = 0.5*default_adhesion;
//	}
//	// without red light, the cellular environment decides on the behaviour
//	else
//	{
		if ( (3 >= cell_number) && (cell_number > 2) )
		{
			phenotype.mechanics.cell_cell_adhesion_strength =  1.2*default_adhesion;
			phenotype.motility.migration_speed = 0.4*default_speed;
		}
		else if ( (cell_number > 4) )
		{
			phenotype.mechanics.cell_cell_adhesion_strength =  1.8*default_adhesion;
			phenotype.motility.migration_speed = 0.1*default_speed;
		}
		else
		{
			phenotype.mechanics.cell_cell_adhesion_strength = default_adhesion;
			phenotype.motility.migration_speed = default_speed;
		}
////	}

//	if( (time >= 1*24*60) && domain1.point_inside( P ) )
//	{
//		phenotype.mechanics.cell_cell_adhesion_strength = 2.8 * default_adhesion;
//		phenotype.motility.migration_speed = defaul_speed;
//	}
//	else
//	{
//		phenotype.mechanics.cell_cell_adhesion_strength = default_adhesion;
//		phenotype.motility.migration_speed = defaul_speed;
//	}

	return;
}

inline void cell_rule_HEK293_B( Cell* pCell, Phenotype& phenotype, double dt )
{

	// Get first some information
	// Neighbors?
	std::vector<Cell*> nearby = pCell->cells_in_my_container();

	// count only the cells not already marked for death
	int cell_number = 0;
	for ( auto iter = nearby.begin(); iter < nearby.end(); ++iter )
	{
		if( !(*iter)->phenotype.death.dead && (*iter)->type == 2 ){ ++cell_number; }
	}

	// get the current position
//	Point P{pCell->position[0], pCell->position[1]};

	// current time
//	double time = PhysiCell_globals.current_time;

	// Default mechanics
	double default_adhesion = hek293_cell_B.phenotype.mechanics.cell_cell_adhesion_strength;
	double default_speed = hek293_cell_B.phenotype.motility.migration_speed;

	// Change the mechanical properties based on the cellular environment.
	// This is either to mimic mechanical contact dependent behavior or
	// secretion of some chemical and concentration dependence. In the latter case
	// should be replaced by a microenvironment sampling function

	/* **************
	 * Optogenetics *
	 ****************/
	/*
	 * Currently the optogenetics is in three places: update_cell, arrest_function and cell_rule.
	 * In phenotype all cell internal states like cell cycle, apoptosis, etc. and in
	 * cell_rule for the mechanical properties, like adhesion.
	 */

	// if there is red light, no  colonial growth
//	if( pCell->custom_data.variables[ri].value > 0.0 )
//	{
//		pCell->phenotype.motility.migration_speed = default_speed;
//		pCell->phenotype.mechanics.cell_cell_adhesion_strength = 0.5*default_adhesion;
//	}
//	// without red light, the cellular environment decides on the behaviour
//	else
//	{
		if ( (3 >= cell_number) && (cell_number > 2) )
		{
			phenotype.mechanics.cell_cell_adhesion_strength =  1.2*default_adhesion;
			phenotype.motility.migration_speed = 0.4*default_speed;
		}
		else if ( (cell_number > 4) )
		{
			phenotype.mechanics.cell_cell_adhesion_strength =  1.8*default_adhesion;
			phenotype.motility.migration_speed = 0.1*default_speed;
		}
		else
		{
			phenotype.mechanics.cell_cell_adhesion_strength = default_adhesion;
			phenotype.motility.migration_speed = default_speed;
		}
////	}

//	if( (time >= 1*24*60) && domain1.point_inside( P ) )
//	{
//		phenotype.mechanics.cell_cell_adhesion_strength = 2.8 * default_adhesion;
//		phenotype.motility.migration_speed = defaul_speed;
//	}
//	else
//	{
//		phenotype.mechanics.cell_cell_adhesion_strength = default_adhesion;
//		phenotype.motility.migration_speed = defaul_speed;
//	}

	return;
}

inline void cell_rule_HEK293_C( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Get first some information
	// Neighbors?
	std::vector<Cell*> nearby = pCell->cells_in_my_container();

	// count only the cells not already marked for death
	int cell_number = 0;
	for ( auto iter = nearby.begin(); iter < nearby.end(); ++iter )
	{
		if( !(*iter)->phenotype.death.dead && (*iter)->type == 3 ){ ++cell_number; }
	}

	// get the current position
//	Point P{pCell->position[0], pCell->position[1]};

	// current time
//	double time = PhysiCell_globals.current_time;

	// Default mechanics
	double default_adhesion = hek293_cell_C.phenotype.mechanics.cell_cell_adhesion_strength;
	double default_speed = hek293_cell_C.phenotype.motility.migration_speed;

	// Change the mechanical properties based on the cellular environment.
	// This is either to mimic mechanical contact dependent behavior or
	// secretion of some chemical and concentration dependence. In the latter case
	// should be replaced by a microenvironment sampling function

	/* **************
	 * Optogenetics *
	 ****************/
	/*
	 * Currently the optogenetics is in three places: update_cell, arrest_function and cell_rule.
	 * In phenotype all cell internal states like cell cycle, apoptosis, etc. and in
	 * cell_rule for the mechanical properties, like adhesion.
	 */

	// if there is red light, no  colonial growth
//	if( pCell->custom_data.variables[ri].value > 0.0 )
//	{
//		pCell->phenotype.motility.migration_speed = default_speed;
//		pCell->phenotype.mechanics.cell_cell_adhesion_strength = 0.5*default_adhesion;
//	}
//	// without red light, the cellular environment decides on the behaviour
//	else
//	{
		if ( (3 >= cell_number) && (cell_number > 2) )
		{
			phenotype.mechanics.cell_cell_adhesion_strength =  1.2*default_adhesion;
			phenotype.motility.migration_speed = 0.4*default_speed;
		}
		else if ( (cell_number > 4) )
		{
			phenotype.mechanics.cell_cell_adhesion_strength =  1.8*default_adhesion;
			phenotype.motility.migration_speed = 0.1*default_speed;
		}
		else
		{
			phenotype.mechanics.cell_cell_adhesion_strength = default_adhesion;
			phenotype.motility.migration_speed = default_speed;
		}
////	}

//	if( (time >= 1*24*60) && domain1.point_inside( P ) )
//	{
//		phenotype.mechanics.cell_cell_adhesion_strength = 2.8 * default_adhesion;
//		phenotype.motility.migration_speed = defaul_speed;
//	}
//	else
//	{
//		phenotype.mechanics.cell_cell_adhesion_strength = default_adhesion;
//		phenotype.motility.migration_speed = defaul_speed;
//	}

	return;
}

// I put the intracelluar update in here because of the shorter update time of cell mechanics (0.1) compared to phenotype (6)
// Everything related to intracellular components should better be updated with the diffusion time (0.01)?
inline void cell_rule_HEK293_antithetic( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Get first some information
	// Neighbors?
	std::vector<Cell*> nearby = pCell->cells_in_my_container();

	// count only the cells not already marked for death
	int non_dead_cells = 0;
	for ( auto iter = nearby.begin(); iter < nearby.end(); ++iter )
	{
		if( !(*iter)->phenotype.death.dead ){ non_dead_cells++; }
	}

	int iZ1 = pCell->custom_data.find_variable_index( "Z1" );
	int iZ2 = pCell->custom_data.find_variable_index( "Z2" );

//	std::cout << "Z1 = " << pCell->custom_data[iZ1] << "Z2 = " << pCell->custom_data[iZ1] << std::endl;
//	std::cout << std::endl;

//	std::cout <<  pCell->custom_data["eta_annihilation"] *
//			pCell->custom_data[iZ1] * pCell->custom_data[iZ2] * mechanics_dt << std::endl;
	// current time
//	double time = PhysiCell_globals.current_time;

	// update Z1
	if( UniformRandom() <= pCell->custom_data["mu_Z1"] * mechanics_dt )
	{
		pCell->custom_data[iZ1]++;
//		std::cout << "Z1 = " << pCell->custom_data[iZ1] << "Z2 = " << pCell->custom_data[iZ1] << std::endl;
//		std::cout << std::endl;
	}
//	if( UniformRandom() <= 0.001 * pCell->custom_data[iZ1] * mechanics_dt )
//	{
//			pCell->custom_data[iZ1]--;
//	}
	// update Z2
	if( UniformRandom() <= pCell->custom_data["theta_Z2"] * (non_dead_cells-1) * mechanics_dt )
	{
		pCell->custom_data[iZ2]++;
	}
//	if( UniformRandom() <= 0.001 * pCell->custom_data[iZ2] * mechanics_dt )
//	{
//				pCell->custom_data[iZ2]--;
//	}
	// annihilation step
	if( UniformRandom() <= pCell->custom_data["eta_annihilation"] *
			pCell->custom_data[iZ1] * pCell->custom_data[iZ2] * mechanics_dt )
	{
		pCell->custom_data[iZ1]--;
		pCell->custom_data[iZ2]--;
	}
//	if( pCell->custom_data[iZ2] > 4 )
//	{
//		std::cout << "Z1 = " << pCell->custom_data[iZ1] << " Z2 = " << pCell->custom_data[iZ2] << std::endl;
//		std::cout << std::endl;
//	}
	return;
}

inline bool arrest_function_HEK293( Cell* pCell, Phenotype& phenotype, double dt )
{
	bool arrest = false;

/*
	// Get first some information
	// Neighbors?
	std::vector<Cell*> nearby = pCell->cells_in_my_container();

	// based on the voxel size 30x30 (don't know how to access this here) and cell's radius
	// I calculate a confluency level
	double r = pCell->phenotype.geometry.radius;
	double cl = nearby.size()*r*r*3.1415/(0.8*900.0);
	std::cout << cl << " " << nearby.size() << std::endl;
	if( cl >= 1 )
	{
//		 	 std::cout << r << " " << cl << " " << nearby.size() << std::endl;
		if( UniformRandom() <= 0.8 ){ _HI_; arrest = true; }
	}
*/

//	// current time
//	double time = PhysiCell_globals.current_time;

	/***************
	 * Optogenetics *
	****************/
	// the update_phenotype function is called first
	// in this the light intensities are updated
	// we do not need to check here
	if ( pCell->custom_data.variables[ri].value > 0.0025)
	{
		if( UniformRandom() <= 0.95*pCell->custom_data.variables[ri].value )
		{
//			std::cout << "Stopping proliferation" << std::endl;
			arrest = true;
		}
	}
	return arrest;
}
