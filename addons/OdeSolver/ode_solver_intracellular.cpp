#include "ode_solver_intracellular.h"

#include <sstream>
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace boost::numeric::odeint;
using namespace PhysiCell;

OdeSolverIntracellular::OdeSolverIntracellular() : Intracellular()
{
	intracellular_type = "ode_solver";
	intracellular_dt = 0.01;
}

OdeSolverIntracellular::OdeSolverIntracellular(pugi::xml_node& node)
{
	intracellular_type = "ode_solver";
	initialize_intracellular_from_pugixml(node);
}

OdeSolverIntracellular::OdeSolverIntracellular(OdeSolverIntracellular* copy)
{
	intracellular_type = copy->intracellular_type;
	mathml_filename = copy->mathml_filename;
	intracellular_dt = copy->intracellular_dt;
}

void OdeSolverIntracellular::initialize_intracellular_from_pugixml(pugi::xml_node& node)
{
	substrate_values.resize(microenvironment.density_names.size(),0);

	pugi::xml_node time = node.child( "intracellular_dt" );
	if ( time )
	{
		intracellular_dt = xml_get_my_double_value(time);
		std::cout << "\n------------- "  << __FUNCTION__ << ": intracellular_dt  = " << intracellular_dt << std::endl;
	}

	int ID = 0;
	pugi::xml_node node_substrate = node.child( "substrate_int" );
	while( node_substrate )
	{
//		ID = xml_get_my_int_value(node_substrate.child("ID"));
		pugi::xml_node node_species = node_substrate.child( "map" );
        // ---------  substrates
		std::string substrate_name = node_species.attribute( "substrate" ).value();
		ID = microenvironment.find_density_index(substrate_name);
		if ( ID == -1 )
		{
			ID = substrate_values.size()+1;
		}
		if( substrate_name != "" )
		{
            //std::cout << "-----------" << node_species.attribute( "sbml_species" ).value() << std::endl;
			std::string species_name = node_species.attribute( "mathml_species" ).value();
			substrate_name_to_species_name[substrate_name] = species_name;
			if ( species_name_to_index.find(species_name) == species_name_to_index.end() )
			{
				species_name_to_index[species_name] = ID;
			}
			else
			{
				std::cout << __FUNCTION__ << " ERROR: Key duplicate in mathml_species: " << species_name << std::endl;
		        throw std::exception();
			}

			if ( substrate_name_to_index.find(substrate_name) == substrate_name_to_index.end() )
			{
				substrate_name_to_index[substrate_name] = ID;
			}
			else
			{
				std::cout << __FUNCTION__ << " ERROR: Key duplicate in substrate: " << substrate_name << std::endl;
		        throw std::exception();
			}

			if ( index_to_substrate_name.find(ID) == index_to_substrate_name.end() )
			{
				index_to_substrate_name[ID] = substrate_name;
			}
			else
			{
				std::cout << __FUNCTION__ << " ERROR: Key duplicate in ID: " << ID << std::endl;
		        throw std::exception();
			}

			std::cout << "\n------------- "  << __FUNCTION__ << ": substrate_name    = " << substrate_name << std::endl;
			std::cout << "              "  << __FUNCTION__ << ": species_name      = " << species_name << std::endl;
			std::cout << "              "  << __FUNCTION__ << ": ID                = " << ID << std::endl;

		}
		else
		{
			std::cout << __FUNCTION__ << " ERROR: No substrate_name specified!" << std::endl;
			throw std::exception();
		}
		pugi::xml_node initial_value = node_substrate.child( "initial_condition_int" );

		double init_val;
		if ( initial_value )
		{
			init_val = xml_get_my_double_value(initial_value);
			std::cout << "              "  << __FUNCTION__ << ": initial_condition = " << substrate_values[ID] << std::endl;
		} else {
			init_val = 0;
			std::cout << "              "  << __FUNCTION__ << ": initial_condition = " << substrate_values[ID] << " (undefined: auto-assign)" << std::endl;
		}
		if (ID > substrate_values.size())
		{
			substrate_values.resize(ID, init_val);
		}
		else
		{
			substrate_values[ID] = init_val;
		}

		node_substrate = node_substrate.next_sibling( "substrate_int" );
	}

    std::cout << "\n------------- " << __FUNCTION__ << ": substrate_species map:" << std::endl;
    for(auto elm : substrate_name_to_species_name)
    {
        std::cout << "              "  << elm.first << " --> " << elm.second << std::endl;
    }

    pugi::xml_node node_mathml = node.child( "mathml_file" );
	if ( node_mathml )
	{
		mathml_filename = xml_get_my_string_value (node_mathml);
		std::cout << "\n------------- "  << __FUNCTION__ << ": mathml_file       = " << mathml_filename << std::endl;
		readMathMLFile(mathml_filename);
	}
	std::cout << std::endl;
	return;
}

void OdeSolverIntracellular::readMathMLFile(std::string filename)
{
	pugi::xml_document ode_solver_math_ml_file;
	pugi::xml_node ode_solver_math_ml_root;
	pugi::xml_parse_result result = ode_solver_math_ml_file.load_file( filename.c_str()  );

	if( result.status != pugi::xml_parse_status::status_ok )
	{
		std::cout << "              Error loading " << filename << "!" << std::endl;
	}
	else
	{
		std::cout << "              Success loading MathML file " << filename << std::endl;
	}
	ode_solver_math_ml_root = ode_solver_math_ml_file.child("ODE-Definition");
	return;
}

void OdeSolverIntracellular::start()
{
	N_int_vals = substrate_values.size();
	initialized = false;
	return;
}

void OdeSolverIntracellular::update_Cell_parameters(PhysiCell::Cell &cell)
{
	// Check if Intracellular was just previously initialized. If so do not copy values. Only used to write initial values.
	if (initialized == false)
	{
		N_ext_vals = cell.phenotype.molecular.internalized_total_substrates.size();
		initialized = true;

		// Just for testing
		std::ofstream outfile;
		std::string filename = "output/logs_" + std::to_string(cell.ID) + ".txt";
		outfile.open(filename, std::ios_base::app); // append instead of overwrite
		outfile << "steps ";
		for (int i=0; i<N_int_vals + N_ext_vals; i++)
		{
			if ( i < N_ext_vals)
			{
				outfile << "ext:" << cell.phenotype.secretion.pMicroenvironment->density_names[i] << " ";
			}
			else if ( i <  2*N_ext_vals)
			{
				outfile << "int:" << cell.phenotype.secretion.pMicroenvironment->density_names[i-N_ext_vals] << " ";
			}
			else
			{
				outfile << "int:" << index_to_substrate_name[i-N_ext_vals+1] << " ";
			}
		}
		outfile << std::endl;
		outfile.close();
		// Testing part end
	}
	else
	{
		if ( default_microenvironment_options.track_internalized_substrates_in_each_agent )
		{
			update_substrate_values(cell.phenotype.molecular.internalized_total_substrates);
		}
	}

	// Combine external and internal substrate values into one large vector and feed it into the ODE Solver
	int voxel_index = cell.get_current_voxel_index();
	std::vector<double> internal_substrate_values = substrate_values;
	std::vector<double> external_substrate_values = (*(cell.phenotype.secretion.pMicroenvironment))(voxel_index);
	// We have to multiply by volume since phenotype gives us only a density
	double voxel_volume = cell.phenotype.secretion.pMicroenvironment->voxels(voxel_index).volume;
	external_substrate_values *= voxel_volume;

	// Create a vector of combined substrates to feed into the update function
	std::vector<double> combined_substrate_values = external_substrate_values;
	combined_substrate_values.insert(combined_substrate_values.end(), substrate_values.begin(), substrate_values.end());

	// Actually integrate the RHS of the equation
	size_t steps = integrate( update_RHS , combined_substrate_values , PhysiCell_globals.current_time , PhysiCell_globals.current_time + diffusion_dt , intracellular_dt );

	// Split the combined vector into the internal and external part and update parameters
	substrate_values = std::vector<double>(combined_substrate_values.begin()+N_ext_vals ,combined_substrate_values.end());
	std::vector<double> external_change(combined_substrate_values.begin(), combined_substrate_values.begin()+N_ext_vals);
	external_change -= external_substrate_values;
	external_change /= voxel_volume;

	// Just for testing
	std::ofstream outfile;
	std::string filename = "output/logs_" + std::to_string(cell.ID) + ".txt";
	outfile.open(filename, std::ios_base::app); // append instead of overwrite
	outfile << steps << " " << combined_substrate_values << std::endl;
	outfile.close();

	// Update newly calculated parameters
	if ( default_microenvironment_options.track_internalized_substrates_in_each_agent )
	{
		update_internalized_substrates(cell.phenotype.molecular.internalized_total_substrates);
	}

	// Update external values and divide by voxel volume

	(*(cell.phenotype.secretion.pMicroenvironment))(voxel_index) += external_change;
	return;
}

void OdeSolverIntracellular::update_internalized_substrates(std::vector<double> &internalized_substrates)
{
	if ( N_ext_vals == internalized_substrates.size() )
	{
		internalized_substrates = std::vector<double>(substrate_values.begin(), substrate_values.begin()+N_ext_vals);
	}
	else
	{
		std::cout << __FUNCTION__ << ": ERROR: internalized_substrates vector does not match received size" << std::endl;
		throw std::exception();
	}
	return;
}

void OdeSolverIntracellular::update_substrate_values(std::vector<double> &internalized_substrates)
{
	std::vector<double> only_internal_vals = std::vector<double>(substrate_values.begin()+N_ext_vals, substrate_values.end());
	if ( only_internal_vals.size() != N_int_vals-N_ext_vals )
	{
		std::cout << __FUNCTION__ << ": ERROR: pure internal values does not match size" << std::endl;
		throw std::exception();
	}
	substrate_values = internalized_substrates;
	substrate_values.insert(substrate_values.end(), only_internal_vals.begin(), only_internal_vals.end());
	return;
}

bool OdeSolverIntracellular::need_update()
{
	if( *update_RHS == NULL )
	{
		return false;
	}
	else
	{
		return true;
	}
}

// solve the intracellular model
void OdeSolverIntracellular::update()
{
	return;
}


double OdeSolverIntracellular::get_parameter_value(std::string name)
{
	return substrate_values[substrate_name_to_index[name]];
}

double OdeSolverIntracellular::get_parameter_value(int index)
{
	return substrate_values[index];
}

// rwh: might consider doing a multi-[species_name, value] "set" method
void OdeSolverIntracellular::set_parameter_value(std::string species_name, double value)
{
	substrate_values[species_name_to_index[species_name]] = value;
    return;
}

void OdeSolverIntracellular::set_parameter_value(int index, double value)
{
	substrate_values[index] = value;
    return;
}

std::string OdeSolverIntracellular::get_state()
{
    return "initialized";
}
bool OdeSolverIntracellular::has_variable(std::string name)
{
	return true;
}

bool OdeSolverIntracellular::get_boolean_variable_value(std::string name)
{
	return false;
}

void OdeSolverIntracellular::set_boolean_variable_value(std::string name, bool value)
{
	return;
}

void OdeSolverIntracellular::print_current_nodes()
{
	return;
}

// ================  specific to "roadrunner" ================
int OdeSolverIntracellular::update_phenotype_parameters(Phenotype& phenotype)
{
	return 1;
}

int OdeSolverIntracellular::validate_PhysiCell_tokens(Phenotype& phenotype)
{
	return 1;
}

int OdeSolverIntracellular::validate_SBML_species()
{
	return 1;
}

int OdeSolverIntracellular::create_custom_data_for_SBML(Phenotype& phenotype)
{
	return 0;
}

// ================  specific to "odeSolver" ================
void OdeSolverIntracellular::setUpdateFunction(update_func update_RHS_func)
{
	update_RHS = update_RHS_func;
	return;
}
