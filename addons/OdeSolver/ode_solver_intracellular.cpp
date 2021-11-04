#include "ode_solver_intracellular.h"

#include <sstream>
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <fstream>

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
	pugi::xml_node node_mathml = node.child( "mathml_file" );
	if ( node_mathml )
	{
		mathml_filename = PhysiCell::xml_get_my_string_value (node_mathml);
		std::cout << "\n------------- "  << __FUNCTION__ << ": mathml_file = " << mathml_filename << std::endl;
		readMathMLFile(mathml_filename);
	}

	pugi::xml_node time = node.child( "intracellular_dt" );
	if ( time )
	{
		intracellular_dt = xml_get_my_double_value(time);
		std::cout << "\n------------- "  << __FUNCTION__ << ": intracellular_dt = " << intracellular_dt << std::endl;
	}

	int i = 0;
	pugi::xml_node node_substrate = node.child( "substrate" );
	while( node_substrate )
	{
		pugi::xml_node node_species = node_substrate.child( "map" );
        // ---------  substrates
		std::string substrate_name = node_species.attribute( "PC_substrate" ).value();
		if( substrate_name != "" )
		{
            //std::cout << "-----------" << node_species.attribute( "sbml_species" ).value() << std::endl;
			std::string species_name = node_species.attribute( "mathml_species" ).value();
			substrate_species[substrate_name] = species_name;
			species_name_to_index[species_name] = i;
			substrate_name_to_index[substrate_name] = i;
			index_to_substrate_name[i] = substrate_name;

			std::cout << "\n------------- "  << __FUNCTION__ << ": species_name= " << species_name << std::endl;

		}
		pugi::xml_node initial_value = node_substrate.child( "initial_condition" );
		if ( initial_value )
		{
			std::cout << "              " << "initial_condition= " << xml_get_my_double_value(initial_value) << std::endl;
			substrate_values.push_back(xml_get_my_double_value(initial_value));
		} else {
			std::cout << "              " << "initial_condition undefined setting = 0" << std::endl;
			substrate_values.push_back(0);
		}
		i++;

		node_substrate = node_substrate.next_sibling( "substrate" );
	}

    std::cout << "\n------------- substrate_species map:"  << std::endl;
    for(auto elm : substrate_species)
    {
        std::cout << "      "  << elm.first << " -> " << elm.second << std::endl;
    }
    std::cout << std::endl;

	return;
}

void OdeSolverIntracellular::readMathMLFile(std::string filename)
{
	pugi::xml_document ode_solver_math_ml_file;
	pugi::xml_node ode_solver_math_ml_root;

	std::cout << "Using config file " << filename << " ... " << std::endl ;
	pugi::xml_parse_result result = ode_solver_math_ml_file.load_file( filename.c_str()  );

	if( result.status != pugi::xml_parse_status::status_ok )
	{
		std::cout << "Error loading " << filename << "!" << std::endl;
	}
	ode_solver_math_ml_root = ode_solver_math_ml_file.child("ODE-Definition");
	std::cout << "Success loading MathML file " << filename << std::endl;
	return;
}

void OdeSolverIntracellular::start()
{
	std::cout << "Substrate_Values at Start()" << std::endl;
	std::cout << substrate_values[0] << "," << substrate_values[1] << std::endl;
	return;
}

bool OdeSolverIntracellular::need_update()
{
	return true;
}

// solve the intracellular model
void OdeSolverIntracellular::update()
{
	size_t steps = integrate( update_RHS , substrate_values , PhysiCell_globals.current_time , PhysiCell_globals.current_time + diffusion_dt , intracellular_dt );
//	std::cout << substrate_values << std::endl;
	std::ofstream outfile;
	outfile.open("logs.txt", std::ios_base::app); // append instead of overwrite
	outfile << substrate_values << " " << steps << std::endl;
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
int OdeSolverIntracellular::update_phenotype_parameters(PhysiCell::Phenotype& phenotype)
{
	return 1;
}

int OdeSolverIntracellular::validate_PhysiCell_tokens(PhysiCell::Phenotype& phenotype)
{
	return 1;
}

int OdeSolverIntracellular::validate_SBML_species()
{
	return 1;
}

int OdeSolverIntracellular::create_custom_data_for_SBML(PhysiCell::Phenotype& phenotype)
{
	return 0;
}

// ================  specific to "odeSolver" ================
void OdeSolverIntracellular::setUpdateFunction(update_func update_RHS_func)
{
	update_RHS = update_RHS_func;
	return;
}
