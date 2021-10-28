#include "ode_solver_intracellular.h"

#include <sstream>
#include <iostream>

OdeSolverIntracellular::OdeSolverIntracellular() : Intracellular()
{
	intracellular_type = "ode_solver";
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

}

void OdeSolverIntracellular::initialize_intracellular_from_pugixml(pugi::xml_node& node)
{
	pugi::xml_node node_sbml = node.child( "mathml_file" );
		if ( node_sbml )
		{
			mathml_filename = PhysiCell::xml_get_my_string_value (node_sbml);
	        std::cout << "\n------------- "  << __FUNCTION__ << ": mathml_file = " << mathml_filename << std::endl;
	    }
	    return;
}


void OdeSolverIntracellular::start()
{
	return;
}

bool OdeSolverIntracellular::need_update()
{
	return true;
}

// solve the intracellular model
void OdeSolverIntracellular::update()
{
	std::cout << "Updating intracellular right now!" << std::endl;
	return;
}

double OdeSolverIntracellular::get_parameter_value(std::string param_name)
{
    return 0;
}
	
// rwh: might consider doing a multi-[species_name, value] "set" method
void OdeSolverIntracellular::set_parameter_value(std::string species_name, double value)
{
	std::cout << "Test 2 "<< species_name << " with value " << value << std::endl;
    return;
}

std::string OdeSolverIntracellular::get_state()
{
    return "hi";
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
