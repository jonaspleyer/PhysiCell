#ifndef _OdeSolver_Intracellular_h_
#define _OdeSolver_Intracellular_h_

#include "../../core/PhysiCell.h"
#include "../../core/PhysiCell_phenotype.h"
#include "../../core/PhysiCell_cell.h"
#include "../../modules/PhysiCell_pugixml.h"

class OdeSolverIntracellular : public PhysiCell::Intracellular
{
 private:
 public:
	OdeSolverIntracellular();
	
	OdeSolverIntracellular(pugi::xml_node& node);

	OdeSolverIntracellular(OdeSolverIntracellular* copy);

	Intracellular* clone() {
		return static_cast<Intracellular*>(new OdeSolverIntracellular(this));
	}
	Intracellular* getIntracellularModel() {
		return static_cast<Intracellular*>(this);
	}

	std::string mathml_filename;

    std::string intracellular_type;  // specified in XML <intracellular type="...">:  "maboss", "sbml", ...

    // ================  generic  ================
	// This function parse the xml cell definition
	void initialize_intracellular_from_pugixml(pugi::xml_node& node);
	
	// This function initialize the model, needs to be called on each cell once created
	void start();
	
	// This function checks if it's time to update the model
	bool need_update();

	// This function update the model for the time_step defined in the xml definition
	void update();

	// Get value for model parameter
	double get_parameter_value(std::string name);
	
	// Set value for model parameter
	void set_parameter_value(std::string name, double value);

	std::string get_state();
	
    // ================  specific to "maboss" ================
	bool has_variable(std::string name);
	bool get_boolean_variable_value(std::string name);
	void set_boolean_variable_value(std::string name, bool value);
	void print_current_nodes();


    // ================  specific to "roadrunner" ================
    int update_phenotype_parameters(PhysiCell::Phenotype& phenotype);
    int validate_PhysiCell_tokens(PhysiCell::Phenotype& phenotype);
    int validate_SBML_species();
    int create_custom_data_for_SBML(PhysiCell::Phenotype& phenotype);

};

#endif
