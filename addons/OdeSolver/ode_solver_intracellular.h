#ifndef _OdeSolver_Intracellular_h_
#define _OdeSolver_Intracellular_h_

#include "../../core/PhysiCell.h"
#include "../../core/PhysiCell_phenotype.h"
#include "../../core/PhysiCell_cell.h"
#include "../../modules/PhysiCell_pugixml.h"

#include <map>

typedef void (*update_func)(const std::vector<double> &X, std::vector<double> &dX, const double dt);

class OdeSolverIntracellular : public PhysiCell::Intracellular
{
 private:
	// The following maps are used to jump between PhysiCell_settings.xml, the mathml_file and internal definitions
	// DEFINITIONS:
	// - substrate_name (PhysiCell_settings.xml): name of secreted substrate
	// - index (internal):
	// - species_name (mathml_file): name of species as specified in mathml_file
	std::map<std::string, std::string> substrate_name_to_species_name;
	std::map<std::string, int> species_name_to_index;
	std::map<std::string, int> substrate_name_to_index;
	std::map<int, std::string> index_to_substrate_name;
	// Actual values of internal substrates
	std::vector<double> substrate_values;

	bool initialized;
	int N_int_vals;
	int N_ext_vals;
 public:
	OdeSolverIntracellular();
	
	OdeSolverIntracellular(pugi::xml_node& node);

	OdeSolverIntracellular(OdeSolverIntracellular* copy);

	Intracellular* clone() {
		OdeSolverIntracellular* clone = new OdeSolverIntracellular(this);
		// Strings
		clone->mathml_filename = this->mathml_filename;
		clone->intracellular_type = this->intracellular_type;
		// Maps
		clone->substrate_name_to_species_name = this->substrate_name_to_species_name;
		clone->species_name_to_index = this->species_name_to_index;
		clone->substrate_name_to_index = this->substrate_name_to_index;
		clone->index_to_substrate_name = this->index_to_substrate_name;
		// Vectors
		clone->substrate_values = this->substrate_values;
		clone->intracellular_dt = this->intracellular_dt;
		// Functions
		clone->update_RHS = this->update_RHS;
		return static_cast<Intracellular*>(clone);
	}
	Intracellular* getIntracellularModel() {
        std::cout << "------ ode_solver_intracellular: getIntracellularModel called\n";
		return static_cast<Intracellular*>(this);
	}

	// Timescale on which to solve the given differential equations. May be altered by methods with non-constant stepsize approaches.
	double intracellular_dt;

	// Filename to obtain RHS of Differential equation from
	std::string mathml_filename;
	// Just description
    std::string intracellular_type;  // specified in XML <intracellular type="...">:  "maboss", "sbml", ...
    
    // ================  generic  ================
	// This function parse the xml cell definition
	void initialize_intracellular_from_pugixml(pugi::xml_node& node);
	
	// Get Information from mathml_file to specify RHS of ODE
	void readMathMLFile(std::string filename);

	// This function initialize the model, needs to be called on each cell once created
	void start();
	
	// This function checks if it's time to update the model
	bool need_update();

	// This function update the model for the time_step defined in the xml definition
	void update();

	// Get value for model parameter
	double get_parameter_value(std::string name);
	
	double get_parameter_value(int index);

	// Set value for model parameter
	void set_parameter_value(std::string name, double value);

	void set_parameter_value(int index, double value);

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

    // ================  specific to "odeSolver" ================
    void setUpdateFunction(update_func update_RHS_func);
    void update_Cell_parameters(PhysiCell::Cell &cell);
    // Pointer to the right-hand-side of the ODE we want to solve.
    void (*update_RHS)(const std::vector<double> &X, std::vector<double> &dX, const double dt);
    void update_internalized_substrates(std::vector<double> &internalized_substrates);
    void update_substrate_values(std::vector<double> &internalized_substrates);
};

#endif
