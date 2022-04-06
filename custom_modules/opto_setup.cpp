#include "opto_setup.h"
#include <cmath>
#include <random>
#include <chrono>

double optogenetics_dt = 0.0;
double optogenetics_next_run_time = 0.0;
double optogenetics_update_dt = 0.0;
double optogenetics_next_update_time = 0.0;
double minimum_time_interval = 0.0;

// Create a supervisor to monitor all controllers
Opto::Controller::Supervisor supervisor;

// *********************************************************************************
// RUN OPTOGENETIC ROUTINES IN MAIN
void shine_light( const double& t) {
    // std::cout << "[Main] Shining light at time " << t << std::endl;
    supervisor.run_all_controllers();
    return;
}


void run_optogenetics ( const double& t ) {
    if (t-optogenetics_next_run_time + 0.01*minimum_time_interval > 0) {
        optogenetics_next_run_time+=optogenetics_dt;
        shine_light(t);

    } else if (t-optogenetics_next_update_time + 0.01*minimum_time_interval > 0) {
        optogenetics_next_update_time += optogenetics_update_dt;
        // Additionally run opdates for optogenetic controllers
        // std::cout << "[Main] Running opto update at time " << t << std::endl;
        supervisor.update_all_controllers();
    }
}


// *********************************************************************************
// DIFFERENTIATION CONTROLLER MODULES
Val Diff_1_ObservableCuboid::measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) {
    int N_cells = std::count_if(
        cells.begin(),
        cells.end(),
        [this](PhysiCell::Cell* cell){
            // Check that cell is not dead
            bool cell_not_dead = !cell->phenotype.death.dead;
            // Check if cells are differentiated
            bool cell_correct_type_1 = cell->type_name == "differentiation_cell";
            bool cell_correct_type_2 = false;
            if (cell_correct_type_1) {
                cell_correct_type_2 = fabs(cell->custom_data["diff"] - diff_enable_1) < 0.1;
            }
            return cell_not_dead && cell_correct_type_1 && cell_correct_type_2;
        }
        );
    // std::cout << "matching cells in domain: " << N_cells << std::endl;
    return N_cells;
}


Val Diff_2_ObservableCuboid::measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) {
    int N_cells = std::count_if(
        cells.begin(),
        cells.end(),
        [this](PhysiCell::Cell* cell){
            // Check that cell is not dead
            bool cell_not_dead = !cell->phenotype.death.dead;
            // Check if cells are differentiated
            bool cell_correct_type_1 = cell->type_name == "differentiation_cell";
            bool cell_correct_type_2 = false;
            if (cell_correct_type_1) {
                cell_correct_type_2 = fabs(cell->custom_data["diff"] - diff_enable_2) < 0.1;
            }
            return cell_not_dead && cell_correct_type_1 && cell_correct_type_2;
        }
        );
    // std::cout << "matching cells in domain: " << N_cells << std::endl;
    return N_cells;
}


Val Diff_3_ObservableCuboid::measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) {
    int N_cells = std::count_if(
        cells.begin(),
        cells.end(),
        [this](PhysiCell::Cell* cell){
            // Check that cell is not dead
            bool cell_not_dead = !cell->phenotype.death.dead;
            // Check if cells are differentiated
            bool cell_correct_type_1 = cell->type_name == "differentiation_cell";
            bool cell_correct_type_2 = false;
            if (cell_correct_type_1) {
                cell_correct_type_2 = fabs(cell->custom_data["diff"] - diff_enable_3) < 0.1;
            }
            return cell_not_dead && cell_correct_type_1 && cell_correct_type_2;
        }
        );
    // std::cout << "matching cells in domain: " << N_cells << std::endl;
    return N_cells;
}


double Diff_Metric::calculate(Val& target, Val& observed) {
    return (target-observed)/target;
}


double PID_Controllfunctor::adjust(std::deque<double> state) {
    // Target is implicitly always 0.0
    // This implements a PID Controller
    double calculated = K_p*state.back();
    if (state.size()>1) {
        calculated += K_d*(state.back()-state[state.size()-1])/update_dt;
    }
    if (state.size()>1) {
        calculated += K_i*update_dt*std::accumulate(
            state.begin(),
            state.end(),
            0.0
        );
    }
    std::cout << "Calculated: " << calculated << "\n";
    return calculated;
}


void Diff_1_Effect::apply(PhysiCell::Cell* cell, const double discrepancy) {
    // Check if cells meets criteria
    bool cell_not_dead = !cell->phenotype.death.dead;
    bool cell_correct_type = cell->type_name == "differentiation_cell";

    if (cell_correct_type) {
        cell->phenotype.intracellular->set_parameter_value(00, std::max(discrepancy*10000.0, 0.0));
        cell->phenotype.intracellular->set_parameter_value(20, -std::min(discrepancy*1000.0, 0.0));
    }
    return;
}


void Diff_2_Effect::apply(PhysiCell::Cell* cell, const double discrepancy) {
    // Check if cells meets criteria
    bool cell_not_dead = !cell->phenotype.death.dead;
    bool cell_correct_type = cell->type_name == "differentiation_cell";

    if (cell_correct_type) {
        cell->phenotype.intracellular->set_parameter_value(10, std::max(discrepancy*10000.0, 0.0));
        cell->phenotype.intracellular->set_parameter_value(20, -std::min(discrepancy*1000.0, 0.0));
    }
    return;
}


void Diff_3_Effect::apply(PhysiCell::Cell* cell, const double discrepancy) {
    // Check if cells meets criteria
    bool cell_not_dead = !cell->phenotype.death.dead;
    bool cell_correct_type = cell->type_name == "differentiation_cell";

    if (cell_correct_type) {
        // TODO this is currently meaningless
        cell->phenotype.intracellular->set_parameter_value(20,std:: max(discrepancy*10000.0, 0.0));
        cell->phenotype.intracellular->set_parameter_value(20, -std::min(discrepancy*1000.0, 0.0));
    }
    return;
}


// *********************************************************************************
// SETUP OPTOGENETIC CONTROL FRAMEWORK (PRE-MAIN ROUTINE)
void setup_optogenetics( void ) {
    optogenetics_dt = PhysiCell::parameters.doubles("optogenetics_dt");
    optogenetics_update_dt = PhysiCell::parameters.doubles("optogenetics_update_dt");
    // TODO check if this works even when setting mechanics dt in xml file lower
    minimum_time_interval = std::min(std::min(optogenetics_dt, optogenetics_update_dt), PhysiCell::mechanics_dt);

	std::cout << "\n\n\n**************************************" << std::endl;
    std::cout << "          Setting up OptoGen          " << std::endl;
    std::cout << "**************************************\n" << std::endl;

    // Create individual domains for controll and effect
    std::vector<std::array<Kernel::Iso_cuboid_3, 2>> diff_domains{};
    
    double x0 = BioFVM::microenvironment.mesh.bounding_box[0];
    double y0 = BioFVM::microenvironment.mesh.bounding_box[1];
    double z0 = BioFVM::microenvironment.mesh.bounding_box[2];

    double x1 = BioFVM::microenvironment.mesh.bounding_box[3];
    double y1 = BioFVM::microenvironment.mesh.bounding_box[4];
    double z1 = BioFVM::microenvironment.mesh.bounding_box[5];

    double frac = PhysiCell::parameters.doubles("fraction_box_height");

    double size=PhysiCell::parameters.doubles("Number_of_domains");
    for (double i=0; i<size; i++) 
    {
        std::array<Kernel::Iso_cuboid_3, 2> _doms = {
            // Observable domain
            Kernel::Iso_cuboid_3(
                Kernel::Point_3(x0 + i/size*(x1-x0), y0 + frac*(y1-y0), z0),
                Kernel::Point_3(x0 + (i+1.0)/size*(x1-x0), y1, z1)
            ),
            // Effect domain
            _doms[1] = 
            Kernel::Iso_cuboid_3(
                Kernel::Point_3(x0 + i/size*(x1-x0), y0, z0),
                Kernel::Point_3(x0 + (i+1.0)/size*(x1-x0), y0 + frac*(y1-y0), z1)
            )
        };
        diff_domains.push_back(_doms);
    };

    int i = 0;
    for (auto const& [_observable_domain, _effect_domain] : diff_domains) {
        auto cont1 = new Diff_1_Controller(
            // Observable and Effect Domain
            // Kernel::Iso_cuboid_3(Kernel::Point_3(-150.0, -150.0, -10.0), Kernel::Point_3(150.0, 150.0, 10.0)),
            _observable_domain,
            _effect_domain,
            // Target: How many differentiated cells do we want to have?
            PhysiCell::parameters.doubles("Target_cells_per_domain_diff_1")
        );
        supervisor.add_controller("Diff_Controller_diff_1_" + std::to_string(i), cont1);

        auto cont2 = new Diff_2_Controller(
            // Observable and Effect Domain
            // Kernel::Iso_cuboid_3(Kernel::Point_3(-150.0, -150.0, -10.0), Kernel::Point_3(150.0, 150.0, 10.0)),
            _observable_domain,
            _effect_domain,
            // Target: How many differentiated cells do we want to have?
            PhysiCell::parameters.doubles("Target_cells_per_domain_diff_2")
        );
        supervisor.add_controller("Diff_Controller_diff_2_" + std::to_string(i), cont2);

        /* auto cont3 = new Diff_3_Controller(
            // Observable and Effect Domain
            // Kernel::Iso_cuboid_3(Kernel::Point_3(-150.0, -150.0, -10.0), Kernel::Point_3(150.0, 150.0, 10.0)),
            _observable_domain,
            _effect_domain,
            // Target: How many differentiated cells do we want to have?
            PhysiCell::parameters.doubles("Target_cells_per_domain_diff_3")
        );
        supervisor.add_controller("Diff_Controller_diff_3_" + std::to_string(i), cont3);*/
        i++;
    }

    // TODO CURRENTLY NOT WORKING
    // Add visitor to save state to csv file
    // supervisor.add_visitor_update("Safe_state", new Visitor_Write_State_to_csv("testfile.csv"));
    // supervisor.update_visit_all_controllers("Safe_state");

    // auto cb1 = Kernel::Iso_cuboid_3(Kernel::Point_3(-300.0, -300.0, -10.0), Kernel::Point_3(300.0, 300.0, 10.0));
    // auto cb2 = Kernel::Iso_cuboid_3(Kernel::Point_3(-100.0, 0.0, -10.0), Kernel::Point_3(250.0, 500.0, 10.0));

    std::cout << "**************************************\n" << std::endl;
}