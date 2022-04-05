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
Val Diff_ObservableCuboid::measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) {
    int N_cells = std::count_if(
        cells.begin(),
        cells.end(),
        [this](PhysiCell::Cell* cell){
            // Check that cell is not dead
            bool cell_not_dead = !cell->phenotype.death.dead;
            // Check if cells are differentiated
            bool cell_correct_type = fabs(cell->custom_data["diff"] - diff_enable_1) < 0.1;
            return cell_not_dead && cell_correct_type;
        }
        );
    // std::cout << "matching cells in domain: " << N_cells << std::endl;
    return N_cells;
}


double Diff_Metric::calculate(Val& target, Val& observed) {
    return (observed-target)/target;
}


double PI_Controllfunctor::adjust(std::deque<double> state) {
    // Target is implicitly always 0.0
    // This implements a PI Controller
    double calculated = K_p*state.back() + K_i*(state.back()-state[state.size()-1])/update_dt;
    return calculated;
}


void Diff_Effect::apply(PhysiCell::Cell* cell, const double discrepancy) {
    // Check if cells meets criteria
    bool cell_diff_disable = fabs(cell->custom_data["diff"] - diff_disable) < 0.1;
    bool cell_in_diff_1 = fabs(cell->custom_data["diff"] - diff_enable_1) < 0.1;
    bool cell_not_dead = !cell->phenotype.death.dead;

    // Get a random parameter between 0 and 1
    double rand = PhysiCell::UniformRandom();
    // Determine threshhold via discrepancy
    double thresh = 1/fabs(discrepancy);
    
    // Check that the cell is not already dead
    if (cell_not_dead) {
        if (discrepancy > 0.0) {
            // If discrepancy is positive (more differentiated cells than targeted) kill some cells
            if (rand > thresh && cell_in_diff_1) {
                cell->start_death(cell->phenotype.death.find_death_model_index(100));
            }
        } else if (
            // Check if the cell is not yet differentiated
            cell_diff_disable && rand > 1-fabs(discrepancy)/(1+fabs(discrepancy))
        ) {
            // If the cell is not differentiated yet, the probability
            // to differentiate is given by the discrepancy
            cell->custom_data["diff"] = diff_enable_1;
            cell->phenotype.motility.migration_speed = new_migr_speed;
            cell->phenotype.mechanics.cell_cell_adhesion_strength = new_adh_strength;
        }
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
    
    int size=1;
    for (int i=0; i<size; i++) 
    {
        std::array<Kernel::Iso_cuboid_3, 2> _doms = {
            // Observable domain
            Kernel::Iso_cuboid_3(
                Kernel::Point_3(
                    BioFVM::microenvironment.mesh.bounding_box[0],
                    BioFVM::microenvironment.mesh.bounding_box[1] + 
                        PhysiCell::parameters.doubles("fraction_box_height")
                        *(BioFVM::microenvironment.mesh.bounding_box[4] - BioFVM::microenvironment.mesh.bounding_box[1]),
                    BioFVM::microenvironment.mesh.bounding_box[2]
                ),
                Kernel::Point_3(
                    BioFVM::microenvironment.mesh.bounding_box[3],
                    BioFVM::microenvironment.mesh.bounding_box[4],
                    BioFVM::microenvironment.mesh.bounding_box[5]
                )
            ),
            // Effect domain
            _doms[1] = 
            Kernel::Iso_cuboid_3(
                Kernel::Point_3(
                    BioFVM::microenvironment.mesh.bounding_box[0],
                    BioFVM::microenvironment.mesh.bounding_box[1],
                    BioFVM::microenvironment.mesh.bounding_box[2]
                ),
                Kernel::Point_3(
                    BioFVM::microenvironment.mesh.bounding_box[3],
                    BioFVM::microenvironment.mesh.bounding_box[1] + 
                        PhysiCell::parameters.doubles("fraction_box_height")
                        *(BioFVM::microenvironment.mesh.bounding_box[4] - BioFVM::microenvironment.mesh.bounding_box[1]),
                    BioFVM::microenvironment.mesh.bounding_box[5]
                )
            )
        };
        diff_domains.push_back(_doms);
    };

    int i = 0;
    for (auto const& [_observable_domain, _effect_domain] : diff_domains) {
        auto cont = new Diff_Controller(
            // Observable and Effect Domain
            // Kernel::Iso_cuboid_3(Kernel::Point_3(-150.0, -150.0, -10.0), Kernel::Point_3(150.0, 150.0, 10.0)),
            _observable_domain,
            _effect_domain,
            // Target: How many differentiated cells do we want to have?
            PhysiCell::parameters.doubles("Target_cells_per_domain_diff_1")
        );
        supervisor.add_controller("Diff_Controller_diff_1_" + std::to_string(i), cont);
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