#include "opto_setup.h"
#include <cmath>
#include <random>
#include <stdio.h>

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
        optogenetics_next_run_time += optogenetics_dt;
        optogenetics_next_update_time += optogenetics_update_dt;
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
    int N_cells_diff_1 = std::count_if(
        cells.begin(),
        cells.end(),
        [this](PhysiCell::Cell* cell) {
            // Check that cell is not dead
            if (cell->phenotype.death.dead) {
                return false;
            } else if (cell->type_name != "differentiation_cell") {
                return false;
            } else {
                double diff = cell->custom_data["diff"];
                return fabs(diff - diff_enable_1) < 0.1;
            }
        }
        );
    return N_cells_diff_1;
}


void Diff_Effect::apply(PhysiCell::Cell* cell, const double discrepancy) {
    // Check if cells meets criteria
    bool cell_not_dead = !cell->phenotype.death.dead;
    bool cell_correct_type = cell->type_name == "default";

    // If discrepancy is positive, the ratio diff_1_cells to diff_2_cells is not large enough
    // Thus we need to increase production of substrate1
    // And reduce production of substrate2
    if (cell_not_dead && cell_correct_type) {
		// Secrete differentiating substrate
        double current_value_1 = cell->phenotype.intracellular->get_parameter_value(",00");
        // std::cout << current_value_1 << " " << current_value_2 << std::endl;
        double new_value_1 = std::max(current_value_1 + discrepancy, 0.0);
        cell->phenotype.intracellular->set_parameter_value(00, new_value_1);
		
		// Secrete killer substrate
        double current_value_killer = cell->phenotype.intracellular->get_parameter_value(",30");
        double new_value_killer = std::max(current_value_killer - discrepancy, 0.0);
		// std::cout << current_value_killer << " " << new_value_killer << "\n";
        cell->phenotype.intracellular->set_parameter_value(30, new_value_killer);
    }
    return;
}


double Diff_Metric::calculate(Val& target, Val& observed) {
    return target-observed;
}


double PID_Controllfunctor::adjust(std::deque<double> state) {
    // Target is implicitly always 0.0
    // This implements a PI Controller
    double calculated = K_p*state.back();
    double prop = calculated;
    double differential = 0;
    double integral = 0;
    
    if (state.size()>1) {
        differential = K_d*(state.back()-state.end()[-2])/update_dt;
        calculated += differential;
    }
    if (state.size()>10) {
        integral = K_i*update_dt*std::accumulate(
            state.begin(),
            state.end(),
            0.0
        );
        calculated += integral;
    }
    auto fp = fopen("controller_logs.txt", "a");
    // fprintf(stderr, "Prop: %E Diff: %E Int: %E Calc: %E sSize: %2d LastState: %2E\n", prop, differential, integral, calculated, state.size(), state.back());
    fprintf(fp, "%2d, %E, %E, %E, %E, %E\n", state.size(), prop, differential, integral, calculated, state.back());
    fclose(fp);
    return calculated;
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
    
    double x0 = BioFVM::microenvironment.mesh.bounding_box[0] - BioFVM::microenvironment.mesh.dx;
    double y0 = BioFVM::microenvironment.mesh.bounding_box[1] - BioFVM::microenvironment.mesh.dy;
    double z0 = BioFVM::microenvironment.mesh.bounding_box[2] - BioFVM::microenvironment.mesh.dz;
    double x1 = BioFVM::microenvironment.mesh.bounding_box[3] + BioFVM::microenvironment.mesh.dx;
    double y1 = BioFVM::microenvironment.mesh.bounding_box[4] + BioFVM::microenvironment.mesh.dy;
    double z1 = BioFVM::microenvironment.mesh.bounding_box[5] + BioFVM::microenvironment.mesh.dz;

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

    Val target = PhysiCell::parameters.doubles("target_diff_3_N_cells");
    int i = 0;
    for (auto const& [_observable_domain, _effect_domain] : diff_domains) {
        auto cont1 = new Diff_Controller(
            // Observable and Effect Domain
            // Kernel::Iso_cuboid_3(Kernel::Point_3(-150.0, -150.0, -10.0), Kernel::Point_3(150.0, 150.0, 10.0)),
            _observable_domain,
            _effect_domain,
            // Target: How many differentiated cells do we want to have?
            target
        );
        supervisor.add_controller("Diff_Controller_diff_" + std::to_string(i), cont1);
        i++;
    }

    auto fp = fopen("controller_logs.txt", "w");
    // fprintf(stderr, "Prop: %E Diff: %E Int: %E Calc: %E sSize: %2d LastState: %2E\n", prop, differential, integral, calculated, state.size(), state.back());
    fprintf(fp, "StateSize, Proportional, Differential, Integral, Total, LastState\n");
    fclose(fp);

    // TODO CURRENTLY NOT WORKING
    // Add visitor to save state to csv file
    // supervisor.add_visitor_update("Safe_state", new Visitor_Write_State_to_csv("testfile.csv"));
    // supervisor.update_visit_all_controllers("Safe_state");

    // auto cb1 = Kernel::Iso_cuboid_3(Kernel::Point_3(-300.0, -300.0, -10.0), Kernel::Point_3(300.0, 300.0, 10.0));
    // auto cb2 = Kernel::Iso_cuboid_3(Kernel::Point_3(-100.0, 0.0, -10.0), Kernel::Point_3(250.0, 500.0, 10.0));

    std::cout << "**************************************\n" << std::endl;
}
