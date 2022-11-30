#include "opto_setup.h"
#include <cmath>
#include <random>
#include <chrono>
#include <fstream>
#include <iostream>

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
// DENSITY CONTROLLER MODULES
Val DensityObservableCuboid::measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) {
    int N_cells = std::count_if(cells.begin(), cells.end(), [](PhysiCell::Cell* cell){return !cell->phenotype.death.dead;});

    // Write this to a file for plots later
    std::ofstream outfile;
    outfile.open("information.csv", std::ios_base::app);
    outfile << id_x << "," << id_y << "," << N_cells << "," << cells.size() << "," << PhysiCell::parameters.doubles("cell_target_boxes") << "\n";

    return N_cells;
}


double DensityMetric::calculate(Val& target, Val& observed) {
    return observed-target;
}


double DensityControllFunctor::adjust(std::deque<double> state) {
    // Target is implicitly alays 0.0
    return state.back();
}


double PD_Controllfunctor::adjust(std::deque<double> state) {
    // Target is implicitly always 0.0
    // This implements a PI Controller
    double calculated = K_p*state.back() + K_d*(state.back()-state[state.size()-1])/update_dt;
    return calculated;
}


void DensityEffect::apply(PhysiCell::Cell* cell, const double discrepancy) {
    double dt = PhysiCell::parameters.doubles("optogenetics_dt");

    double increase = std::min(
        PhysiCell::parameters.doubles("light_ion_exposure_increment") * discrepancy,
        PhysiCell::parameters.doubles("light_ion_exposure_increment_max")
    );

    // Increase the ion cocnentration by the amount of light exposure which is given
	cell->custom_data["light_ion_concentration"] += dt * increase;
    return;
}



// *********************************************************************************
// LOGGING CONTROLLER MODULES
Val LoggingObservableCuboid::measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) {
    int N_cells = std::count_if(cells.begin(), cells.end(), [](PhysiCell::Cell* cell){return !cell->phenotype.death.dead;});

    // Write this to a file for plots later
    std::ofstream outfile;
    outfile.open("information.csv", std::ios_base::app);
    outfile << id_x << "," << id_y << "," << N_cells << "," << cells.size() << "," << "NaN" << "\n";

    return N_cells;
}


void LoggingEffect::apply(PhysiCell::Cell* cell, const double discrepancy) { return; }


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


    // Create boxes with individual controllers
    int n_boxes = PhysiCell::parameters.ints("n_boxes_along_axis");

    // Get the size of the total simulation domain
    double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = microenvironment.mesh.bounding_box[3];
	double Ymax = microenvironment.mesh.bounding_box[4];
	double Zmax = microenvironment.mesh.bounding_box[5];

    double dX = (Xmax-Xmin)/n_boxes/2.0;
    double dY = (Ymax-Ymin)/n_boxes/2.0;
    
    for (int i=0; i<2*n_boxes; i++) {
        for (int j=0; j<2*n_boxes; j++) {
            // Define the box
            auto cuboid = Kernel::Iso_cuboid_3(
                Kernel::Point_3(Xmin + dX * i,      Ymin + dY * j,      Zmin),
                Kernel::Point_3(Xmin + dX * i + dX, Ymin + dY * j + dY, Zmax)
            );
            
            // If we are at an even vertical position, we take the even boxes horizontally
            if ((i+j) % 2==1) {
                // Define the controller for this box
                auto cont = new DensityController(
                    cuboid,
                    PhysiCell::parameters.doubles("cell_target_boxes"),
                    i,
                    j
                );

                // Controller is named by location on grid
                std::string name = "Density_controller_" + std::to_string(i) + "_" + std::to_string(j);

                // Add the controller to the supervisor such that is automatically run
                supervisor.add_controller(name, cont);
            }
            // Otherwise we simply use a dummy controller for logging purposes
            else {
                // Define the controller for this box
                auto cont = new LoggingController(cuboid, 0, i, j);

                // Controller is named by location on grid
                std::string name = "Logging_controller_" + std::to_string(i) + "_" + std::to_string(j);

                // Add the controller to the supervisor such that is automatically run
                supervisor.add_controller(name, cont);
            }
        }
    }

    std::ofstream outfile;
    outfile.open("information.csv");
    outfile << "id_x,id_y,N_cells,total cells,setpoint\n";

    std::cout << "**************************************\n" << std::endl;
}
