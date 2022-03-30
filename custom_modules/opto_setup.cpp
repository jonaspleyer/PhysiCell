#include "opto_setup.h"

double optogenetics_dt = 0.0;
double optogenetics_next_run_time = 0.0;
double optogenetics_update_dt = 0.0;
double optogenetics_next_update_time = 0.0;
double minimum_time_interval = 0.0;

// Create a supervisor to monitor all controllers
Opto::Controller::Supervisor supervisor;


void shine_light( const double& t) {
    // std::cout << "[Main] Shining light at time " << t << std::endl;
    supervisor.run_all_controllers();
    return;
}


void run_optogenetics ( const double& t ) {
    if (t-optogenetics_next_update_time + 0.01*minimum_time_interval > 0) {
        optogenetics_next_update_time += optogenetics_update_dt;
        // Additionally run opdates for optogenetic controllers
        // std::cout << "[Main] Running opto update at time " << t << std::endl;
        supervisor.update_all_controllers();
    } else if (t-optogenetics_next_run_time + 0.01*minimum_time_interval > 0) {
        optogenetics_next_run_time+=optogenetics_dt;
        shine_light(t);
    }
}


struct Red_LED : public Opto::Light::LightSource {
    double free_intensity = 0;
    double free_frequency = 500.0;
};


Val MyObservableSphere::measure(Kernel::Sphere_3& _domain, std::vector<PhysiCell::Cell*> cells) {
    return 1;
}


Val MyObservableTetrahedon::measure(Kernel::Tetrahedron_3& _domain, std::vector<PhysiCell::Cell*> cells) {
    return 1;
}


double MyMetric::calculate(Val& v1, Val& v2) {
    return fabs(v1-v2);
}


double MyControllFunctor::adjust(std::deque<double> state) {
    return state.back();
}


void MyEffect::apply(PhysiCell::Cell* cell, const double discrepancy) {
    return;
}


void setup_optogenetics( void ) {
    optogenetics_dt = PhysiCell::parameters.doubles("optogenetics_dt");
    optogenetics_update_dt = PhysiCell::parameters.doubles("optogenetics_udpate_dt");
    // TODO check if this works even when setting diffusion dt in xml file lower
    minimum_time_interval = std::min(std::min(optogenetics_dt, optogenetics_update_dt), PhysiCell::mechanics_dt);

	std::cout << "\n\n\n**************************************" << std::endl;
    std::cout << "          Setting up OptoGen          " << std::endl;
    std::cout << "**************************************\n" << std::endl;

    auto cont1 = new DensityController(
        Kernel::Sphere_3(CGAL::ORIGIN, 50.0),
        // Set density target to 10.0 cells/volume(Sphere)
        10.0
    );

    DensityController cont2(
        Kernel::Sphere_3(CGAL::ORIGIN, 50.0),
        12.0
    );
    
    supervisor.add_controller("MyFirstController", cont1);
    // supervisor.add_controller("MySecondController", cont2);

    supervisor.remove_controller("ThirdController");
    
    // supervisor.add_visitor_run("Default", new Opto::Controller::Visitor_Run);

    std::cout << "**************************************\n" << std::endl;
}