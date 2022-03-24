#include "opto_setup.h"

double optogenetics_dt = 0.0;
double optogenetics_next_dt = 0.0;
double optogenetics_update_dt = 0.0;
double optogenetics_next_update_time = 0.0;
double minimum_time_interval = 0.0;

// Create a supervisor to monitor all controllers
Opto::Controller::Supervisor supervisor;


void shine_light( const double& t) {
    std::cout << "Shining light at time " << t << std::endl;
    return;
}


void run_optogenetics ( const double& t ) {
    if (t-optogenetics_next_update_time + 0.01*minimum_time_interval > 0) {
        optogenetics_next_update_time += optogenetics_update_dt;
        // Additionally run opdates for optogenetic controllers
        std::cout << "Running opto update at time " << t << std::endl;
    }
    if (t-optogenetics_next_dt + 0.01*minimum_time_interval > 0) {
        optogenetics_next_dt+=optogenetics_dt;
        shine_light(t);
    }
}


struct Red_LED : public Opto::Light::LightSource {
    double free_intensity = 0;
    double free_frequency = 500.0;
};


Val MyObservable::calculate(Kernel::Sphere_3& _domain, std::vector<PhysiCell::Cell*> cells) {
    return 1;
}


double MyMetric::calculate(Val& v1, Val& v2) {
    return fabs(v1-v2);
}


double MyControllFunctor::adjust(double discrepancy) {
    return discrepancy;
}


void MyEffect::apply(PhysiCell::Cell* cell, const double discrepancy) {
    std::cout << "Adjusting cell with discrepancy " << discrepancy << std::endl;
    return;
}


void setup_optogenetics( void ) {
    optogenetics_dt = PhysiCell::parameters.doubles("optogenetics_dt");
    optogenetics_update_dt = PhysiCell::parameters.doubles("optogenetics_udpate_dt");
    minimum_time_interval = std::min(std::min(optogenetics_dt, optogenetics_update_dt), PhysiCell::diffusion_dt);

	std::cout << "\n\n\n**************************************" << std::endl;
    std::cout << "          Setting up OptoGen          " << std::endl;
    std::cout << "**************************************\n" << std::endl;
    
    MyMetric _metric;
    // Initialize the controller
    Opto::Controller::Controller<Kernel::Sphere_3, Val> cont(
        // Define which light to shine on that domain
        new Red_LED,
        // Define on which domain the controller will work
        Kernel::Sphere_3(CGAL::ORIGIN, 2.0),
        // Define a observable which we want to look at
        new MyObservable,
        // Now define a target to let the Controller aim for
        1.0,
        // Finally specify a metric to compare observable and target
        new MyMetric,
        // Define how to controll (IE integral, Differentiation, etc.)
        new MyControllFunctor,
        // Specify the effect (ie. the actual function of the controller)
        new MyEffect
    );

    supervisor.add_controller("MyFirstController", cont);

	// std::cout << cont.domain << std::endl;
    // Kernel::Point_3 test = Kernel::Point_3(0.0,0.0,1.0);
    /* std::cout << test << std::endl;
    // std::cout << CGAL::ORIGIN << std::endl;
    std::cout << cont.domain << std::endl;
    std::cout << cont.is_point_in_domain(test) << std::endl;*/
    std::cout << "\n" << std::endl;
}