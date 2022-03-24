#include "./optogenetics/OptoGen.h"
#include <math.h>


void shine_light( const double& t);

void run_optogenetics ( const double &t );


typedef double Val;


class MyControllFunctor : public Opto::Controller::ControllFunctor {
	public:
		double adjust(double discrepancy);
};


class MyObservable : public Opto::Controller::Observable<Kernel::Sphere_3, Val> {
	public:
		Val calculate(Kernel::Sphere_3& _domain, std::vector<PhysiCell::Cell*> cells);
};


class MyMetric : public Opto::Controller::Metric<Val> {
	public:
		double calculate(Val& v1, Val& v2);
};


class MyEffect : public Opto::Controller::Effect {
    public:
        void apply(PhysiCell::Cell* cell, const double discrepancy);
};


void setup_optogenetics( void );