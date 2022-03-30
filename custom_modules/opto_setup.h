#include "./optogenetics/OptoGen.h"
#include "./optogenetics/OptoGen_controller.h"
#include <math.h>


void shine_light( const double& t);

void run_optogenetics ( const double &t );


typedef double Val;


class MyControllFunctor : public Opto::Controller::ControllFunctor {
	public:
		double adjust(std::deque<double> state);
};


class MyObservableSphere : public Opto::Controller::Observable<Kernel::Sphere_3, Val> {
	public:
		Val measure(Kernel::Sphere_3& _domain, std::vector<PhysiCell::Cell*> cells) override;
};


class MyObservableTetrahedon : public Opto::Controller::Observable<Kernel::Tetrahedron_3, Val> {
	public:
		Val measure(Kernel::Tetrahedron_3& _domain, std::vector<PhysiCell::Cell*> cells) override;
};


class MyMetric : public Opto::Controller::Metric<Val> {
	public:
		double calculate(Val& v1, Val& v2);
};


class MyEffect : public Opto::Controller::Effect {
    public:
        void apply(PhysiCell::Cell* cell, const double discrepancy);
};


class DensityController : public Opto::Controller::Controller<Kernel::Sphere_3, Val, DensityController> {
public:
	DensityController(
		Kernel::Sphere_3 _domain,
		Val _target
	) {
		domain = _domain;
		target = _target;
	}
    std::unique_ptr<Opto::Light::LightSource> lightsource = std::make_unique<Opto::Light::LightSource>();
    Kernel::Sphere_3 domain{};
    std::unique_ptr<MyObservableSphere> observable = std::make_unique<MyObservableSphere>();
    Val target{};
    std::unique_ptr<MyMetric> metric = std::make_unique<MyMetric>();
    std::unique_ptr<MyControllFunctor> controllfunctor = std::make_unique<MyControllFunctor>();
    std::unique_ptr<MyEffect> effect = std::make_unique<MyEffect>();
};


void setup_optogenetics( void );