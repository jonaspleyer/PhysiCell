#include "./optogenetics/OptoGen.h"
#include "./optogenetics/OptoGen_controller.h"
#include "../core/PhysiCell_utilities.h"
#include <math.h>

// *********************************************************************************
// RUN OPTOGENETIC ROUTINES IN MAIN
void shine_light( const double& t);

void run_optogenetics ( const double &t );


typedef double Val;


// *********************************************************************************
// LIGHT MODULES
struct Red_LED : public Opto::Light::LightSource {
    double free_intensity = 0;
    double free_frequency = 500.0;
};


// *********************************************************************************
// DENSITY CONTROLLER MODULES
class DensityControllFunctor : public Opto::Controller::ControllFunctor {
	public:
		double adjust(std::deque<double> state);
};


class PD_Controllfunctor : public Opto::Controller::ControllFunctor {
	public:
		// Proportional multiplication constant
		double K_p = 5.0;
		// Differential multiplication constant
		double K_d = 5.0;
		// Time constant between update steps
		double update_dt = PhysiCell::parameters.doubles("optogenetics_update_dt");
		double adjust(std::deque<double> state);
};


class DensityObservableCuboid : public Opto::Controller::Observable<Val, Kernel::Iso_cuboid_3> {
	public:
		Val measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) override;
		int id_x;
		int id_y;
};


class DensityMetric : public Opto::Controller::Metric<Val> {
	public:
		double calculate(Val& v1, Val& v2);
};


class DensityEffect : public Opto::Controller::Effect<Kernel::Iso_cuboid_3> {
    public:
		std::unique_ptr<Opto::Light::LightSource> lightsource = std::make_unique<Red_LED>();
		double current_light_intensity = 0.0;
        void apply(PhysiCell::Cell* cell, const double discrepancy);
};


class DensityController : public Opto::Controller::Controller<Val, DensityController, Kernel::Iso_cuboid_3, Kernel::Iso_cuboid_3> {
public:
	DensityController(
		Kernel::Iso_cuboid_3 _domain,
		Val _target,
		int id_x,
		int id_y
	) {
		observable->observable_domain = _domain;
		effect->effect_domain = _domain;
		target = _target;

		observable->id_x = id_x;
		observable->id_y = id_y;
	}

	Kernel::Iso_cuboid_3 domain;
	std::unique_ptr<DensityObservableCuboid> observable = std::make_unique<DensityObservableCuboid>();
    Val target{};
    std::unique_ptr<DensityMetric> metric = std::make_unique<DensityMetric>();
    std::unique_ptr<PD_Controllfunctor> controllfunctor = std::make_unique<PD_Controllfunctor>();
    std::unique_ptr<DensityEffect> effect = std::make_unique<DensityEffect>();
};


class LoggingObservableCuboid : public Opto::Controller::Observable<Val, Kernel::Iso_cuboid_3> {
	public:
		Val measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) override;
		int id_x;
		int id_y;
};



class LoggingEffect : public Opto::Controller::Effect<Kernel::Iso_cuboid_3> {
    public:
		std::unique_ptr<Opto::Light::LightSource> lightsource = std::make_unique<Red_LED>();
		double current_light_intensity = 0.0;
        void apply(PhysiCell::Cell* cell, const double discrepancy);
};


class LoggingController : public Opto::Controller::Controller<Val, LoggingController, Kernel::Iso_cuboid_3, Kernel::Iso_cuboid_3> {
public:
	LoggingController(
		Kernel::Iso_cuboid_3 _domain,
		Val _target,
		int id_x,
		int id_y
	) {
		observable->observable_domain = _domain;
		effect->effect_domain = _domain;
		target = _target;

		observable->id_x = id_x;
		observable->id_y = id_y;
	}

	Kernel::Iso_cuboid_3 domain;
	std::unique_ptr<LoggingObservableCuboid> observable = std::make_unique<LoggingObservableCuboid>();
    Val target{};
    std::unique_ptr<DensityMetric> metric = std::make_unique<DensityMetric>();
    std::unique_ptr<PD_Controllfunctor> controllfunctor = std::make_unique<PD_Controllfunctor>();
    std::unique_ptr<LoggingEffect> effect = std::make_unique<LoggingEffect>();
};


// *********************************************************************************
// NEW VISITOR FOR SAVING STATE OF CONTROLLER TO CSV FILE
// TODO CURRENTLY NOT WORKING
class Visitor_Write_State_to_csv : public Opto::Controller::Visitor_Update {
public:
    Visitor_Write_State_to_csv(std::string _csv_state_file_name) {
        csv_state_file_name = _csv_state_file_name;
    }
    std::string csv_state_file_name;

    template<typename ObservableType, class DerivedT, class ObservableDomain, class EffectDomain=ObservableDomain>
    void visit(Opto::Controller::Controller<ObservableType, DerivedT, ObservableDomain, EffectDomain>& ci) {
        std::cout << "Running custom update function" << std::endl;
    }
};


// *********************************************************************************
// SETUP OPTOGENETIC CONTROL FRAMEWORK (PRE-MAIN ROUTINE)
void setup_optogenetics( void );
