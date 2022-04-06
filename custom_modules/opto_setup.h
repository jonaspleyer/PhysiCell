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
// DIFFERENTIATION CONTROLLER MODULES
class PI_Controllfunctor : public Opto::Controller::ControllFunctor {
	public:
		// Proportional multiplication constant
		double K_p = 0.1;
		// Differential multiplication constant
		double K_i = 0.1;
		// Time constant between update steps
		double update_dt = PhysiCell::parameters.doubles("optogenetics_update_dt");
		double adjust(std::deque<double> state);
};


class Diff_1_ObservableCuboid : public Opto::Controller::Observable<Val, Kernel::Iso_cuboid_3> {
	public:
		double diff_disable = PhysiCell::parameters.doubles("diff_disable");
		double diff_enable_1 = PhysiCell::parameters.doubles("diff_enable_1");
		Kernel::Iso_cuboid_3 domain;
		Val measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) override;
};


class Diff_2_ObservableCuboid : public Opto::Controller::Observable<Val, Kernel::Iso_cuboid_3> {
	public:
		double diff_disable = PhysiCell::parameters.doubles("diff_disable");
		double diff_enable_2 = PhysiCell::parameters.doubles("diff_enable_2");
		Kernel::Iso_cuboid_3 domain;
		Val measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) override;
};


class Diff_3_ObservableCuboid : public Opto::Controller::Observable<Val, Kernel::Iso_cuboid_3> {
	public:
		double diff_disable = PhysiCell::parameters.doubles("diff_disable");
		double diff_enable_3 = PhysiCell::parameters.doubles("diff_enable_3");
		Kernel::Iso_cuboid_3 domain;
		Val measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) override;
};


class Diff_Metric : public Opto::Controller::Metric<Val> {
	public:
		double calculate(Val& v1, Val& v2);
};


class Diff_1_Effect : public Opto::Controller::Effect<Kernel::Iso_cuboid_3> {
    public:
		double diff_disable = PhysiCell::parameters.doubles("diff_disable");
		double diff_enable_1 = PhysiCell::parameters.doubles("diff_enable_1");
		std::unique_ptr<Opto::Light::LightSource> lightsource = std::make_unique<Red_LED>();
        void apply(PhysiCell::Cell* cell, const double discrepancy);
};


class Diff_2_Effect : public Opto::Controller::Effect<Kernel::Iso_cuboid_3> {
    public:
		double diff_disable = PhysiCell::parameters.doubles("diff_disable");
		double diff_enable_2 = PhysiCell::parameters.doubles("diff_enable_2");
		std::unique_ptr<Opto::Light::LightSource> lightsource = std::make_unique<Red_LED>();
        void apply(PhysiCell::Cell* cell, const double discrepancy);
};


class Diff_3_Effect : public Opto::Controller::Effect<Kernel::Iso_cuboid_3> {
    public:
		double diff_disable = PhysiCell::parameters.doubles("diff_disable");
		double diff_enable_3 = PhysiCell::parameters.doubles("diff_enable_3");
		std::unique_ptr<Opto::Light::LightSource> lightsource = std::make_unique<Red_LED>();
        void apply(PhysiCell::Cell* cell, const double discrepancy);
};


class Diff_1_Controller : public Opto::Controller::Controller<Val, Diff_1_Controller, Kernel::Iso_cuboid_3, Kernel::Iso_cuboid_3> {
public:
	Diff_1_Controller(
		Kernel::Iso_cuboid_3 _observable_cuboid,
		Kernel::Iso_cuboid_3 _effect_cuboid,
		Val _target
	) {
		observable->observable_domain = _observable_cuboid;
		effect->effect_domain = _effect_cuboid;
		target = _target;
	}

	Kernel::Iso_cuboid_3 domain;
	std::unique_ptr<Diff_1_ObservableCuboid> observable = std::make_unique<Diff_1_ObservableCuboid>();
    Val target{};
    std::unique_ptr<Diff_Metric> metric = std::make_unique<Diff_Metric>();
    std::unique_ptr<PI_Controllfunctor> controllfunctor = std::make_unique<PI_Controllfunctor>();
    std::unique_ptr<Diff_1_Effect> effect = std::make_unique<Diff_1_Effect>();
};


class Diff_2_Controller : public Opto::Controller::Controller<Val, Diff_2_Controller, Kernel::Iso_cuboid_3, Kernel::Iso_cuboid_3> {
public:
	Diff_2_Controller(
		Kernel::Iso_cuboid_3 _observable_cuboid,
		Kernel::Iso_cuboid_3 _effect_cuboid,
		Val _target
	) {
		observable->observable_domain = _observable_cuboid;
		effect->effect_domain = _effect_cuboid;
		target = _target;
	}

	Kernel::Iso_cuboid_3 domain;
	std::unique_ptr<Diff_2_ObservableCuboid> observable = std::make_unique<Diff_2_ObservableCuboid>();
    Val target{};
    std::unique_ptr<Diff_Metric> metric = std::make_unique<Diff_Metric>();
    std::unique_ptr<PI_Controllfunctor> controllfunctor = std::make_unique<PI_Controllfunctor>();
    std::unique_ptr<Diff_2_Effect> effect = std::make_unique<Diff_2_Effect>();
};


class Diff_3_Controller : public Opto::Controller::Controller<Val, Diff_3_Controller, Kernel::Iso_cuboid_3, Kernel::Iso_cuboid_3> {
public:
	Diff_3_Controller(
		Kernel::Iso_cuboid_3 _observable_cuboid,
		Kernel::Iso_cuboid_3 _effect_cuboid,
		Val _target
	) {
		observable->observable_domain = _observable_cuboid;
		effect->effect_domain = _effect_cuboid;
		target = _target;
	}

	Kernel::Iso_cuboid_3 domain;
	std::unique_ptr<Diff_3_ObservableCuboid> observable = std::make_unique<Diff_3_ObservableCuboid>();
    Val target{};
    std::unique_ptr<Diff_Metric> metric = std::make_unique<Diff_Metric>();
    std::unique_ptr<PI_Controllfunctor> controllfunctor = std::make_unique<PI_Controllfunctor>();
    std::unique_ptr<Diff_3_Effect> effect = std::make_unique<Diff_3_Effect>();
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