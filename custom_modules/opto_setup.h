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
class PID_Controllfunctor : public Opto::Controller::ControllFunctor {
	public:
		// Proportional multiplication constant
		double K_p = PhysiCell::parameters.doubles("pid_controller_time_proportional");
		// Differential multiplication constant K_d = K_p * T_d
		// where T_d is a derivation time constant with which the controller will attempt to
		// approach the set point.
		double K_d = K_p * PhysiCell::parameters.doubles("pid_controller_time_differential");
		// Integral multiplication constant K_i = K_d  / T_i
		// where T_i is a integration time constant which describes how long the controller 
		// will consistently "permit" values higher/lower than expected.
		double K_i = K_p / PhysiCell::parameters.doubles("pid_controller_time_integral");
		// Time constant between update steps
		double update_dt = PhysiCell::parameters.doubles("optogenetics_update_dt");
		double adjust(std::deque<double> state);
		std::deque<std::array<double, 3>> calculated_responses;
		int calculated_responses_size = 200;
		double low_pass_cutoff = 0.02/PhysiCell::parameters.doubles("optogenetics_update_dt");
		int diff_index = -1;
};


class Diff_ObservableCuboid : public Opto::Controller::Observable<Val, Kernel::Iso_cuboid_3> {
	public:
		double diff_disable = PhysiCell::parameters.doubles("diff_disable");
		double diff_enable = 0.0;
		int diff_index = 0;
		Kernel::Iso_cuboid_3 domain;
		Val measure(Kernel::Iso_cuboid_3& _domain, std::vector<PhysiCell::Cell*> cells) override;
};


class Diff_Metric : public Opto::Controller::Metric<Val> {
	public:
		double calculate(Val& target, Val& observed);
};


class Diff_Effect : public Opto::Controller::Effect<Kernel::Iso_cuboid_3> {
    public:
		double diff_disable = PhysiCell::parameters.doubles("diff_disable");
		double diff_enable = 0.0;
		int diff_index = 0;
		std::unique_ptr<Opto::Light::LightSource> lightsource = std::make_unique<Red_LED>();
        void apply(PhysiCell::Cell* cell, const double discrepancy);
};


class Diff_Controller : public Opto::Controller::Controller<Val, Diff_Controller, Kernel::Iso_cuboid_3, Kernel::Iso_cuboid_3> {
public:
	Diff_Controller(
		Kernel::Iso_cuboid_3 _observable_cuboid,
		Kernel::Iso_cuboid_3 _effect_cuboid,
		Val _target,
		double diff_enable,
		int diff_index
	) {
		observable->observable_domain = _observable_cuboid;
		effect->effect_domain = _effect_cuboid;
		target = _target;
		effect->diff_enable = diff_enable;
		effect->diff_index = diff_index;
		observable->diff_enable = diff_enable;
		observable->diff_index = diff_index;
		controllfunctor->diff_index = diff_index;
	}

	int state_max_size = 20000;
	Kernel::Iso_cuboid_3 domain;
	std::unique_ptr<Diff_ObservableCuboid> observable = std::make_unique<Diff_ObservableCuboid>();
    Val target{};
    std::unique_ptr<Diff_Metric> metric = std::make_unique<Diff_Metric>();
    std::unique_ptr<PID_Controllfunctor> controllfunctor = std::make_unique<PID_Controllfunctor>();
    std::unique_ptr<Diff_Effect> effect = std::make_unique<Diff_Effect>();
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