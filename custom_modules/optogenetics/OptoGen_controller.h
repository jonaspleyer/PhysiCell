#ifndef OPTOGENETICS_CONTROLLER
#define OPTOGENETICS_CONTROLLER

// Custom light libraries used to define cell behaviour
#include "OptoGen_light.h"
#include "OptoGen_utilities.h"

// Include CGal libraries for 2D and 3D geometry
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Include modules for container definition to store template classes
#include <deque>

// Include cell definitions for processing
#include "../../core/PhysiCell_cell.h"

namespace Opto::Controller {


// Interface for defining observable quantities
template<typename ObservableType, class ObservableDomain>
class Observable {
    public:
        ObservableDomain observable_domain;
        std::vector<PhysiCell::Cell*> cells_in_observable_domain{};
        std::deque<ObservableType> state{};
        virtual ObservableType measure(ObservableDomain& _domain, std::vector<PhysiCell::Cell*> cells) = 0;
};


// Used to meaningfully compare a observed to a target value
template<typename ObservableType>
class Metric {
    public:
        std::deque<double> state{};
        virtual double calculate(ObservableType& target, ObservableType& observed) = 0;
};


// Used in describing how regulation should work.
// TODO This should probably include information about the state in the future
class ControllFunctor {
    public:
        virtual double adjust(std::deque<double> state) = 0;
};


// Used to alter cell behaviour
template<class EffectDomain>
class Effect {
    public:
        EffectDomain effect_domain;
        std::vector<PhysiCell::Cell*> cells_in_effect_domain{};
        std::unique_ptr<Opto::Light::LightSource> lightsource{};
        virtual void apply(PhysiCell::Cell *cell, double const discrepancy) = 0;
};


// Only specify that class exists to already reference in some functions
template<typename ObservableType, class DerivedT, class ObservableDomain, class EffectDomain=ObservableDomain>
class Controller;


// *********************************************************************************
// Template functions for actually running and updating controllers.
// These functions will be used in the visitor pattern outlined below
template<typename ObservableType, class DerivedT, class ObservableDomain, class EffectDomain=ObservableDomain>
void run_single_controller(Controller<ObservableType, DerivedT, ObservableDomain, EffectDomain>& controller) {
    // Perform static cast to have access to derived class values
    auto controller_derived = static_cast<DerivedT*>(&controller);
    update_controller<ObservableType, DerivedT, ObservableDomain, EffectDomain>(*controller_derived);
    
    // Calculate the discrepancy which is only dependent on the controller instance
    double discrepancy = controller_derived->controllfunctor->adjust(controller_derived->metric->state);
    
    // Apply effect for every cell in domain of effect
    std::for_each(
        controller_derived->effect->cells_in_effect_domain.begin(),
        controller_derived->effect->cells_in_effect_domain.end(),
        [&controller_derived, &discrepancy](PhysiCell::Cell* cell){
            controller_derived->effect->apply(
                cell,
                discrepancy
            );
        }
    );
}


template<typename ObservableType, class DerivedT, class ObservableDomain, class EffectDomain=ObservableDomain>
void update_controller(Controller<ObservableType, DerivedT, ObservableDomain, EffectDomain>& controller) {
    auto controller_derived = static_cast<DerivedT*>(&controller);
    controller_derived->observable->cells_in_observable_domain = Opto::Utils::get_cells_in_domain(controller_derived->observable->observable_domain);
    controller_derived->effect->cells_in_effect_domain = Opto::Utils::get_cells_in_domain(controller_derived->effect->effect_domain);
    controller_derived->update_state();
}


// *********************************************************************************
// Visitors for running (see steps 1-4 for controller)
// and updating (see steps 1)
class Visitor_Run {
public:
    template<typename ObservableType, class DerivedT, class ObservableDomain, class EffectDomain=ObservableDomain>
    void visit(Controller<ObservableType, DerivedT, ObservableDomain, EffectDomain>& ci) {
        run_single_controller(ci);
    }
};


class Visitor_Update {
public:
    template<typename ObservableType, class DerivedT, class ObservableDomain, class EffectDomain=ObservableDomain>
    void visit(Controller<ObservableType, DerivedT, ObservableDomain, EffectDomain>& ci) {
        update_controller(ci);
    }
};


class ControllerElement {
public:
    virtual void accept(Visitor_Run& v) = 0;
    virtual void accept(Visitor_Update& v) = 0;
};


// *********************************************************************************
// Compactly stores all necessary information related to
// 1. Measure                           ==> observable->observable_domain, observable
// 2. Compare to target                 ==> metric, target
// 3. Quantify control response         ==> controllfunctor
// 4. Apply optogenetic effect          ==> effect->effect_domain, effect
//    on cells
template<typename ObservableType, class DerivedT, class ObservableDomain, class EffectDomain>
class Controller : public ControllerElement {
private:
    int state_max_size = 500;
public:
    std::unique_ptr<Observable<ObservableDomain, ObservableType>> observable{};
    ObservableType target{};
    std::unique_ptr<Metric<ObservableType>> metric{};
    std::unique_ptr<ControllFunctor> controllfunctor{};
    std::unique_ptr<Effect<EffectDomain>> effect{};

    void accept(Visitor_Run& v) override {
        v.visit(*this);
    };
    
    void accept(Visitor_Update& v) override {
        v.visit(*this);
    };

    void update_state() {
        // Perform static cast to gain access to derived class member variables
        auto derived = static_cast<DerivedT*>(this);

        // Update Observable state
        if (derived->observable->state.size() >= derived->state_max_size) {
            derived->observable->state.erase(derived->observable->state.begin());
        }
        derived->observable->state.push_back(derived->observable->measure(derived->observable->observable_domain, derived->observable->cells_in_observable_domain));

        // Update metric state
        if (derived->metric->state.size() >= derived->state_max_size) {
            derived->metric->state.erase(derived->metric->state.begin());
        }
        derived->metric->state.push_back(derived->metric->calculate(derived->target, derived->observable->state.back()));
    }
};


// *********************************************************************************
// Used to store and manage controllers and visitors for updating and running controllers
class Supervisor {
    private:
        // Store Controllers here
        std::map<std::string, std::unique_ptr<ControllerElement>> controllers_by_name{};

        // Store visitors here
        std::map<std::string, std::unique_ptr<Visitor_Update>> visitors_update_by_name{};// {{"Default", new Visitor_Update}, };
        std::map<std::string, std::unique_ptr<Visitor_Run>> visitors_run_by_name{};// {{"Default", std::make_unique<Visitor_Run>()}};

        void raw_update_controller_with_visitor(ControllerElement& controller, Visitor_Update& v);
        void raw_run_controller_with_visitor(ControllerElement& controller, Visitor_Run& v);

        bool controller_name_is_present(std::string const& controller_name);
        bool visitor_run_name_is_present(std::string const& visitor_run_name);
        bool visitor_update_name_is_present(std::string const& visitor_update_name);
    public:
        Supervisor() {
            add_visitor_run("Default", new Visitor_Run);
            add_visitor_update("Default", new Visitor_Update);
        }

        // Update controllers
        void update_controller_with_visitor(std::string const& controller_name, std::string const& visitor_update_name);
        void update_controller_with_all_visitors(std::string const& controller_name);
        void update_visit_all_controllers(std::string const& visitor_update_name);
        void update_all_controllers();

        // Run controllers
        void run_controller_with_visitor(std::string const& controller_name, std::string const& visitor_run_name);
        void run_controller_with_all_visitors(std::string const& controller_name);
        void run_visit_all_controllers(std::string const& visitor_run_name);
        void run_all_controllers();
        
        // Add/Remove Controllers/Visitors
        template<typename ObservableType, class DerivedT, class ObservableDomain, class EffectDomain=ObservableDomain>
        void add_controller(std::string const& name, Controller<ObservableType, DerivedT, ObservableDomain, EffectDomain>* controller)
        {
            if (!controllers_by_name.contains(name)) {
                controllers_by_name[name] = std::unique_ptr<ControllerElement>(controller);
                std::cout << "[Opto] Added a controller with the name " << name << std::endl;
            } else {
                throw std::runtime_error("[Opto] Error: Controller with name " + name + " already present");
            }
        };
        void remove_controller(std::string const& name);
        void add_visitor_run(std::string const& visitor_run_name, Visitor_Run* v);
        void remove_visitor_run(std::string const& visitor_run_name);
        void add_visitor_update(std::string const& visitor_update_name, Visitor_Update* v);
        void remove_visitor_update(std::string const& visitor_update_name);
};


}

#endif