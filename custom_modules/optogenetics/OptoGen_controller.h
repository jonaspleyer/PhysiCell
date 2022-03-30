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
template<class Domain, typename Value>
class Observable {
    public:
        std::deque<Value> state{};
        virtual Value measure(Domain& _domain, std::vector<PhysiCell::Cell*> cells) = 0;
};


// Used to meaningfully compare a observed to a target value
template<typename Value>
class Metric {
    public:
        std::deque<double> state{};
        virtual double calculate(Value& v1, Value& v2) = 0;
};


// Used in describing how regulation should work.
// TODO This should probably include information about the state in the future
class ControllFunctor {
    public:
        virtual double adjust(std::deque<double> state) = 0;
};


// Used to alter cell behaviour
class Effect {
    public:
        std::unique_ptr<Opto::Light::LightSource> lightsource{};
        virtual void apply(PhysiCell::Cell *cell, double const discrepancy) = 0;
};


// Only specify that class exists to already reference in some functions
template<class Domain, typename Value, class DerivedT>
class Controller;


// *********************************************************************************
// Template functions for actually running and updating controllers.
// These functions will be used in the visitor pattern outlined below
template<class Domain, typename Value, class DerivedT>
void run_single_controller(Controller<Domain, Value, DerivedT>& controller) {
    // Perform static cast to have access to derived class values
    auto controller_derived = static_cast<DerivedT*>(&controller);
    update_controller<Domain, Value, DerivedT>(*controller_derived);
    // Apply effect for every cell in domain
    std::for_each(
        controller_derived->cells_in_domain.begin(),
        controller_derived->cells_in_domain.end(),
        [&controller_derived](PhysiCell::Cell* cell){
            controller_derived->effect->apply(
                cell,
                controller_derived->controllfunctor->adjust(controller_derived->metric->state)
            );
        }
    );
}


template<class Domain, typename Value, class DerivedT>
void update_controller(Controller<Domain, Value, DerivedT>& controller) {
    auto controller_derived = static_cast<DerivedT*>(&controller);
    controller_derived->cells_in_domain = Opto::Utils::get_cells_in_domain(controller_derived->domain);
    controller_derived->update_state();
}


// *********************************************************************************
// Visitors for running (see steps 1-4 for controller)
// and updating (see steps 1)
class Visitor_Run {
public:
    template<class Domain, typename Value, class DerivedT>
    void visit(Controller<Domain, Value, DerivedT>& ci) {
        run_single_controller(ci);
    }
};


class Visitor_Update {
public:
    template<class Domain, typename Value, class DerivedT>
    void visit(Controller<Domain, Value, DerivedT>& ci) {
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
// 1. Measure                           ==> domain, observable
// 2. Compare to target                 ==> metric, target
// 3. Quantify control response         ==> controllfunctor
// 4. Determine optogenetic effect      ==> effect
//    on cells
template<class Domain, typename Value, class DerivedT>
class Controller : public ControllerElement {
private:
    int state_max_size = 500;
public:
    std::vector<PhysiCell::Cell*> cells_in_domain{};
    
    Domain domain{};
    std::unique_ptr<Observable<Domain, Value>> observable{};
    Value target{};
    std::unique_ptr<Metric<Value>> metric{};
    std::unique_ptr<ControllFunctor> controllfunctor{};
    std::unique_ptr<Effect> effect{};
    
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
        derived->observable->state.push_back(derived->observable->measure(derived->domain, derived->cells_in_domain));

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
            visitors_run_by_name["Default"] = std::make_unique<Visitor_Run>();
            visitors_update_by_name["Default"] = std::make_unique<Visitor_Update>();
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
        template<class Domain, typename Value, class DerivedT>
        void add_controller(std::string const& name, Controller<Domain, Value, DerivedT>* controller)
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