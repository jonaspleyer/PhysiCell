#ifndef OPTOGENETICS_CONTROLLER
#define OPTOGENETICS_CONTROLLER

// Custom light libraries used to define cell behaviour
#include "OptoGen_light.h"
#include "OptoGen_utilities.h"

// Include CGal libraries for 2D and 3D geometry
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;


// Include cell definitions for processing
#include "../../core/PhysiCell_cell.h"

namespace Opto::Controller {



template<class Domain, typename Value>
class Observable {
    public:
        virtual Value calculate(Domain& _domain, std::vector<PhysiCell::Cell*> cells) = 0;
};


// Used to meaningfully compare a observed to a target value
template<typename Value>
class Metric {
    public:
        virtual double calculate(Value& v1, Value& v2) = 0;
};


// Used in describing how regulation should work.
// TODO This should probably include information about the state in the future
class ControllFunctor {
    public:
        virtual double adjust(double discrepancy) = 0;
};


// Used to alter cell behaviour
class Effect {
    public:
        virtual void apply(PhysiCell::Cell* cell, const double discrepancy) = 0;
};


// Base class of all controllers
// Currently only used to store multiple controller instances in a single vector
struct BaseController {
};


template<class Domain, typename Value>
struct Controller : public BaseController {
    // TODO figure out a way to capture state
    std::vector<PhysiCell::Cell*> cells_in_domain{};
    // Give all necessary information upon instantiation
    Controller(
        Opto::Light::LightSource* _lightsource,
        Domain _domain,
        Observable<Domain, Value>* _observable,
        Value _target,
        Metric<Value>* _metric,
        ControllFunctor* _controllfunctor,
        Effect* _effect
    ) {
        this->lightsource = _lightsource;
        this->domain = _domain;
        this->observable = _observable;
        this->target = _target;
        this->metric = _metric;
        this->controllfunctor = _controllfunctor;
        this->effect = _effect;
    }
    // TODO include this in constructor
    Opto::Light::LightSource* lightsource;
    Domain domain;
    Observable<Domain, Value>* observable;
    Value target;
    Metric<Value>* metric;
    ControllFunctor* controllfunctor;
    Effect* effect;
};


template<class Domain, typename Value>
void run_controller(Controller<Domain, Value>& controller) {
    Value v{};
    Value observed = controller.observable->calculate(controller.domain, *PhysiCell::all_cells);
    double discrepancy = controller.metric->calculate(controller.target, observed);

    std::for_each(
        controller.cells_in_domain.begin(),
        controller.cells_in_domain.end(),
        [&controller, &discrepancy](PhysiCell::Cell* cell){controller.effect->apply(cell, controller.controllfunctor->adjust(discrepancy));}
    );
    // controller.effect->apply(new PhysiCell::Cell, discrepancy);
}


// TODO implement supervisor class for controller management
class Supervisor {
    private:
        std::map<std::unique_ptr<BaseController>, std::vector<PhysiCell::Cell*>> cells_by_controller;
        std::map<std::string, std::unique_ptr<BaseController>> controllers;
    public:
        void run_all_controllers();
        void distribute_cells_to_controllers(std::vector<PhysiCell::Cell*> cells);

        void add_controller(const std::string& name, BaseController& controller);
        void remove_controller(const std::string& name);
};


}

#endif