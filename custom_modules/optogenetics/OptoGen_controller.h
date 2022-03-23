#ifndef OPTOGENETICS_CONTROLLER
#define OPTOGENETICS_CONTROLLER

// Custom light libraries used to define cell behaviour
#include "OptoGen_light.h"

// Include CGal libraries for 2D and 3D geometry
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;


namespace Opto::Controller {


template<class Domain, typename Value>
class Observable {
    private:

    public:
        int test = 0;
        // double operator()(Domain& _domain);
};


// Used to meaningfully compare a observed to a target value
template<typename Value>
class Metric {
    public:
        virtual double operator()(Value& v1, Value& v2);
};


// TODO add custom classes with help of CGAL (possibly in seperate file)
template<class Domain, typename Value>
// template<class Domain>
class Controller {
    private:
        
    public:
        // We want to always construct from a instance of the domain object
        Controller(Domain _domain, Observable<Domain, Value> _observable, Value _target, Metric<Value> _metric) {
            domain = _domain;
            observable = _observable;
            target = _target;
            metric = _metric;
        }
        Domain domain{};

        Opto::Light::LightSource lightsource{};

        bool is_point_in_domain(Kernel::Point_3& point);
        double get_volume();

        Observable<Domain, Value> observable;
        Value target;
        Metric<Value> metric;

        void run(const double& t);
};



}

#endif