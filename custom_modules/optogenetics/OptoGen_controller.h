#ifndef OPTOGENETICS_MODULE
#endif

// Custom light libraries used to define cell behaviour
#include "OptoGen_light.h"

// Include CGal libraries for 2D and 3D geometry
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;


namespace Opto::Controller {

void display_intersection();


template<class Domain, typename Value>
class Observable {
    private:

    public:
        int test = 0;
        // double operator()(Domain& _domain);
};


// TODO add custom classes with help of CGAL (possibly in seperate file)
template<class Domain, typename Value>
// template<class Domain>
class Controller {
    private:
        
    public:
        // We want to always construct from a instance of the domain object
        Controller(Domain _domain, Observable<Domain, Value> _observable) {
            domain = _domain;
            observable = _observable;
        }
        Domain domain{};

        bool is_point_in_domain(Kernel::Point_3& point);
        double get_volume();

        Observable<Domain, Value> observable;
        Value target;

        void run(const double& t);
};



}