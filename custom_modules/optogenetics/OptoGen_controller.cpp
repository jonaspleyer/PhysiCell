// CGal libraries used for geometry handling
#include "OptoGen_controller.h"
#include <CGAL/enum.h>

namespace Opto::Controller {

void display_intersection()
{
    std::cout << "\n\n\n**************************************" << std::endl;
    std::cout << "          Setting up OptoGen          " << std::endl;
    std::cout << "**************************************\n" << std::endl;

    // Observable _observable<Kernel::Sphere_3>{};
    Observable<Kernel::Sphere_3, int> _observable{};
    Controller<Kernel::Sphere_3, int> cont(Kernel::Sphere_3(CGAL::ORIGIN, 2.0), _observable);

    // std::cout << cont.domain << std::endl;
    Kernel::Point_3 test = Kernel::Point_3(0.0,0.0,1.0);
    std::cout << test << std::endl;
    // std::cout << CGAL::ORIGIN << std::endl;
    std::cout << cont.domain << std::endl;
    std::cout << cont.is_point_in_domain(test) << std::endl;
}

template<class Domain, typename Value>
bool Controller<Domain, Value>::is_point_in_domain(Kernel::Point_3& point) {
    return domain.has_on_bounded_side(point);
}

template<class Domain, typename Value>
double Controller<Domain, Value>::get_volume() {
    return 0.0;
}

}