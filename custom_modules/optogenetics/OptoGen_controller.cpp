// CGal libraries used for geometry handling
#include "OptoGen_controller.h"
#include <CGAL/enum.h>

namespace Opto::Controller {


template<class Domain, typename Value>
bool Controller<Domain, Value>::is_point_in_domain(Kernel::Point_3& point) {
    return domain.has_on_bounded_side(point);
}

template<class Domain, typename Value>
double Controller<Domain, Value>::get_volume() {
    return 0.0;
}

}