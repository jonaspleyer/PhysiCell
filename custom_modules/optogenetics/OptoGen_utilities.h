#ifndef OPTOGENETICS_UTILITIES
#define OPTOGENETICS_UTILITIES

#include "OptoGen_light.h"
#include "OptoGen_controller.h"

#include "../../core/PhysiCell_cell.h"

// Include CGal libraries for 2D and 3D geometry
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

#include <vector>
#include <algorithm>

namespace Opto::Utils {


// TODO this probably needs to be extended depending on which CGAL modules will be used in the future
template<class Domain>
bool point_is_in_domain(Domain& domain, std::vector<double> point) {
    return domain.has_on_bounded_side(Kernel::Point_3(point[0], point[1], point[2]));
}


template<class Domain>
std::vector<PhysiCell::Cell*> get_cells_in_domain(Domain& domain) {
    std::vector<PhysiCell::Cell*> cells{};
    
    std::copy_if(
        PhysiCell::all_cells->begin(),
        PhysiCell::all_cells->end(),
        std::back_inserter(cells),
        [&domain](PhysiCell::Cell* cell){
            return point_is_in_domain(domain, cell->position);
        }
    );
    
    return cells;
}


}

#endif