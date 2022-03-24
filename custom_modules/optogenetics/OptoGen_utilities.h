#ifndef OPTOGENETICS_UTILITIES
#define OPTOGENETICS_UTILITIES

#include "OptoGen_light.h"
#include "OptoGen_controller.h"

#include "../../core/PhysiCell_cell.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

#include <vector>
#include <algorithm>

namespace Opto::Utils {


// TODO implement this in correct way
template<class Domain>
std::vector<PhysiCell::Cell*> get_cells_in_domain(Domain& domain) {
    std::vector<PhysiCell::Cell*> cells{};

    std::copy_if(
        PhysiCell::all_cells->begin(),
        PhysiCell::all_cells->end(),
        std::back_inserter(cells),
        [&domain](PhysiCell::Cell* cell){
            // TODO Here is the problem
            return true;
        }
    );
    
    return cells;
}


// TODO this needs! to be fixed
template<class Domain>
bool point_is_in_domain(Domain& domain, Kernel::Point_3& point) {
    return true;
}

}

#endif