#ifndef OPTOGENETICS_UTILITIES
#define OPTOGENETICS_UTILITIES

#include "OptoGen_light.h"
#include "OptoGen_controller.h"

#include "../../core/PhysiCell_cell.h"

#include <vector>
#include <algorithm>

namespace Opto::Utils {


// TODO this needs! to be fixed
template<class Domain>
bool point_is_in_domain(Domain& domain, std::vector<double> point) {
    return true;
}


// TODO check this implementation
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