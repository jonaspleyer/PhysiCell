// CGal libraries used for geometry handling
#include "OptoGen_controller.h"
#include <CGAL/enum.h>

namespace Opto::Controller {


// TODO implement this correctly
void Supervisor::add_controller(const std::string& name, BaseController& controller) {
    // cells_by_controller[controller] = std::vector<PhysiCell::Cell*>{};
    // controllers.push_back(controller);
    std::cout << "Added a controller with the name " << name << std::endl;
    return;
}


// TODO implement this
void Supervisor::remove_controller(const std::string& name) {

}


// TODO implement this
void Supervisor::run_all_controllers() {

}


// TODO implement this
void Supervisor::distribute_cells_to_controllers(std::vector<PhysiCell::Cell*> cells) {

}

}