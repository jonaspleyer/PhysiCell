// CGal libraries used for geometry handling
#include "OptoGen_controller.h"

namespace Opto::Controller {


// *********************************************************************************
// SUPERVISOR
// *********************************************************************************
// Tests if controllers/visitors_(run/update) are present
bool Supervisor::controller_name_is_present(std::string const& controller_name) {
    return controllers_by_name.contains(controller_name);
}


bool Supervisor::visitor_run_name_is_present(std::string const& visitor_run_name) {
    return visitors_run_by_name.contains(visitor_run_name);
}


bool Supervisor::visitor_update_name_is_present(std::string const& visitor_name) {
    return visitors_update_by_name.contains(visitor_name);
}


// *********************************************************************************
// Methods to update controllers with visitors in different configurations
void Supervisor::update_controller_with_visitor(std::string const& controller_name, std::string const& visitor_update_name) {
    if (!controller_name_is_present(controller_name)) {
        throw std::runtime_error("[Opto] ERROR: Controller with name " + controller_name + " does not exist.");
    }
    if (!visitor_update_name_is_present(visitor_update_name)) {
        throw std::runtime_error("[Opto] ERROR: Visitor_Update with name " + visitor_update_name + " does not exist.");
    }
    controllers_by_name[controller_name]->accept(*visitors_update_by_name[visitor_update_name]);
}


void Supervisor::update_controller_with_all_visitors(std::string const& controller_name) {
    if (!controller_name_is_present(controller_name)) {
        throw std::runtime_error("[Opto] ERROR: Controller with name " + controller_name + " does not exist.");
    }
    for (auto const& [visitor_update_name, visitor_update] : visitors_update_by_name) {
        controllers_by_name[controller_name]->accept(*visitor_update);
    }
}


void Supervisor::update_visit_all_controllers(std::string const& visitor_update_name) {
    if (!visitor_update_name_is_present(visitor_update_name)) {
        throw std::runtime_error("[Opto] ERROR: Visitor_Update with name " + visitor_update_name + " does not exist.");
    }
    for (auto const& [controller_name, controller] : controllers_by_name) {
        controller->accept(*visitors_update_by_name[visitor_update_name]);
    }
}


void Supervisor::update_all_controllers() {
    for (auto const& [controller_name, controller] : controllers_by_name) {
        for (auto const& [visitor_update_name, visitor_update] : visitors_update_by_name) {
            controller->accept(*visitor_update);
        }
    }
}


// *********************************************************************************
// Methods to run controllers with visitors in different configurations 
void Supervisor::run_controller_with_visitor(std::string const& controller_name, std::string const& visitor_run_name) {
    if (!controller_name_is_present(controller_name)) {
        throw std::runtime_error("[Opto] ERROR: Controller with name " + controller_name + " does not exist.");
    }
    if (!visitor_run_name_is_present(visitor_run_name)) {
        throw std::runtime_error("[Opto] ERROR: Visitor_Run with name " + visitor_run_name + " does not exist.");
    }
    controllers_by_name[controller_name]->accept(*visitors_run_by_name[visitor_run_name]);
}


void Supervisor::run_controller_with_all_visitors(std::string const& controller_name) {
    if (!controller_name_is_present(controller_name)) {
        throw std::runtime_error("[Opto] ERROR: Controller with name " + controller_name + " does not exist.");
    }
    for (auto const& [visitor_run_name, visitor_run] : visitors_run_by_name) {
        controllers_by_name[controller_name]->accept(*visitor_run);
    }
}


void Supervisor::run_visit_all_controllers(std::string const& visitor_run_name) {
    if (!visitor_run_name_is_present(visitor_run_name)) {
        throw std::runtime_error("[Opto] ERROR: Visitor_Run with name " + visitor_run_name + " does not exist.");
    }
    for (auto const& [controller_name, controller] : controllers_by_name) {
        controller->accept(*visitors_run_by_name[visitor_run_name]);
    }
}


void Supervisor::run_all_controllers() {
    for (auto const& [controller_name, controller] : controllers_by_name) {
        for (auto const& [visitor_run_name, visitor_run] : visitors_run_by_name) {
            controller->accept(*visitor_run);
        }
    }
}


// *********************************************************************************
// Methods to manage visitors (ie. add and remove them)
void Supervisor::add_visitor_run(std::string const& visitor_run_name, Visitor_Run* v) {
    if (visitor_run_name_is_present(visitor_run_name)) {
        throw std::runtime_error("[Opto] ERROR: Visitor_Run with name " + visitor_run_name + " already exists.");
    }
    visitors_run_by_name[visitor_run_name] = std::unique_ptr<Visitor_Run>(v);
}


void Supervisor::remove_visitor_run(std::string const& visitor_run_name) {
    if (!visitor_run_name_is_present(visitor_run_name)) {
        std::cerr << "[Opto] WARNING: Visitor_Run with name " << visitor_run_name << " is not stored. Skip deleting Visitor_Run." << std::endl;
    }
    visitors_run_by_name.erase(visitor_run_name);
}


void Supervisor::add_visitor_update(std::string const& visitor_update_name, Visitor_Update* v) {
    if (visitor_update_name_is_present(visitor_update_name)) {
        throw std::runtime_error("[Opto] ERROR: Visitor_Update with name " + visitor_update_name + " already exists.");
    }
    visitors_update_by_name[visitor_update_name] = std::unique_ptr<Visitor_Update>(v);
}


void Supervisor::remove_visitor_update(std::string const& visitor_update_name) {
    if (!visitor_update_name_is_present(visitor_update_name)) {
        std::cerr << "[Opto] WARNING: Visitor_update with name " << visitor_update_name << " is not stored. Skip deleting Visitor_update." << std::endl;
    }
    visitors_update_by_name.erase(visitor_update_name);
}


// *********************************************************************************
// Methods to remove controller. Notice that a method for adding controllers was 
// implemented in header file as template.
void Supervisor::remove_controller(std::string const& controller_name) {
    if (!controller_name_is_present(controller_name)) {
        std::cerr << "[Opto] WARNING: Controller of name " << controller_name << " is not stored. Skip deleting Controller." << std::endl;
    }
    controllers_by_name.erase(controller_name);
}

}