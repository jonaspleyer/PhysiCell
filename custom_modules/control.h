/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
/*
 * control.h
 *
 *  Created on: 17.01.2020
 *  Updated on: 21.08.2021
 *  Author: cfleck
 */

#ifndef CUSTOM_MODULES_CONTROL_H_
#define CUSTOM_MODULES_CONTROL_H_

#include "../core/PhysiCell.h"
#include "domain.h"
#include "light.h"

#include <vector>
#include<cmath>

// For debugging
#ifndef _HI_
#define _HI_ std::cout << "HERE I AM" << std::endl;
#endif

using namespace PhysiCell;

struct Density
{
	double time;
	double density;
	int number_of_cells;
	double area;
	std::string domain_type;
	std::string domain_name;
	Point domain_center;
	std::vector<Cell*> cells;
};

class Controller_Base;

// Virtual base class for the controller functor
class ControlFunctor_Base
{
protected:
	double t;
	double prev_t;
	double err;
	double prev_err;
public:
	bool erase_light_pattern = true;
	ControlFunctor_Base(): t{0}, prev_t{0}, err{0}, prev_err{0} {;}
	virtual ~ControlFunctor_Base() {;}
	virtual void operator()( Controller_Base* controller, double _err, double t ) = 0;
};

class P_Control : public ControlFunctor_Base
{
private:
	double cn;
public:
	P_Control() : Kp{0}, cn{0} {;}
	double Kp;
	void operator()( Controller_Base* controller, double _err, double _t );
};

class PI_Control : public ControlFunctor_Base
{
protected:
	double cn;
	double cn_i;
	double int_cn;

public:
	bool integration = false;
	PI_Control() : Kp{0}, Ki{0}, cn{0}, cn_i{0}, int_cn{0} {;}
	double Kp;
	double Ki;
	void switch_integration_on( void ){ integration = true; }
	void switch_integration_off( void ){ integration = false; int_cn = 0; }
	void operator()( Controller_Base* controller, double _err, double _t ) = 0;
};

class PI_Control_1 : public PI_Control
{
public:
	PI_Control_1() {;}
	void operator()( Controller_Base* controller, double _err, double _t );
};

class PI_Control_2 : public PI_Control
{
public:
	PI_Control_2() {;}
	void operator()( Controller_Base* controller, double _err, double _t );
};

//class Antithetic_Control : public ControllFunctor_Base
//{
//private:
//	double cn;
//public:
//	Antithetic_Control() : Kp{0}, Ki{0}, cn{0} {;}
//	double Kp;
//	double Ki;
//	void operator()( Controller_Base* controller, double _err, double _t );
//};

class Antithetic_Control : public ControlFunctor_Base
{
private:
	int Z1;
	int Z2;
	double cn;
public:
	double mu;
	double theta;
	double eta;

	Antithetic_Control() : mu{1.0}, theta{1.0}, eta{0.1}, Z1{0}, Z2{0}, cn{0} {;}
	void operator()( Controller_Base* controller, double _err, double _t );
};

class Controller_Base
{
protected:
	bool set = false;
	Domain* domain = NULL;
public:
	std::vector<Light*> lights;
	ControlFunctor_Base* control_functor = NULL;

	bool on = false;
	Controller_Base() : lights(0) {;}
	virtual ~Controller_Base() {;}
	virtual void set_up( std::vector<Light*> _lights, ControlFunctor_Base* _controll_function ) = 0;
	virtual void set_domain( Domain* _domain ) { domain = _domain; }
	virtual Point get_center( void ) { return domain->get_center(); }
	virtual void run( double t ) = 0;
};

// Determines the cell density in a predefined domain
class DensityFunctor
{
private:
	Point P;
	Density current_density;
	Domain* domain = NULL;
	Point center;
public:
	std::vector<Density> history;
	DensityFunctor() : history(0) {;}
	void set_domain( Domain* _domain ) { domain = _domain; center = (*_domain).get_center(); }
	std::string get_domain_type ( void ) { return (*domain).get_type(); }
	std::string get_domain_name ( void ) { return (*domain).get_name(); }
	Point get_center( void ) { return center; }
	Density operator()( void );
	Density operator()( int cell_type );
};

class Local_DensityFunctor
{
private:
	bool grid_set = false;
	bool measurement_set = false;

	// it is expected that the measurement domain is a subset of the grid domain
	Rectangle* grid_domain = NULL; // domain on which the underlying grid is created
	Point grid_domain_center;
	std::string grid_domain_type = "Rectangle";
	std::string grid_domain_name = "NoName";

	std::vector<Rectangle*> density_domains;
	int number_of_density_domains = 0;
	std::vector<DensityFunctor> density_functors;

	Domain* measurement_domain = NULL; // domain on which to measure the local densities
	std::string measurement_domain_type = "NoType";
	std::string measurement_domain_name = "NoName";
	double measurement_domain_area = 0;

	std::vector<DensityFunctor*> measure_density_functors;
	std::vector<Density> measured_densities;

	int total_cell_number = 0;
	double mean = 0;
	double variance = 0;
public:
	double control_density = 0;

	Local_DensityFunctor() : density_functors(0), density_domains(0), measure_density_functors(0),
			measured_densities(0) {;}

	bool set_grid_domain( Rectangle* _domain, int n );
	Point get_grid_domain_center( void ) { return grid_domain_center; }

	bool set_measurement_domain( Domain* _domain, bool invert = false );
	std::string get_measurement_domain_type( void ) { return (*measurement_domain).get_type(); }
	std::string get_measurement_domain_name( void ) { return measurement_domain_name; }
	double get_measurement_domain_area( void ) { return measurement_domain_area; }
	bool point_inside_measurement_domain( Point P ){ return measurement_domain ? (*measurement_domain)(P) : false; }

	int get_total_cell_number( void ){ return total_cell_number; }
	double get_mean( void ){ return mean; }
	double get_variance( void ){ return variance; }

	std::vector<Density> operator()( void );
	std::vector<Density> operator()( int cell_type_index );
};

// The project specific controller
class Density_Controller : public Controller_Base
{
private:
	int cell_type_index;
	std::string cell_type_name = "NoName";
	bool cell_type_set = false;

	bool PI_control = false;
public:
	DensityFunctor density;
	Density current_density;

	double control_density; // target density
	double error;

	Density_Controller() : cell_type_index{-1}, control_density{0.001}, error{0.0} {;};
	void set_up( std::vector<Light*> _lights, ControlFunctor_Base* _control_functor);
	void set_domain( Domain* _domain ) override { domain = _domain; density.set_domain( _domain ); }
	void set_cell_type( int _cell_type_index, std::string _cell_type_name = "NoName" );
	int get_cell_type_index( void ){ return cell_type_set ? cell_type_index : -1; }
	std::string get_cell_type_name( void ){ return cell_type_set ? cell_type_name : "NoName"; }
	void unset_cell_type( void ){ cell_type_name = "NoName"; cell_type_set = false; }

	void run( double t );
};

struct Control_Lattice
{
public:
	Control_Lattice() : LED_types(0), controller_param(0), density_domains(0), led_lattices(0), controller(0),
	cell_type_index(0), cell_type_name(0) {;}
	Rectangle* domain = NULL; // domain on which the lattice is to be build
	int number_of_density_domains = 0;
	int number_of_sides = 0;

	bool cell_types_defined = false;
	int number_of_cell_types = 1;
	std::vector<int> cell_type_index;
	std::vector<std::string> cell_type_name;
	std::vector<double> cell_type_density;

	std::vector<Light_Type> LED_types;
	std::vector<double> controller_param;

	std::vector<Rectangle*> density_domains;
	std::vector<LED_Lattice*> led_lattices;
	std::vector<Controller_Base*> controller;
};

//ToDo: implement proper error handling
class Light_Control
{
private:
	bool is_there_light = false;
	bool state_changed = false;
	double last_update_time;
	int number_of_lights;
	int number_of_controllers;

	int counter = 0;
public:
	std::vector<Local_DensityFunctor> density;
	std::vector<Density> local_densities;

	std::vector<Light*> lights;
	std::vector<Controller_Base*> controllers;
	double update_dt = 10.0;
	double update_tolerance = 0.0001;

	Light_Control() : number_of_lights{0}, lights(0), number_of_controllers{0}, controllers(0),
			last_update_time{0.0}, density(0), local_densities(0) { };

	void add_light( Light* light ) { lights.push_back( light ); ++number_of_lights; }
	void remove_light( int n ) { lights.erase( lights.begin() + n ); --number_of_lights; }
	int get_number_of_lights( void ) { return number_of_lights; }

	void add_controller( Controller_Base* controller );
	void remove_controller( int n );
	Light& operator[]( int n ) { return *lights[n]; }

	// check whether there is light
	bool operator()( Point x );
	bool operator()( Point x, int m );
	// check whether there is light of LED type m
	bool operator()( Point x, double t );
	bool operator()( Point x, double t, int m );
	// get light parameters
	double get_density( Point x, double w, int m = 0 );
	double get_density( Point x, double t, double w, int m = 0 );
	double get_intensity( Point x, int m = 0 );
	double get_intensity( Point x, double t, int m = 0 );
	// set the intensity of the LEDs at point x
	void set_intensity( double intensity, Point x, int m = -1 );

	void update ( double t );
// TODO: change the return of the () operator; instead of just true/false return the indices of the lights switch on
};
#endif /* CUSTOM_MODULES_CONTROL_H_ */
