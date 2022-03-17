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
 * optogenetics.cpp
 *
 *  Created on: 28.10.2019
 *  Updated on: 21.08.2021
 *  Author: Christian Fleck
 *  Institution: ETH Zürich
 *
 *  Project: CyGenTiG
 *
 *
 *  Comment: Here are the functionality to represent optogenetic manipulation on the cellular tissue.
 *  Note that we employ a count-clockwise counting for the nodes.
 *
 *  ToDo:Error Handling; Time dependent domains;
 *  Optogenetic class encompassing domain and light specification (intensity, spectrum, etc.)
 */

#include  <cmath>
#include "./optogenetics.h"

// The governing class for optogenetics
Light_Control opto;
Light_Type blue_led("BLUE LED", {500,10});
Light_Type red_led("RED LED", {660,10});
Light_Type green_led("GREEN LED", {500,10});
int index_red;
int index_blue;
int index_green;

/*
This function creates a control lattice. A control lattice consists of nxn density domains on which the
density is to be measured and controlled. On each of these density domains a LED lattice is defined,
consisting of m LEDs of given types.
*/
Control_Lattice create_control_lattice( Control_Lattice control_lattice )
{
	int n = control_lattice.number_of_density_domains;
	int m = control_lattice.number_of_sides;

	// check some things first
	if(!control_lattice.domain)
	{
		std::cout << "Create control lattice: No domain given." << std::endl;
		return control_lattice;
	}
	if(control_lattice.LED_types.size() == 0)
	{
		std::cout << "Create control lattice: No LEDs defined." << std::endl;
		return control_lattice;
	}
	std::cout << "Creating a density control lattice on " << control_lattice.domain->get_name() << " with " << n*n << " density domains." << std::endl;
	std::cout << "Area of each domain is " << control_lattice.domain->get_area()/(n*n)  << " um^2." << std::endl;
	std::cout << "Each density domain has " << m*m << " lattices sides." << std::endl;
	std::cout << "On each lattice side are " << control_lattice.LED_types.size() << " LED(s) defined." << std::endl;

	// create the density lattice
	control_lattice.density_domains.resize(n*n);
	for( int i = 0; i < n*n; ++i )
	{
		control_lattice.density_domains[i] = new Rectangle;
	}
	Point P0, P1;
	Rectangle primitive_domain;
	Rectangle* wd = control_lattice.domain;

	std::vector<Point> nodes(2);
	std::vector<double> translation(2);
	std::vector<double> d1(2);
	std::vector<double> d2(2);
	d1 = {wd->AB[0]/n, wd->AB[1]/n};
	d2 = {wd->AD[0]/n, wd->AD[1]/n};
//	std::cout << d1[0] << " " << d1[1] << std::endl;
//	std::cout << d2[0] << " " << d2[1] << std::endl;

	// construct the primitive lattice side
	P0 = wd->get_nodes()[0];
	P1 = P0; P1.move(d1);
//	std::cout << P1.tuple() << std::endl;
	nodes = {P0,P1}; primitive_domain.set_up(nodes, sqrt(wd->AD_AD)/n);

	auto iter = control_lattice.density_domains.begin();
	for( int i = 0; i < n; ++i )
	{
		for( int j = 0; j < n; ++j )
		{
			*(*iter) = primitive_domain;
			translation = j*d2 + i*d1;
//			std::cout << translation[0] << " " << translation[1] << std::endl;
			(*iter)->move( translation );
//			std::cout << (*iter)->get_center().tuple() << std::endl;
			(*iter)->set_name(std::to_string(i*n+j));
//			std::cout << "translation vector = " << translation << std::endl;
//			std::cout << (*iter)->get_nodes()[0].tuple() << std::endl;
//			std::cout << (*iter)->get_nodes()[1].tuple() << std::endl;
			++iter;
		}
	}
	// create the LED lattices
	control_lattice.led_lattices.resize(n*n);
	for( int i = 0; i < n*n; ++i )
	{
		control_lattice.led_lattices[i] = new LED_Lattice;
		control_lattice.led_lattices[i]->set_domain(control_lattice.density_domains[i]);
		control_lattice.led_lattices[i]->set_lattice(m);
		for( auto iter = control_lattice.LED_types.begin(); iter < control_lattice.LED_types.end(); ++iter )
		{
			control_lattice.led_lattices[i]->add_LED(*iter);
		}
		// switches the lattice on
		control_lattice.led_lattices[i]->switch_on();
	}
	// density controller
	control_lattice.controller.resize(control_lattice.number_of_cell_types*n*n);

	Density_Controller* pC;
	PI_Control* pCF;
//	Antithetic_Control* pCF;
	bool erase;
	auto ll = control_lattice.led_lattices.begin();
	auto dc = control_lattice.controller.begin();
	for( auto iter = control_lattice.density_domains.begin(); iter < control_lattice.density_domains.end(); ++iter )
	{
		erase = true;
		for( int i = 0; i < control_lattice.number_of_cell_types; ++i )
		{
			pC = new Density_Controller;
			if( control_lattice.cell_types_defined )
			{
				pC->set_cell_type( control_lattice.cell_type_index[i], control_lattice.cell_type_name[i] );
				pC->control_density = control_lattice.cell_type_density[i];
				if( control_lattice.cell_type_index[i] == 1 )
				{
					pCF = new PI_Control_1;
//					pCF = new Antithetic_Control;
				}
				else if( control_lattice.cell_type_index[i] == 2 )
				{
					pCF = new PI_Control_1;
//					pCF = new Antithetic_Control;
				}
				else if( control_lattice.cell_type_index[i] == 3 )
				{
					pCF = new PI_Control_1;
//					pCF = new Antithetic_Control;
				}
				pCF->erase_light_pattern = erase;
				erase = erase ? !erase : erase; // only the first controller erases the pattern
			}
			else
			{
				pC->control_density = control_lattice.controller_param[0];
				pCF = new PI_Control_1;
//				pCF = new Antithetic_Control;
			}
			// For PI Control
			pCF->Kp = control_lattice.controller_param[1];
			pCF->Ki = control_lattice.controller_param[2];

			// For Antithetic Control
//			pCF->mu = control_lattice.controller_param[1];
//			pCF->theta = control_lattice.controller_param[2];
//			pCF->eta = control_lattice.controller_param[3];

			pC->set_up( {(*ll)}, pCF );
			pC->set_domain((*iter));
			pC->on = true;
			(*dc) = pC;
			++dc;
		}
		++ll;
	}
	return control_lattice;
}

void setup_controll_system( void )
{
	std::vector<Point> nodes(4);
	Point P0, P1;
	double la, lb;
	double control_density, Kp,Ki;
	int cell_number;

	// set domains for spatial control
	Circle* circle0 {new Circle};
	Circle* circle1 {new Circle};
	Circle* circle2 {new Circle};
	Circle* circle3 {new Circle};
	Circle* circle4 {new Circle};
	Circle* circle5 {new Circle};

	Domain_Intersection* rest {new Domain_Intersection};
	Domain_Intersection* ring0 {new Domain_Intersection};
	Domain_Intersection* ring1 {new Domain_Intersection};
	Domain_Intersection* ring2 {new Domain_Intersection};
	Domain_Intersection* ring3 {new Domain_Intersection};
	Domain_Intersection* ring4 {new Domain_Intersection};
	rest->set_name("Rest");
	ring1->set_name("Ring0");
	ring1->set_name("Ring1");
	ring2->set_name("Ring2");
	ring3->set_name("Ring3");
	ring4->set_name("Ring4");

	Rectangle* rec1 {new Rectangle};
	Rectangle* rec2 {new Rectangle};
	Rectangle* rec3 {new Rectangle};
	Rectangle* rec4 {new Rectangle};
	Rectangle* rec5 {new Rectangle};

	Default_Domain* default_domain {new Default_Domain}; // default domain is the entire simulation domain
	Domain_Union* cross {new Domain_Union};
	cross->set_name("Cross");
	Domain_Intersection* flag_without_cross {new Domain_Intersection};
	flag_without_cross->set_name("Flag_without_cross");

	// Define a circle
	// Note the radii are swt such that they fit on a 25x25 grid on the 1400x1400 domain
	Point P{0.0,0.0};
	double radius = 644;
	circle0->set_up( P, radius, "Disk0" );

	radius = 592;
	circle1->set_up( P, radius, "Disk1" );

	// Define another circle
	radius = 532;
	circle2->set_up( P, radius, "Disk2" );

	// Define another circle
	radius = 364;
	circle3->set_up( P, radius, "Disk3" );

	// Define another circle
	radius = 252;
	circle4->set_up( P, radius, "Disk4" );

	// Define another circle
	radius = 196;
	circle5->set_up( P, radius, "Central_Disk" );

	// This is the outside
	rest->add_domain(circle0, true);
	rest->set_area(default_domain->get_area() - circle0->get_area());

	// Define now the first ring
	ring0->add_domain(circle0);
	ring0->add_domain(circle1, true);
	ring0->set_area(circle0->get_area() - circle1->get_area());

	// Define now the first ring
	ring1->add_domain(circle1);
	ring1->add_domain(circle3, true);
	ring1->set_area(circle1->get_area() - circle3->get_area());

	// Define now the second ring
	ring2->add_domain(circle3);
	ring2->add_domain(circle5, true);
	ring2->set_area(circle3->get_area() - circle5->get_area());

	// Define now the third ring
	ring3->add_domain(circle3);
	ring3->add_domain(circle4, true);
	ring3->set_area(circle3->get_area() - circle4->get_area());

	// Define now the fourth ring
	ring4->add_domain(circle4);
	ring4->add_domain(circle5, true);
	ring4->set_area(circle4->get_area() - circle5->get_area());

	P0.set(-140,-350);
	P1.set(140,-350);
	rec1->set_up({P0,P1}, 700 );
	rec2->set_up({P0,P1}, 700 );
	rec2->rotate(90);

	P0.set(-630,-490);
	P1.set(630,-490);
	rec3->set_up({P0,P1}, 980, "Flag" );
	cross->add_domain(rec1);
	cross->add_domain(rec2);
	cross->set_area( 313600 );
	flag_without_cross->add_domain(rec3);
	flag_without_cross->add_domain(rec1, true);
	flag_without_cross->add_domain(rec2,true);
	flag_without_cross->set_area( 921200 );


	// Setting for local in the low half and global control in the upper half
	P0.set(-700, 0);
	P1.set(700, 0);
	rec4->set_up({P0,P1}, 700, "Upper half" );
	P0.set(-700, -700);
	P1.set(700, -700);
	rec5->set_up({P0,P1}, 700, "Lower half" );

	Kp = 0.9;
	Ki = 0.002;
//	cell_number = floor(1400*1400/(0.9*6*6*3.1415));
	cell_number = 4000;

	// This builds a control lattice for local control
	Control_Lattice control_lattice;

	// First define the domain
	control_lattice.domain = rec5;
	// how many lattice sides should be constructed (nxn)
	control_lattice.number_of_density_domains = 20;
	// how many LED slots should be on each domain
	control_lattice.number_of_sides = 5;
	// which types of LEDs should be there
	control_lattice.LED_types = {blue_led,red_led,green_led};
	control_density = cell_number/default_domain->get_area();

	control_lattice.controller_param = {control_density, Kp, Ki};

	// definition for controlling different cell types
	control_lattice.number_of_cell_types = 2;
	control_lattice.cell_type_index.resize(control_lattice.number_of_cell_types);
	control_lattice.cell_type_name.resize(control_lattice.number_of_cell_types);
	control_lattice.cell_type_density.resize(control_lattice.number_of_cell_types);
	control_lattice.cell_type_index[0] = 1;
	control_lattice.cell_type_name[0] = "HEK293_cell_A";
	control_lattice.cell_type_density[0] = control_density;
	control_lattice.cell_type_index[1] = 2;
	control_lattice.cell_type_name[1] = "HEK293_cell_B";
	control_lattice.cell_type_density[1] = 0.5*control_density;
	control_lattice.cell_types_defined = true;

	// this creates now the lattice. It takes the struct control_lattice and fills it further
	control_lattice = create_control_lattice( control_lattice );

	// plug it into the light control system
	for( auto iter = control_lattice.led_lattices.begin(); iter < control_lattice.led_lattices.end(); ++iter )
	{
			opto.add_light((*iter));
	}
	for( auto iter = control_lattice.controller.begin(); iter < control_lattice.controller.end(); ++iter )
	{
		opto.add_controller((*iter));
	}

	LED_Lattice* ll;
	// create LED lattices
	ll = new LED_Lattice;
	// couple it to a domain
	ll->set_domain(rec4);
	// set the number of lattice sides nxn
	ll->set_lattice(100);
	//	add LEDs to each lattice side
	ll->add_LED(blue_led);
	ll->add_LED(red_led);
	ll->add_LED(green_led);
	// Switches only the device on, not the individual LEDs
	ll->switch_on();
	// add it to the control system
	opto.add_light( ll );

	// Set up the controller
	Density_Controller* pC;
	PI_Control* pCF;
	// For the first cell type
	pC =  new Density_Controller;
	pCF = new PI_Control_1;
	pCF->Kp = 0.9;
	pCF->Ki = 0.002;
	pCF->erase_light_pattern = true;
	pC->set_up( {ll}, pCF );
	pC->set_domain(rec4);
	pC->set_cell_type( 1, "HEK293_cell_A" );
	pC->control_density = control_density;
	pC->on = true;
	opto.add_controller( pC );

	// For the second cell type
	pC =  new Density_Controller;
	pCF = new PI_Control_1;
	pCF->Kp = 0.9;
	pCF->Ki = 0.002;
	pCF->erase_light_pattern = false;
	pC->set_up( {ll}, pCF );
	pC->set_domain( rec4 );
	pC->set_cell_type( 2, "HEK293_cell_B" );
	pC->control_density = 0.5*control_density;
	pC->on = true;
	opto.add_controller( pC );
	// done

//	std::cout << "Setting up illumination domain 1: " << circle->get_type() << std::endl;
//	std::cout << "center = " << circle->get_center().tuple() << ", Radius = " << radius << std::endl;
//	std::cout << "Area = " << circle->get_area() << std::endl;
//	std::cout << std::endl;

	// set up the lights
	// single LEDs
//	LED* led_1 {new LED( "BLUE LED", {500,10} )};
//	LED led_1("BLUE LED", {500,10});
//	LED* led_2 {new LED( "RED LED", {660,10} )};
//	LED* led_3 {new LED( "RED LED", {660,10} )};

	// Link the domains to the light sources
//	led_1.set_domain(default_domain);
//	led_2->set_domain(rectangle1);
//	led_3->set_domain(rectangle2);

//	std::cout << led_lattice->get_type() << " defined on the " <<
//			led_lattice->get_domain_type() << ": " << led_lattice->get_domain_name() << std::endl;
//	std::cout << "with " << led_lattice->lattice.size() << " " <<
//			led_lattice->lattice[0].get_type() << "s" << std::endl;
//	std::cout << "Each " << led_lattice->lattice[0].get_type() << " is a " <<
//			led_lattice->lattice[0].domain.get_type() << std::endl;
//	std::cout << "and has " << led_lattice->lattice[0].leds.size() << " " <<
//			led_lattice->lattice[0].leds[0].get_type() << std::endl;
//	std::cout << std::endl;

	// produce a checker board pattern (needs an odd number of lattice sides, e.g., 3x3 or 7x7)
//	bool c = false;
//	for( auto iter = led_lattice->lattice.begin(); iter < led_lattice->lattice.end(); ++iter )
//	{
//		if(c)
//		{
//			(*iter).switch_on();
//			c = !c;
//		}
//		else
//		{
//			c = !c;
//		}
//	}

	// Set internal time control
//	led_1->time = []( double t ){ return t > 720 ? true : false; };

	// switch the lights on (off by default)
//	led_1->switch_on();
//	led_1->invert = true;
//	led_2->switch_off();
//	led_3->switch_off();

	// Combine it to the illumination system
//	opto.add_light(led_1);
//	opto.add_light(led_2);
//	opto.add_light(led_3);

	/* This builds controllers for the density
	 * on specified domains
	 */
/*
	// Set up the LED lattice
	// create LED lattices
	LED_Lattice* led_lattice {new LED_Lattice};
	// couple it to a domain
	led_lattice->set_domain(default_domain);
	// set the number of lattice sides nxn
	led_lattice->set_lattice(130);
	//	add LEDs to each lattice side
	led_lattice->add_LED(blue_led);
	led_lattice->add_LED(red_led);
	led_lattice->add_LED(green_led);
	// Switches only the device on, not the individual LEDs
	led_lattice->switch_on();
//	led_lattice->switch_all_sides_on();
	// add it to the control system
	opto.add_light( led_lattice );

	// we use this more often it should be better formalized
	// Now, set up the controller
	std::vector<Domain*> control_domains;
	control_domains.push_back(rest);
	control_domains.push_back(ring0);
	control_domains.push_back(ring1);
	control_domains.push_back(ring2);
//	control_domains.push_back(ring3);
//	control_domains.push_back(ring4);
	control_domains.push_back(circle5);

	std::vector<std::vector<double>> densities;
	std::vector<double> density_1(control_domains.size());
	std::vector<double> density_2(control_domains.size());
	std::vector<double> density_3(control_domains.size());
	density_1 = {0, 0, 0.7, 0, 0};
	densities.push_back(density_1);
	density_2 = {0, 0, 0, 0.7, 0};
	densities.push_back(density_2);
	density_3 = {0, 1, 0.3, 0.3, 1};
	densities.push_back(density_3);

	Density_Controller* pC;
	PI_Control_1* pCF;
	control_density = floor(1400*1400/(0.9*6*6*3.1415))/default_domain->get_area();
	int i = 0;
	bool erase = true;
	for( auto iter = control_domains.begin(); iter < control_domains.end(); ++iter )
	{
		int j=1;
		for( auto ds = densities.begin(); ds < densities.end(); ++ds )
		{
			pC =  new Density_Controller;
			pCF = new PI_Control_1;
			pCF->Kp = 0.9;
			pCF->Ki = 0.002;
			pCF->erase_light_pattern = erase;
			pC->set_domain( (*iter) );
			pC->set_up( {opto.lights[0]}, pCF );
			pC->set_cell_type(j);
//			std::cout << (*ds)[i] << std::endl;
			pC->control_density = (*ds)[i]*control_density;
			pC->on = true;
			erase = erase ? !erase : erase; // only the first controller erases the pattern
			++j;
			// add it to the control system
			opto.add_controller( pC );
		}
		++i;
	}
*/




/*
	// This builds a control lattice on the whole simulation domain
	Control_Lattice control_lattice;
	double control_density, Kp,Ki;
	int cell_number;
	double mu =  parameters.doubles("HEK293_cell_antithetic_mu_Z1");
	double theta = parameters.doubles("HEK293_cell_antithetic_theta_Z2");
	double eta = parameters.doubles("HEK293_cell_antithetic_eta_annihilation");
	Kp = 0.9;
	Ki = 0.002;
	cell_number = floor(1400*1400/(0.9*6*6*3.1415));
	// First define the domain
	control_lattice.domain = default_domain;
	// how many lattice sides should be constructed (nxn)
	control_lattice.number_of_density_domains = 35;
	// how many LED slots should be on each domain
	control_lattice.number_of_sides = 4;
	// which types of LEDs should be there
	control_lattice.LED_types = {blue_led,red_led,green_led};

	control_density = cell_number/default_domain->get_area();

	// for the antithetic control it is important to adjust the parameters mu and theta
	// according to the desired density. I keep theta fixed and adjust mu.
	mu = control_density * control_lattice.domain->get_area() * theta /
			(control_lattice.number_of_density_domains*control_lattice.number_of_density_domains);


//	control_lattice.controller_param = {control_density, Kp, Ki};
	control_lattice.controller_param = {control_density, mu, theta, eta};

	// definition for controlling different cell types
	control_lattice.number_of_cell_types = 3;
	control_lattice.cell_type_index.resize(control_lattice.number_of_cell_types);
	control_lattice.cell_type_name.resize(control_lattice.number_of_cell_types);
	control_lattice.cell_type_index[0] = 1;
	control_lattice.cell_type_name[0] = "HEK293_cell_A";
	control_lattice.cell_type_index[1] = 2;
	control_lattice.cell_type_name[1] = "HEK293_cell_B";
	control_lattice.cell_type_index[2] = 3;
	control_lattice.cell_type_name[2] = "HEK293_cell_C";
	control_lattice.cell_types_defined = true;

	// this creates now the lattice. It takes the struct control_lattice and fills it further
	control_lattice = create_control_lattice( control_lattice );

	for( auto iter = control_lattice.led_lattices.begin(); iter < control_lattice.led_lattices.end(); ++iter )
	{
		opto.add_light((*iter));
	}

	// set up different density regions
	Density_Controller* dc;
	for( auto iter = control_lattice.controller.begin(); iter < control_lattice.controller.end(); ++iter )
	{
		dc = static_cast<Density_Controller*>((*iter));
		// different target densities in different domains
		dc->control_density = 0.0*control_density;

		Swiss Cross
		if( (*rec3)((*iter)->get_center()) )
		{
			if( dc->get_cell_type_index()  == 2 ){ dc->control_density = control_density; }
			else{ dc->control_density = 0.0*control_density; }
		}
		if( (*rec1)((*iter)->get_center()) || (*rec2)((*iter)->get_center()) )
		{
			if( dc->get_cell_type_index() == 1 ){ dc->control_density = 1.5*control_density; }
			else{ dc->control_density = 0.0*control_density; }
		}

		if( (*circle0)((*iter)->get_center())  && !(*circle1)((*iter)->get_center()) )
		{
			if( dc->get_cell_type_index()  == 3 ){ dc->control_density = control_density; }
			else{ dc->control_density = 0.0*control_density; }
		}
		else if( (*circle1)((*iter)->get_center())  && !(*circle2)((*iter)->get_center()) )
		{
			if( dc->get_cell_type_index()  == 1 ){ dc->control_density = 0.7*control_density; }
			else if( dc->get_cell_type_index()  == 3 ){ dc->control_density = 0.3*control_density; }
			else{ dc->control_density = 0.0*control_density; }
		}
		else if( (*circle2)((*iter)->get_center()) && !(*circle3)((*iter)->get_center()) )
		{
			if( dc->get_cell_type_index()  == 1 ){ dc->control_density = 0.7*control_density; }
			else if( dc->get_cell_type_index()  == 3 ){ dc->control_density = 0.3*control_density; }
			else{ dc->control_density = 0.0*control_density; }
		}
		else if( (*circle3)((*iter)->get_center()) && !(*circle4)((*iter)->get_center()) )
		{
			if( dc->get_cell_type_index()  == 2 ){ dc->control_density = 0.7*control_density; }
			else if( dc->get_cell_type_index()  == 3 ){ dc->control_density = 0.3*control_density; }
			else{ dc->control_density = 0.0*control_density; }
		}
		else if( (*circle4)((*iter)->get_center()) && !(*circle5)((*iter)->get_center()) )
		{
			if( dc->get_cell_type_index()  == 2 ){ dc->control_density = 0.7*control_density; }
			else if( dc->get_cell_type_index()  == 3 ){ dc->control_density = 0.3*control_density; }
			else{ dc->control_density = 0.0*control_density; }
		}
		else if( (*circle5)((*iter)->get_center()) )
		{
			if( dc->get_cell_type_index()  == 3 ){ dc->control_density = control_density; }
			else{ dc->control_density = 0.0*control_density; }
		}

		//plug in the controllers
		opto.add_controller((*iter));
	}

*/

//	std::cout << "bi= " << index_blue << ", ri= " << index_red << ", gi= " << index_green << std::endl;
//	static_cast<LED_Lattice*>(opto.lights[0])->switch_all_sides_on(index_red);
//	static_cast<LED_Lattice*>(opto.lights[0])->switch_all_sides_on(index_blue);


//	LED_Lattice* pl = static_cast<LED_Lattice*>(opto.lights[0]);
//	std::cout << pl->get_led_type(1) << std::endl;
//	std::cout << pl->get_led_index("BLUE LED") << std::endl;
//	std::cout << pl->lattice[0].get_led_type(1) << std::endl;
//	std::cout << pl->lattice[0].get_led_index("BLUE LED") << std::endl;


//	Rectangle* rectangle1 {new Rectangle};
//	Rectangle* rectangle2 {new Rectangle};
//	Rectangle* rectangle3 {new Rectangle};
//	Rectangle* rectangle4 {new Rectangle};
//	Rectangle* rectangle5 {new Rectangle};
//
//	/*
//	// create density control lattices
//	int ln=18;
//	cell_number = 0;
//	Rectangle* lattice{new Rectangle[ln]};
//	la = 140;
//	P0.set(-la/2,la/2); P1.set(-la/2,-la/2); nodes = {P0,P1};
//	for( int i=0; i<ln; ++i )
//	{
//		lattice[i].set_up( nodes, la );
//	}
//	// define the letter C
//	lattice[0].move({-3*la,0});
//	lattice[1].move({-3*la,la});
//	lattice[2].move({-3*la,2*la});
//	lattice[3].move({-2*la,2*la});
//	lattice[4].move({-1*la,2*la});
//	lattice[5].move({-3*la,-la});
//	lattice[6].move({-3*la,-2*la});
//	lattice[7].move({-2*la,-2*la});
//	lattice[8].move({-1*la,-2*la});
//	// define the F
//	lattice[9].move({1*la,-2*la});
//	lattice[10].move({1*la,-1*la});
//	lattice[11].move({1*la,0*la});
//	lattice[12].move({1*la,1*la});
//	lattice[13].move({1*la,2*la});
//	lattice[14].move({2*la,0*la});
//	lattice[15].move({3*la,0*la});
//	lattice[16].move({2*la,2*la});
//	lattice[17].move({3*la,2*la});
//
////		rectangle1->set_up( nodes, 600 );
////		rectangle2->set_up( nodes, 600 );
////		rectangle2->move({600,0}); // move right
////		rectangle3->set_up( nodes, 600 );
////		rectangle3->move({600,-600}); // move right and down
////		rectangle4->set_up( nodes, 600 );
////		rectangle4->move({0,-600}); // move down
//
//	for( int n=0; n<ln; ++n)
//	{
//		cell_number = 2000 + 0*n*500;
//		std::cout << n << " " << cell_number << " " <<
//				cell_number/default_domain->get_area() << " " << lattice[n].get_center().tuple() << std::endl;
//		control_lattice.domain = &(lattice[n]);
//		control_lattice.number_of_density_domains = 2;
//		control_lattice.number_of_sides = 8;
//		control_lattice.LED_types = {red_led};
//		control_density = cell_number/default_domain->get_area();
//		control_lattice.controller_param = {control_density, Kp, Ki};
//		control_lattice = create_control_lattice( control_lattice );
//		for( int i = 0; i < control_lattice.number_of_density_domains*control_lattice.number_of_density_domains; ++i )
//		{
//			opto.add_light(control_lattice.led_lattices[i]);
//			opto.add_controller(control_lattice.controller[i]);
//		}
//	}
//*/
//
//	P0.set(-700,700); P1.set(-700,0); nodes = {P0,P1};
//	rectangle1->set_up( nodes, 1400 );

//
//	P0.set(-700,0); P1.set(-700,-700); nodes = {P0,P1};
//	rectangle2->set_up( nodes, 1400 );
//	cell_number = 3000;
//	control_lattice.domain = rectangle2;
//	control_lattice.number_of_density_domains = 16;
//	control_lattice.number_of_sides = 8;
//	control_lattice.LED_types = {red_led};
//	control_density = cell_number/default_domain->get_area();
//	control_lattice.controller_param = {control_density, Kp, Ki};
//	control_lattice = create_control_lattice( control_lattice );
//	for( int i = 0; i < control_lattice.number_of_density_domains*control_lattice.number_of_density_domains; ++i )
//	{
//		opto.add_light(control_lattice.led_lattices[i]);
//		opto.add_controller(control_lattice.controller[i]);
//	}
//
//
//	LED* led {new LED( "BLUE LED", {500,10} )};
//	P0.set(-650,650); P1.set(-650,-650); nodes = {P0,P1};
//	rectangle5->set_up( nodes, 1300 );
//	led->set_domain(rectangle5);
//	led->switch_off();
//	led->invert = true;
//	opto.add_light(led);

	opto.density.resize(2);
	opto.density[0].set_grid_domain( default_domain, 40 );
	opto.density[0].set_measurement_domain( rec4, false );
	opto.density[0].control_density = control_density;
	opto.density[1].set_grid_domain( default_domain, 40 );
	opto.density[1].set_measurement_domain( rec5, false );
	opto.density[1].control_density = control_density;

//	opto.density[2].set_grid_domain( default_domain, 40 );
//	opto.density[2].set_measurement_domain( ring2, false );
//	opto.density[2].control_density = control_density;
//	opto.density[3].set_grid_domain( default_domain, 40 );
//	opto.density[3].set_measurement_domain( ring3, false );
//	opto.density[3].control_density = control_density;
//	opto.density[4].set_grid_domain( default_domain, 40 );
//	opto.density[4].set_measurement_domain( ring4, false );
//	opto.density[4].control_density = control_density;
//	opto.density[3].set_grid_domain( default_domain, 40 );
//	opto.density[3].set_measurement_domain( circle5, false );
//	opto.density[3].control_density = control_density;

/* Swiss Cross
	opto.density.resize(3);
	opto.density[0].set_grid_domain( default_domain, 20 );
	opto.density[0].set_measurement_domain( rec3, true );
	opto.density[0].control_density = 0.0*control_density;

	opto.density[1].set_grid_domain( default_domain, 20 );
	opto.density[1].set_measurement_domain( flag_without_cross, false );
	opto.density[1].control_density = control_density;

	opto.density[2].set_grid_domain( default_domain, 20 );
	opto.density[2].set_measurement_domain( cross, false	 );
	opto.density[2].control_density = 0.4*control_density;
*/

	// this simplifies the access, but should be finally removed to avoid the global variables
	index_red = static_cast<LED_Lattice*>(opto.lights[0])->get_led_index("RED LED");
	index_blue = static_cast<LED_Lattice*>(opto.lights[0])->get_led_index("BLUE LED");
	index_green = static_cast<LED_Lattice*>(opto.lights[0])->get_led_index("GREEN LED");

	// Set the update intervall
	opto.update_dt = PhysiCell::phenotype_dt;
}
