/*
 * control.cpp
 *
 *  Created on: 17.01.2020
 *  Updated on: 21.08.2021
 *  Author: cfleck
 */

#include "control.h"
#include "custom.h"

bool Local_DensityFunctor::set_grid_domain( Rectangle* _domain, int n )
{
	if( (*_domain).get_type() != "Rectangle" )
	{
		std::cout << "Cannot set up density measurement grid. Domain for density measurement grid is of type "
				<< (*_domain).get_type() << ", but must be of type Rectangle." << std::endl;
		grid_set = false;
		return grid_set;
	}
	grid_domain = _domain;
	grid_domain_center = _domain->get_center();
	grid_domain_name = _domain->get_name();
	std::cout << std::endl;
	std::cout << "Setting up domain for density measurement grid " << (*_domain).get_name() <<  std::endl;
	std::cout << "Creating underlying grid on " << grid_domain->get_name() << " for density measurement." << std::endl;
	std::cout << "Area of the grid domain = " << _domain->get_area() << " um^2." << std::endl;
	std::cout << "Number of grid sites = " << n*n << "" << "." << std::endl;
	std::cout << "Area of grid site = " << grid_domain->get_area()/(n*n) << " um^2." << std::endl;

	// delete previous grid
	for( auto iter = density_domains.begin();  iter < density_domains.end(); ++iter)
	{
		delete (*iter);
	}

	// build grid
	density_domains.resize(n*n);
	number_of_density_domains = n*n;
	density_functors.resize(n*n);

	// the reason I use an array of pointers to domain instead of an array of rectangles
	// is that the set_up function of the density functor expects a pointer to Domain
	for( int i = 0; i < n*n; ++i )
	{
		density_domains[i] = new Rectangle;
	}
	Point P0, P1;
	Rectangle primitive_domain;
	Rectangle* wd = grid_domain;

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

	auto ld = density_domains.begin();
	auto lf = density_functors.begin();
	for( int i = 0; i < n; ++i )
	{
		for( int j = 0; j < n; ++j )
		{
			// move domain
			*(*ld) = primitive_domain;
			translation = j*d2 + i*d1;
//			std::cout << translation[0] << " " << translation[1] << std::endl;
			(*ld)->move( translation );
//			std::cout << (*ld)->get_center().tuple() << std::endl;
			(*ld)->set_name(std::to_string(i*n+j));

			// connect to density functor
			(*lf).set_domain( (*ld) );
//			std::cout << "translation vector = " << translation << std::endl;
//			std::cout << (*iter)->get_nodes()[0].tuple() << std::endl;
//			std::cout << (*iter)->get_nodes()[1].tuple() << std::endl;
			++ld; ++lf;
		}
	}
	grid_set = true;
	return grid_set;
}

bool Local_DensityFunctor::set_measurement_domain( Domain* _domain, bool invert )
{
	if( !grid_domain )
	{
		std::cout << "Grid domain not set. Cannot set measurement domain" << std::endl;
		measurement_set = false;
		return measurement_set;
	}
	if( grid_domain->get_area() < _domain->get_area() )
	{
		std::cout << "Grid domain smaller as measurement domain. Cannot set measurement domain" << std::endl;
		measurement_set = false;
		return measurement_set;
	}

	measurement_domain = _domain;
	if( invert )
	{
		measurement_domain_type = "Inverted " + (*_domain).get_type();
		measurement_domain_area = grid_domain->get_area() - _domain->get_area();
		measurement_domain_name = "Inverse of " + (*_domain).get_name();
	}
	else
	{
		measurement_domain_area = _domain->get_area();
		measurement_domain_name = (*_domain).get_name();
	}

	std::cout << std::endl;
	std::cout << "Creating measurement domain on " << measurement_domain_name << " based on " <<
				grid_domain->get_name() << " for density measurement." << std::endl;
	std::cout << "Area of the measurement domain = " << measurement_domain_area << " um^2." << std::endl;

	measure_density_functors.erase( measure_density_functors.begin(), measure_density_functors.end() );

	for( auto iter = density_functors.begin(); iter < density_functors.end(); ++iter )
	{
		if( (*_domain).get_name() == "DefaultDomain" )
		{
			measure_density_functors.push_back( &(*iter) );
		}
		else if( invert ^ _domain->point_inside( (*iter).get_center() ) )
		{
//			std::cout << "Point inside " << (*iter).get_center().tuple() << std::endl;
			measure_density_functors.push_back( &(*iter) );
		}
	}
	std::cout << "Number of lattice sites = " << measure_density_functors.size() << "." << std::endl;
	measured_densities.resize( measure_density_functors.size() );
	measurement_set = true;
	return measurement_set;
}

std::vector<Density> Local_DensityFunctor::operator()( void )
{
	int tc = 0;
	double m = 0;
	double m2 = 0;
	double tmp;

	auto md = measured_densities.begin();
	for( auto iter = measure_density_functors.begin(); iter < measure_density_functors.end(); ++iter )
	{
		(*md) = (**iter)();
		tmp = (*md).density;
		m += tmp;
		m2 += tmp*tmp;
		tc += (*md).number_of_cells;
		++md;
//		std::cout << "Density in domain " << (*md).domain_name << " centered at "
//				<< (*md).domain_center.tuple() << " " << (*md).density << std::endl;

	}
	total_cell_number = tc;
	mean = m/measured_densities.size();
	variance = m2/(measured_densities.size()-1) - m*m/(measured_densities.size()*(measured_densities.size()-1));

	return measured_densities;
}

std::vector<Density> Local_DensityFunctor::operator()( int cell_type_index )
{
	int tc = 0;
	double m = 0;
	double m2 = 0;
	double tmp;

	auto md = measured_densities.begin();
	for( auto iter = measure_density_functors.begin(); iter < measure_density_functors.end(); ++iter )
	{
		(*md) = (**iter)( cell_type_index );
		tmp = (*md).density;
		m += tmp;
		m2 += tmp*tmp;
		tc += (*md).number_of_cells;
		++md;
//		std::cout << "Density in domain " << (*md).domain_name << " centered at "
//				<< (*md).domain_center.tuple() << " " << (*md).density << std::endl;

	}
	total_cell_number = tc;
	mean = m/measured_densities.size();
	variance = m2/(measured_densities.size()-1) - m*m/(measured_densities.size()*(measured_densities.size()-1));

	return measured_densities;
}

void Light_Control::add_controller( Controller_Base* controller )
{
	controllers.push_back( controller );
	++number_of_controllers;
	local_densities.resize(number_of_controllers);
	return;

}

void Light_Control::remove_controller( int n )
{
	controllers.erase( controllers.begin() + n );
	--number_of_controllers;
	local_densities.resize(number_of_controllers);
	return;
}

bool Light_Control::operator()( Point x )
{
	is_there_light = false; auto iter = lights.begin();
	do {
		is_there_light = (*iter)->check(x);
		++iter;
	}
	while(!is_there_light && (iter < lights.end()));
	return is_there_light;
}

bool Light_Control::operator()( Point x, double t )
{
	is_there_light = false; auto iter = lights.begin();
	do {
		is_there_light = (*iter)->check(x,t);
		++iter;
	}
	while(!is_there_light && (iter < lights.end()));
	return is_there_light;
}

bool Light_Control::operator()( Point x, int m )
{
	is_there_light = false; auto iter = lights.begin();
	do {
		is_there_light = (*iter)->check(x,m);
		++iter;
	}
	while(!is_there_light && (iter < lights.end()));
	return is_there_light;
}

bool Light_Control::operator()( Point x, double t, int m )
{
	is_there_light = false; auto iter = lights.begin();
	do {
		is_there_light = (*iter)->check(x,t,m); ++iter;
	}
	while(!is_there_light && (iter < lights.end()));
	return is_there_light;
}

double Light_Control::get_density( Point x, double w, int m )
{
	double light = 0.0;
	for( auto iter = lights.begin(); iter < lights.end(); ++iter )
	{
		light += (*iter)->check(x) ? (*iter)->get_density(x,w,m) : 0.0;
	}
	return light;
}

double Light_Control::get_density( Point x, double t, double w, int m )
{
	double light = 0.0;
	for( auto iter = lights.begin(); iter < lights.end(); ++iter )
	{
		light += ((*iter)->check(x,t)) ? (*iter)->get_density(x,t,w,m) : 0.0;
	}
	return light;
}

void Light_Control::set_intensity( double intensity, Point x, int m )
{
	for( auto iter = lights.begin(); iter < lights.end(); ++iter )
	{
		if((*iter)->check(x))
		{
			(*iter)->set_intensity(intensity,x,m);
		}
	}
}

double Light_Control::get_intensity( Point x, int m )
{
	double light = 0.0;
	for( auto iter = lights.begin(); iter < lights.end(); ++iter )
	{
		light += ((*iter)->check(x)) ? (*iter)->get_intensity(x,m) : 0.0;
	}
	return light;
}

double Light_Control::get_intensity( Point x, double t, int m )
{
	double light = 0.0;
	for( auto iter = lights.begin(); iter < lights.end(); ++iter )
	{
		light += ((*iter)->check(x,t,m)) ? (*iter)->get_intensity(x,m) : 0.0;
	}
	return light;
}

void Light_Control::update( double t )
{
	double delta = t - last_update_time;
	int write = 5;

	// Time to update the light controller?
	if ( fabs(delta-update_dt) >= update_tolerance*update_dt ) { return; }

	// Do the updating
//	std::cout << "Updating the light equipment. Time = " << t << std::endl;
	last_update_time = t;

//	Control is only on between t1 and t2
	if( t >360 && number_of_controllers > 0 )
	{
		for( auto iter = controllers.begin(); iter < controllers.end(); ++iter)
		{
			(*iter)->run(t);
		}
	}

	// document the density
	counter++;
	if ( counter == write )
	{
		for( int i = 0; i < 2; ++i )
		{
			cell_number_file[i] << t;
			density_file[i] << t;
			variance_file[i] << t;
			for( auto iter = density.begin(); iter < density.end(); ++iter )
			{
				local_densities = (*iter)(i+1);
				std::cout << "Number of cells of type " << i+1 << " on domain \"" << (*iter).get_measurement_domain_name() << "\" = " << (*iter).get_total_cell_number() <<
//					" (should be: "<< floor((*iter).control_density * (*iter).get_measurement_domain_area()) << ")" <<
				std::endl;

				cell_number_file[i] << ", " << (*iter).get_total_cell_number();
				density_file[i] << ", " << (*iter).get_mean();
				variance_file[i] << ", " << (*iter).get_variance();
			}
			cell_number_file[i] << std::endl;
			density_file[i] << std::endl;
			variance_file[i] << std::endl;
			std::cout << std::endl;
		}
		std::cout << std::endl;
		counter = 0;
	}
/*
	// this is an example how parameters could be adjusted based on measured values/features
	bool announced = true;
	static double reduced_density_below_threshold = false;
	static double CV_threshold = 0.45;

	Local_DensityFunctor& pD = density[0];
	pD(1);
	Density_Controller* pC;
	if( (sqrt(pD.get_variance())/pD.get_mean() < CV_threshold ) && !reduced_density_below_threshold && (t > 360) )
	{
		CV_threshold *= 0.999;
		for( auto iter = controllers.begin(); iter < controllers.end(); ++iter)
		{
			pC = static_cast<Density_Controller*>(*iter);
			if( pD.point_inside_measurement_domain(pC->get_center()) && pC->get_cell_type_index() == 2 )
			{
				if( pC->control_density < 0.00025 ){ reduced_density_below_threshold = true; }
				if( announced )
				{
					std::cout << "Reducing control density for cell type 2 from " << pC->control_density << " to " <<
												0.99*pC->control_density << std::endl;
					announced = !announced;
				}
				pC->control_density *= 0.99;
			}
		}
	}
*/


	// switching something on and off
//	if( t > 1440 && t < 2880 && state_changed == false )
//	{
//		static_cast<LED_Lattice*>(lights[0])->switch_all_sides_on(index_red);
//		state_changed = true;
//	}
//	else if( t > 4320 && t < 5760 && state_changed == true )
//	{
//		static_cast<LED_Lattice*>(lights[0])->switch_all_sides_on(index_blue);
//		state_changed = false;
//	}
//	else if( t > 5760 && t < 7200 && state_changed == false )
//	{
//		static_cast<LED_Lattice*>(lights[0])->switch_all_sides_off(index_red);
//		state_changed = true;
//	}
//	else if( t > 7200 && state_changed == true )
//	{
//		static_cast<LED_Lattice*>(lights[0])->switch_all_sides_off(index_blue);
//		state_changed = false;
//	}

//	this changes the density outside of the circle after 2 days
//	if( t> 2880 && state_changed == false )
//	{
//		std::cout << "Changing desity" << std::endl;
//		Point P{0,0};
//		double radius = 300;
//		Circle* circle {new Circle};
//		circle->set_up( P, radius, "centre" );
//		Density_Controller* pC;
//		for( auto iter = controllers.begin(); iter < controllers.end(); ++iter)
//		{
//			pC = static_cast<Density_Controller*>(*iter);
//			if(!(*circle)(pC->get_center()))
//			{
//				pC->control_density =0.002;
//			}
//		}
//		state_changed = true;
//	}

//	if( t > 360 && !lights.back()->check() )
//	{
//		lights.back()->switch_on();
//		lights.back()->invert = true;
//	}
}

void Density_Controller::set_cell_type( int _cell_type_index, std::string _cell_type_name )
{
	cell_type_index = _cell_type_index;
	cell_type_name = _cell_type_name;
	cell_type_set = true;
	return;
}

void Density_Controller::set_up( std::vector<Light*> _lights, ControlFunctor_Base* _controll_functor )
{
	lights = _lights;
	control_functor = _controll_functor;
	set= true;
	return;
}

void Density_Controller::run( double t )
{
	if(!on){ return; }
	if( cell_type_set ){ current_density = density( cell_type_index ); }
	else { current_density = density(); }
	error = (current_density.density - control_density)*current_density.area;
//	error = control_density > 0 ? (current_density.density - control_density)/control_density :
//			current_density.density*current_density.area;

	if( !PI_control && (fabs(error) < 0.2*control_density*current_density.area || fabs(error) < 2) )
	{
		static_cast<PI_Control_1*>(control_functor)->switch_integration_on();
		PI_control = true;
	}

	// do the actual control step
	(*control_functor)( this, error, t );
	return;
}

Density DensityFunctor::operator()( void )
{
	if(!domain)
	{
		std::cout << "Domain for cell density measurements not defined" << std::endl;
		return current_density;
	}

	// set current time
	current_density.time = PhysiCell_globals.current_time;
	current_density.area = (*domain).get_area();
	current_density.domain_name = (*domain).get_name();
	current_density.domain_type = (*domain).get_type();
	current_density.domain_center = (*domain).get_center();
	current_density.number_of_cells = 0;
	current_density.cells.clear();

	if( current_density.domain_name == "DefaultDomain" )
	{
		for( auto iter = (*all_cells).begin(); iter < (*all_cells).end(); ++iter )
		{
			if( !((*iter)->phenotype.death.dead) )
			{
				++current_density.number_of_cells;
				current_density.cells.push_back(*iter);
			}
		}
//		std::cout << "Number of cells = " << current_density.number_of_cells << std::endl;
		current_density.density = static_cast<double>(current_density.number_of_cells/current_density.area);
	}
	else
	{
		for( auto iter = (*all_cells).begin(); iter < (*all_cells).end(); ++iter )
		{
			P = {(*iter)->position[0], (*iter)->position[1]};
			if( !((*iter)->phenotype.death.dead) && (*domain).point_inside( P ) )
			{
				++current_density.number_of_cells;
				current_density.cells.push_back(*iter);
			}
		}
		current_density.density = static_cast<double>(current_density.number_of_cells/current_density.area);
	}
	history.push_back( current_density );
	return current_density;
}

Density DensityFunctor::operator()( int cell_type_index )
{
	if(!domain)
	{
		std::cout << "Domain for cell density measurements not defined" << std::endl;
		return current_density;
	}
//	std::cout << domain->get_name() << std::endl;
	// set current time
	current_density.time = PhysiCell_globals.current_time;

	current_density.area = (*domain).get_area();
	current_density.domain_name = (*domain).get_name();
	current_density.domain_type = (*domain).get_type();
	current_density.domain_center = (*domain).get_center();
	current_density.number_of_cells = 0;
	current_density.cells.clear();

	if( current_density.domain_name == "DefaultDomain" )
	{
		for( auto iter = (*all_cells).begin(); iter < (*all_cells).end(); ++iter )
		{
			if( !((*iter)->phenotype.death.dead) && ((*iter)->type == cell_type_index) )
			{
				++current_density.number_of_cells;
				current_density.cells.push_back(*iter);
			}
		}
//		std::cout << "Number of cells = " << current_density.number_of_cells << std::endl;
		current_density.density = static_cast<double>(current_density.number_of_cells/current_density.area);
	}
	else
	{
		for( auto iter = (*all_cells).begin(); iter < (*all_cells).end(); ++iter )
		{
			P = {(*iter)->position[0], (*iter)->position[1]};
			if( !((*iter)->phenotype.death.dead) && ((*iter)->type == cell_type_index) && (*domain).point_inside( P ) )
			{
				++current_density.number_of_cells;
				current_density.cells.push_back(*iter);
			}
		}
//		std::cout << "Number of cells = " << current_density.number_of_cells << std::endl;
		current_density.density = static_cast<double>(current_density.number_of_cells/current_density.area);
	}
	history.push_back( current_density );
	return current_density;
}


void P_Control::operator()( Controller_Base* _controller, double _err, double _t )
{
	Density_Controller* controller = static_cast<Density_Controller*>(_controller);
	prev_t = t; t = _t;
	prev_err = err; err = _err;
	cn = floor( Kp*err );

	LED_Lattice* ll = static_cast<LED_Lattice*>(controller->lights[0]);
	if ( erase_light_pattern ){ ll->switch_all_sides_off(); }

	if( cn <= 0 )
	{
		std::cout << "Controller function: Nothing to do" << std::endl;
		std::cout << std::endl;
		return;
	}

	std::cout << "Controller function: control number = " << cn << std::endl;
	std::cout << std::endl;

	for( int i=1; i<=cn; ++i )
	{
//		floor((ll->get_lattice_size())*UniformRandom());
		ll->switch_side_on(floor((ll->get_lattice_size())*UniformRandom()));
	}
}

void PI_Control_1::operator()( Controller_Base* _controller, double _err, double _t )
{
	Density_Controller* controller = static_cast<Density_Controller*>(_controller);
	LED_Lattice* ll = static_cast<LED_Lattice*>(controller->lights[0]);
	std::vector<Cell*> cells = controller->current_density.cells;
	int max = floor( controller->current_density.number_of_cells );

	// erase previous light pattern
	if ( erase_light_pattern ){

//		std::cout << "Controller for " << controller->current_density.domain_name << " and cell type " <<
//				controller->get_cell_type_index() << " erasing light pattern." << std::endl;
		ll->switch_all_sides_off();
	}
//	else
//	{
//		std::cout << "Controller for " << controller->current_density.domain_name << " and cell type " <<
//						controller->get_cell_type_index() << " NOT erasing light pattern." << std::endl;
//	}

	prev_t = t; t = _t;
	prev_err = err; err = _err;

	if( integration ) { int_cn += Ki*(t-prev_t)*_err; }
	cn_i =  Kp*err + int_cn;

//	if(integration){ std::cout << "PI_CONTROL_1" << std::endl; }
//			else{ std::cout << "P_CONTROL_1" << std::endl; }
//	if( controller->current_density.domain_name == "Ring1" && controller->get_cell_type_index() == 3 )
//	{
//		std:: cout << Kp*err << " " << int_cn << " " << cn_i <<
//				" " << controller->control_density*controller->current_density.area << " " << max << std::endl;
//	}
//	if( controller->control_density == 0)
//	{
//		std:: cout << Kp*err << " " << int_cn << " " << cn_i <<
//					" " << controller->control_density*controller->current_density.area << " " << max << std::endl;
//	}
	if( cn_i < 0 )
	{
//		for( auto iter = cells.begin(); iter < cells.end(); ++iter )
//		{
//			ll->set_intensity(1.0,{(*iter)->position[0], (*iter)->position[1]}, index_green);
//		}
	}
	else if( (cn_i >= 0) && (cn_i < 2) )
	{
		// no proliferation
		for( auto iter = cells.begin(); iter < cells.end(); ++iter )
		{
			ll->set_intensity(0.5*cn_i,{(*iter)->position[0], (*iter)->position[1]}, index_red);
		}
	//		ll->set_all_intensity(0.5*cn_i,index_red); // no proliferation
		}
	// density too high, remove cells
	else if( cn_i >= 2 )
	{
//		std::cout << "Controller function: control number = " << cn_i << ", removing cells and switching off proliferation" << std::endl;
//		std::cout << std::endl;
//
		auto pcell = cells.begin();
		cn = (cn_i <= max) ? floor(cn_i) : max;
		for(int i = 0; i < cn; ++i)
		{
			pcell = cells.begin() + floor(cells.size()*UniformRandom());
//			ll->switch_side_on({(*pcell)->position[0], (*pcell)->position[1]}, index_blue);
			ll->set_intensity(1.0, {(*pcell)->position[0], (*pcell)->position[1]}, index_blue);
			cells.erase(pcell);
		}
		// no proliferation
		for( auto iter = cells.begin(); iter < cells.end(); ++iter )
		{
			ll->set_intensity(cn_i,{(*iter)->position[0], (*iter)->position[1]}, index_red);
			cells.erase(iter);
		}
//		for( auto iter = cells.begin(); iter < cells.end(); ++iter )
//		{
//			ll->set_intensity(1.0,{(*iter)->position[0], (*iter)->position[1]}, index_green);
//		}
	}
/*
	// density too low
	if( cn_i < 0 )
	{
//		std::cout << "Controller function: control number = " << cn_i << ", doing nothing" << std::endl;
//		std::cout << std::endl;
		return;
	}
	// density just right, inhibit proliferation
	else if( (cn_i >= 0) && (cn_i < 2) )
	{
//		std::cout << "Controller function: control number = " << cn_i << ", switching off proliferation" << std::endl;
//		std::cout << std::endl;
//		ll->switch_all_sides_on(index_red); // no proliferation
		// no proliferation
		for( auto iter = cells.begin(); iter < cells.end(); ++iter )
		{
			ll->set_intensity(0.5*cn_i,{(*iter)->position[0], (*iter)->position[1]}, index_red);
		}
//		ll->set_all_intensity(0.5*cn_i,index_red); // no proliferation
	}
	// density too high, remove cells
	else if( cn_i >= 2 )
	{
//		std::cout << "Controller function: control number = " << cn_i << ", removing cells and switching off proliferation" << std::endl;
//		std::cout << std::endl;

		auto pcell = cells.begin();
//		cn = controller->control_density > 0 ? cn_i*controller->control_density*controller->current_density.area
//				: cn_i;
		cn = (cn_i <= max) ? floor(cn_i) : max;
		for(int i = 0; i < cn; ++i)
		{
			pcell = cells.begin() + floor(cells.size()*UniformRandom());
//			ll->switch_side_on({(*pcell)->position[0], (*pcell)->position[1]}, index_blue);
			ll->set_intensity(1.0, {(*pcell)->position[0], (*pcell)->position[1]}, index_blue);
			cells.erase(pcell);
		}
		// no proliferation
		for( auto iter = cells.begin(); iter < cells.end(); ++iter )
		{
			ll->set_intensity(0.5*cn_i,{(*iter)->position[0], (*iter)->position[1]}, index_red);
		}
//		ll->switch_all_sides_on(index_red); // no proliferation
//		ll->set_all_intensity(1.0,index_red); // no proliferation
	}
*/
	return;
}

void PI_Control_2::operator()( Controller_Base* _controller, double _err, double _t )
{
	Density_Controller* controller = static_cast<Density_Controller*>(_controller);
	LED_Lattice* ll = static_cast<LED_Lattice*>(controller->lights[0]);
	std::vector<Cell*> cells = controller->current_density.cells;
	int max = floor( controller->current_density.number_of_cells );

	// erase previous light pattern
	if ( erase_light_pattern ){ ll->switch_all_sides_off(); }

	prev_t = t; t = _t;
	prev_err = err; err = _err;

	if( integration ) { int_cn += Ki*(t-prev_t)*_err; }
	cn_i =  Kp*err + int_cn;

//	if(integration){ std::cout << "PI_CONTROL_2" << std::endl; }
//			else{ std::cout << "P_CONTROL_2" << std::endl; }
//	std:: cout << Kp*err << " " << int_cn << " " << cn_i <<
//			" " << controller->control_density*controller->current_density.area << " " << max << std::endl;

	// density too low
	if( cn_i < 0 )
	{
//		std::cout << "Controller function: control number = " << cn_i << ", doing nothing" << std::endl;
//		std::cout << std::endl;
		return;
	}
	// density just right, inhibit proliferation
	else if( (cn_i >= 0) && (cn_i < 2) )
	{
//		std::cout << "Controller function: control number = " << cn_i << ", switching off proliferation" << std::endl;
//		std::cout << std::endl;
//		ll->switch_all_sides_on(index_red); // no proliferation
		// no proliferation
		for( auto iter = cells.begin(); iter < cells.end(); ++iter )
		{
			ll->set_intensity(cn_i,{(*iter)->position[0], (*iter)->position[1]}, index_red);
		}
//		ll->set_all_intensity(0.5*cn_i,index_red); // no proliferation
	}
	// density too high, remove cells
	else if( cn_i >= 2 )
	{
//		std::cout << "Controller function: control number = " << cn_i << ", removing cells and switching off proliferation" << std::endl;
//		std::cout << std::endl;

		auto pcell = cells.begin();
//		cn = controller->control_density > 0 ? cn_i*controller->control_density*controller->current_density.area
//				: cn_i;
		cn = (cn_i <= max) ? floor(cn_i) : max;
		for(int i = 0; i < cn; ++i)
		{
			pcell = cells.begin() + floor(cells.size()*UniformRandom());
//			ll->switch_side_on({(*pcell)->position[0], (*pcell)->position[1]}, index_blue);
			ll->set_intensity(1.0, {(*pcell)->position[0], (*pcell)->position[1]}, index_blue);
			cells.erase(pcell);
		}
		// no proliferation
		for( auto iter = cells.begin(); iter < cells.end(); ++iter )
		{
			ll->set_intensity(0.5*cn_i,{(*iter)->position[0], (*iter)->position[1]}, index_red);
		}
//		ll->switch_all_sides_on(index_red); // no proliferation
//		ll->set_all_intensity(1.0,index_red); // no proliferation
	}
	return;
}

void Antithetic_Control::operator()( Controller_Base* _controller, double _err, double _t )
{
	Density_Controller* controller = static_cast<Density_Controller*>(_controller);
	LED_Lattice* ll = static_cast<LED_Lattice*>(controller->lights[0]);
	std::vector<Cell*> cells = controller->current_density.cells;
	double max = controller->current_density.number_of_cells;
	double control_number = controller->control_density * controller->current_density.area;
	double intensity;
//	max = max - control_number > 0 ? max - control_number : 0;

	// erase previous light pattern
	if ( erase_light_pattern ){ ll->switch_all_sides_off(); }

	// update Z1
	if( UniformRandom() <= mu * 0.1 )
	{
		Z1++;
	}
//	if( UniformRandom() <= 0.001 * Z1 * 0.1 )
//	{
//		Z1--;
//	}
	// update Z2
	if( UniformRandom() <= theta * controller->current_density.number_of_cells * 0.1 )
	{
		Z2++;
	}
//	if( UniformRandom() <= 0.001 * Z2 * 0.1 )
//	{
//			Z2--;
//	}
	// annihilation step
	if( UniformRandom() <= eta * Z1 * Z2 * 0.1 )
	{
		Z1--;
		Z2--;
	}

//	if( _err < - 0.1 )
//	{
////		std::cout << "Controller function: control number = " << cn_i << ", doing nothing" << std::endl;
////		std::cout << std::endl;
////		return;
//	}
//	// density just right, inhibit proliferation
//	else if( (_err >= - 0.1) && (_err <= 0.1) )
//	{
////		std::cout << "Controller function: control number = " << cn_i << ", switching off proliferation" << std::endl;
////		std::cout << std::endl;
////		ll->switch_all_sides_on(index_red); // no proliferation
//		ll->set_all_intensity(0.5,index_red); // no proliferation
//	}
//	// density too high, remove cells
//	else if( _err > 0.1 )
//	{
//		std::cout << "Controller function: control number = " << cn_i << ", removing cells and switching off proliferation" << std::endl;
//		std::cout << std::endl;
		auto pcell = cells.begin();
//		std::cout << control_number << std::endl;
//		cn = Z2 < max ? floor( Z2 ) : max;

//		intensity = pow(static_cast<double>(Z2),3)/(pow(control_number,3) + pow(static_cast<double>(Z2),3));
		intensity = control_number > 0 ? static_cast<double>(Z2)/(control_number + static_cast<double>(Z2)) : 1.0 ;

//		std::cout << "domain: " << controller->current_density.domain_name <<
//				" centered at " << controller->current_density.domain_center.tuple() << std::endl;
//		ll->set_all_intensity(cn, index_blue);

		cn = floor( intensity * max );
//		std::cout << cn << " " << max << " " << Z1 << " " << Z2 << " " << intensity << std::endl;
//		std::cout << cn << std::endl;
		for(int i = 0; i < cn; ++i)
		{
//			if( control_number == 0)
//			{
//				std::cout << " Erasing cells" << std::endl;
//			}
			pcell = cells.begin() + floor(cells.size()*UniformRandom());
//			ll->switch_side_on({(*pcell)->position[0], (*pcell)->position[1]}, index_blue);
			ll->set_intensity(intensity, {(*pcell)->position[0], (*pcell)->position[1]}, index_blue);
			cells.erase(pcell);
		}
		// no proliferation
		for( auto iter = cells.begin(); iter < cells.end(); ++iter )
		{
			ll->set_intensity(0.5*intensity,{(*iter)->position[0], (*iter)->position[1]}, index_red);
		}
//		ll->switch_all_sides_on(index_red); // no proliferation
//		ll->set_all_intensity(1.0,index_red); // no proliferation
//	}
	return;
}

//void Antithetic_Control::operator()( Controller_Base* _controller, double _err, double _t )
//{
//	Antithetic_Controller* controller = static_cast<Antithetic_Controller*>(_controller);
//	LED_Lattice* ll = static_cast<LED_Lattice*>(controller->lights[0]);
//	std::vector<Cell*> cells = controller->current_density.cells;
//	ll->switch_all_sides_off();
//	int iZ1 = cells[0]->custom_data.find_variable_index( "Z1" );
//	int iZ2 = cells[0]->custom_data.find_variable_index( "Z2" );
//	double intensity;
//	int counter = 0;
//
//	for( auto iter = cells.begin(); iter < cells.end(); ++iter)
//	{
//		if( !(*iter)->phenotype.death.dead )
//		{
//			intensity += (*iter)->custom_data[iZ1];
//			++counter;
//		}
//		intensity /= counter;
//		intensity = 1.0/(1+intensity);
////		std::cout << intensity << std::endl;
//		ll->set_all_intensity(intensity, index_blue);
//	}
//	return;
//}

//void Antithetic_Control::operator()( Controller_Base* _controller, double _err, double _t )
//{
//	Antithetic_Controller* controller = static_cast<Antithetic_Controller*>(_controller);
//	LED_Lattice* ll = static_cast<LED_Lattice*>(controller->lights[0]);
//	std::vector<Cell*> cells = controller->current_density.cells;
//	int max = floor( controller->current_density.number_of_cells );
//	ll->switch_all_sides_off();
//	int iZ1 = cells[0]->custom_data.find_variable_index( "Z1" );
//	int iZ2 = cells[0]->custom_data.find_variable_index( "Z2" );
//	double intensity;
//
//	for( auto iter = cells.begin(); iter < cells.end(); ++iter)
//	{
//		intensity = 1.0/(1+(*iter)->custom_data[iZ1]);
////		std::cout << (*iter)->custom_data[iZ1] << std::endl;
//		ll->set_intensity(intensity, {(*iter)->position[0], (*iter)->position[1]}, index_blue);
////		if( (*iter)->custom_data[iZ1] < 10 )
////		{
////			std::cout << "Z1 = " << (*iter)->custom_data[iZ1] << " Z2 = " << (*iter)->custom_data[iZ2] << std::endl;
////			std::cout << std::endl;
////		}
//	}
//	return;
//}

/*
 * old function
void PI_Control::operator()( double _err, double _t, std::vector<Light*> lights, std::vector<Cell*> _cells )
{
	LED_Lattice* ll = static_cast<LED_Lattice*>(lights[0]);
	int max = floor( _cells.size() );
	ll->switch_all_sides_off();

	prev_t = t; t = _t;
	prev_err = err; err = _err;

	int_cn += Ki*(t-prev_t)*_err;

	cn_i =  Kp*err + int_cn;
	cn = (cn_i <= max) ? cn_i : max;

//	std:: cout << Kp*err << " " << int_cn << " " << max << " " << cn << std::endl;

	if( cn < 1 )
	{
		std::cout << "Controller function: Nothing to do. cn = " << cn << std::endl;
		std::cout << std::endl;
		return;
	}


//	std::cout << "Controller function: control number = " << cn << std::endl;
//	std::cout << std::endl;

	std::vector<Cell*> cells = _cells;
	auto pcell = cells.begin();
	std::vector<Cell*> sub_set(floor(cn));
	for(int i = 0; i < cn; ++i)
	{
		pcell = cells.begin() + floor(cells.size()*UniformRandom());
		ll->switch_side_on({(*pcell)->position[0], (*pcell)->position[1]}, index_blue);
//		sub_set[i] = *pcell;
		cells.erase(pcell);
	}
}

*/
