/*
 * light.cpp
 *
 *  Created on: 17.01.2020
 *  Updated on: 21.08.2021
 *  Author: cfleck
 */

#include "light.h"

LED::LED( const LED& led ) : gauss(led.get_spectral_param())
{
//	std::cout << "This is the copy constructor" << std::endl;
	type = led.type;
	spectral_type = led.spectral_type;
	on = led.on;
	intensity = led.intensity;
	invert = led.invert;
	domain = led.domain;
	spectrum = &gauss;
}

LASER::LASER( const LASER& laser ) : delta(laser.get_spectral_param())
{
//	std::cout << "This is the copy constructor" << std::endl;
	type = laser.type;
	spectral_type = laser.spectral_type;
	on = laser.on;
	intensity = laser.intensity;
	invert = laser.invert;
	domain = laser.domain;
	spectrum = &delta;
}

BAND::BAND( const BAND& band ) : uniform(band.get_spectral_param())
{
//	std::cout << "This is the copy constructor" << std::endl;
	type = band.type;
	spectral_type = band.spectral_type;
	on = band.on;
	intensity = band.intensity;
	invert = band.invert;
	domain = band.domain;
	spectrum = &uniform;
}

inline void Light::set_intensity( double _intensity )
{
		_intensity < max_intensity ? intensity = _intensity : intensity = max_intensity;
		// in case intensity of larger then 0, set the device on, otherwise leave the state untouched
		on = _intensity > 0.0 ? true : on;
}

inline void Light::set_all_intensity( double _intensity, int m )
{
	set_intensity( _intensity );
}

inline void Light::set_intensity( double _intensity, Point x, int m )
{
	if( (*domain)(x) )
	{
		set_intensity( _intensity );
	}
}

inline void Light::set_intensity( double _intensity, int n, int m )
{
	set_intensity( _intensity );
}

LED_Lattice::LED_Lattice() : lattice(0), led_types(0), led_index()
{
	type = "LED_Lattice";
}

bool LED_Lattice::set_lattice( int n )
{
	if( !domain )
	{
		std::cout << "Cannot set up lattice. No domain for lattice defined." << std::endl;
		set = false;
		return set;
	}
	if( (*domain).get_type() != "Rectangle" )
	{
		std::cout << "Cannot set up lattice. Domain for LED Lattice is of type " << (*domain).get_type() <<
				" but must be of type Rectangle." << std::endl;
		set = false;
		return set;
	}
	lattice.resize(n*n);
	Point P0, P1;
	Rectangle primitive_domain;
	Rectangle* wd = static_cast<Rectangle*>(domain);
	std::vector<Point> nodes(4);
	std::vector<double> translation(2);
	std::vector<double> d1(2);
	std::vector<double> d2(2);
	d1 = {wd->AB[0]/n, wd->AB[1]/n};
	d2 = {wd->AD[0]/n, wd->AD[1]/n};
//	std::cout << d1 << std::endl;
//	std::cout << d2 << std::endl;

	// construct the primitive lattice side
	P0 = wd->get_nodes()[0];
	P1 = P0; P1.move(d1);
	nodes = {P0,P1}; primitive_domain.set_up(nodes, sqrt(wd->AD_AD)/n);

//	P2 = P1; P2.move(d2);
//	P3 = P2; P3.move((-1)*d1);
//	std::cout << P0.tuple() << std::endl;
//	std::cout << P1.tuple() << std::endl;
//	std::cout << P2.tuple() << std::endl;
//	std::cout << P3.tuple() << std::endl;

	auto iter = lattice.begin();
	for( int i=0; i<n; ++i )
	{
		for( int j=0; j<n; ++j )
		{
			(*iter).domain = primitive_domain;
			translation = j*d2 + i*d1;
			(*iter).domain.move( translation );
//			std::cout << "translation vector = " << translation << std::endl;
//			std::cout << (*iter).domain.get_nodes()[0].tuple() << std::endl;
//			std::cout << (*iter).domain.get_nodes()[2].tuple() << std::endl;
			++iter;
		}
	}
	set = true;
	return set;
}

void LED_Lattice::add_LED( std::string type, std::vector<double> param )
{
	LED led(type,param);
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter)
	{
		// attach the led domain with the domain of the LEDs_Domain
		led.set_domain( &((*iter).domain) );
		(*iter).add_LED(led);
//		std::cout << (*iter).leds.back().get_spectral_type() << std::endl;
	}
	led_types.push_back(type);
	led_index[type] = number_leds;
	number_leds++;
}

void LED_Lattice::add_LED( Light_Type definition )
{
	add_LED( definition.type, definition.param );
}

inline bool LED_Lattice::check( Point x )
{
	is_there_light = false;
	if( !on ){ return is_there_light; }

	auto iter = lattice.begin();
	auto l = (*iter).leds.begin();
	do {
		do {
			is_there_light = l->check(x); ++l;
		}while(!is_there_light && l < (*iter).leds.end());
		++iter;
		l = (*iter).leds.begin();
	}
	while(!is_there_light && iter < lattice.end());
//	std::cout << x.tuple() << " " << is_there_light << std::endl;
	return is_there_light;
}

inline bool LED_Lattice::check( Point x, double t )
{
	is_there_light = false;
	if( !( (*time)(t) && on ) ){ return is_there_light; }

	auto iter = lattice.begin();
	auto l = (*iter).leds.begin();
	do {
		do {
			is_there_light = l->check(x,t); ++l;
		}while(!is_there_light && l < (*iter).leds.end());
		++iter;
		l = (*iter).leds.begin();
	}
	while(!is_there_light && iter<lattice.end());
	return is_there_light;
}

inline bool LED_Lattice::check( Point x, int m )
{
	is_there_light = false;
	if( !on ){ return is_there_light; }
	auto iter = lattice.begin();
	do {
		is_there_light = (*iter).leds[m].check(x);
		++iter;
	}
	while(!is_there_light && iter < lattice.end());
//	std::cout << x.tuple() << " " << is_there_light << std::endl;
	return is_there_light;
}

inline bool LED_Lattice::check( Point x, double t, int m )
{
	is_there_light = false;
	if( !( (*time)(t) && on ) ){ return is_there_light; }

	auto iter = lattice.begin();
	do {
		is_there_light = (*iter).leds[m].check(x,t);
		++iter;
	}
	while(!is_there_light && iter<lattice.end());
	return is_there_light;
}

LEDs_Domain* LED_Lattice::search_domain( Point x )
{
	LEDs_Domain* domain = NULL;
	auto iter = lattice.begin();
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter )
	{
		if((*iter).domain.point_inside(x))
		{
			domain = &(*iter);
			break;
		}
	}
	return domain;
}

void LEDs_Domain::add_LED( LED& led )
{
	led_types.push_back(led.get_type());
	led_index[led.get_type()] = number_leds;
	leds.push_back(led);
	number_leds++;
}

void LEDs_Domain::switch_on( int m )
{
	if( (m == -1) || (m + 1 > leds.size()) )
	{
		for( auto iter = leds.begin(); iter < leds.end(); ++iter)
		{
			(*iter).switch_on();
		}
	}
	else
	{
		leds[m].switch_on();
	}
}

void LEDs_Domain::switch_off( int m )
{
	if( (m == -1) || (m + 1 > leds.size()) )
	{
		for( auto iter = leds.begin(); iter < leds.end(); ++iter)
		{
			(*iter).switch_off();
		}
	}
	else
	{
		leds[m].switch_off();
	}
}

void LEDs_Domain::toggle( int m )
{
	if( (m == -1) || (m + 1 > leds.size()) )
		{
		for( auto iter = leds.begin(); iter < leds.end(); ++iter)
			{
				(*iter).switch_off();
			}
		}
		else
		{
			leds[m].switch_off();
		}
}

void LEDs_Domain::set_intensity( double _intensity, int m )
{
	if( (m == -1) || (m + 1 > leds.size()) )
	{
		for( auto iter = leds.begin(); iter < leds.end(); ++iter)
		{
			(*iter).set_intensity(_intensity);
		}
	}
	else
	{
		leds[m].set_intensity(_intensity);
	}
}

double LEDs_Domain::get_intensity( int m )
{
	return leds[m].get_intensity(m);
}

double LEDs_Domain::get_density( double w, int m )
{
	return leds[m].get_density(w);
}

double LEDs_Domain::get_density( Point x, double w, int m )
{
	return leds[m].get_density(x,w);
}

double LEDs_Domain::get_density( Point x, double t, double w, int m )
{
	return leds[m].get_density(x,t,w);
}

void LED_Lattice::switch_side_on( int n, int m )
{
	lattice[n].switch_on( m );
}

void LED_Lattice::switch_side_off( int n, int m )
{
	lattice[n].switch_off( m );
}

void LED_Lattice::switch_side_on( Point x, int m )
{
	auto iter = lattice.begin();
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter )
	{
		if((*iter).domain.point_inside(x))
		{
			(*iter).switch_on(m);
			break;
		}
	}
}

void LED_Lattice::switch_side_off( Point x, int m )
{
	auto iter = lattice.begin();
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter )
	{
		if((*iter).domain.point_inside(x))
		{
			(*iter).switch_off(m);
			break;
		}
	}
}

void LED_Lattice::toggle_side( Point x, int m )
{
	auto iter = lattice.begin();
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter )
	{
		if((*iter).domain.point_inside(x))
		{
			(*iter).toggle(m);
			break;
		}
	}
}

void LED_Lattice::switch_all_sides_on( int m )
{
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter)
	{
		(*iter).switch_on(m);
	}
}

void LED_Lattice::switch_all_sides_off( int m )
{
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter)
		{
			(*iter).switch_off(m);
		}
}

void LED_Lattice::toggle_side( int n, int m )
{
	lattice[n].toggle( m );
}

void LED_Lattice::set_all_intensity( double intensity, int m )
{
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter)
	{
		(*iter).set_intensity(intensity,m);
	}
}

void LED_Lattice::set_intensity( double intensity, int n, int m )
{
	lattice[n].set_intensity(intensity,m);
}

void LED_Lattice::set_intensity( double intensity, Point x, int m )
{
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter )
	{
		if((*iter).domain.point_inside(x))
		{
			(*iter).set_intensity(intensity,m);
			break;
		}
	}
}

double LED_Lattice::get_intensity( int n, int m )
{
	return lattice[n].get_intensity(m);
}

double LED_Lattice::get_intensity( Point x, int m )
{
	double light = 0.0;
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter )
	{
		if((*iter).domain.point_inside(x))
		{
			light = (*iter).get_intensity(m);
			break;
		}
	}
	return light;
}

double LED_Lattice::get_density( double w, int m )
{
//	the LEDs are the same for the whole grid,
//	it is therefore sufficient to check the first LEDs
	double light = 0.0;
	light = lattice[0].leds[m].check() ? lattice[0].leds[m].get_density(w) : 0.0;
	return light;
}

double LED_Lattice::get_density( Point x, double w, int m )
{
	LEDs_Domain* domain = NULL;
	double light = 0.0;
	auto iter = lattice.begin();
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter )
	{
		if((*iter).domain.point_inside(x))
		{
			light = (*iter).get_density(x,w,m);
			break;
		}
	}
	return light;
}

double LED_Lattice::get_density( Point x, double t, double w, int m )
{
	LEDs_Domain* domain = NULL;
	double light = 0.0;
	auto iter = lattice.begin();
	for( auto iter = lattice.begin(); iter < lattice.end(); ++iter )
	{
		if((*iter).domain.point_inside(x))
		{
			light = (*iter).get_density(x,t,w,m);
			break;
		}
	}
	return light;
}
