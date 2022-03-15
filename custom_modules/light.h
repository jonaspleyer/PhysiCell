/*
 * light.h
 *
 *  Created on: 17.01.2020
 *  Updated on: 21.08.2021
 *  Author: cfleck
 */

#ifndef CUSTOM_MODULES_LIGHT_H_
#define CUSTOM_MODULES_LIGHT_H_

#include "../core/PhysiCell.h"
#include "domain.h"
#include "spectrum.h"

#include <vector>
#include<map>

// For debugging
#ifndef _HI_
#define _HI_ std::cout << "HERE I AM" << std::endl;
#endif

class Light
{
protected:
	std::string type;
	std::string spectral_type;
	Spectral_Density *spectrum = NULL;
	Domain* domain = NULL;
	bool on = false; // a simple way to switch the light on and off
	double intensity;
	const double max_intensity = 1.0;
public:
	Light() : intensity{1.0} {};
	Light( std::string _type ) : intensity(0), type{_type} {};
	virtual ~Light() {};
	// TODO: add function f(x,y,z) to represent the spatial illumination instead of the boolean approach using Domain

	// TODO: add spatial density function and intensity parameter

	void set_sepctrum( Spectral_Density* _spectrum ) { spectrum = _spectrum; spectral_type = _spectrum->get_type(); }

	std::string get_type( void ) const { return type; }
	std::string get_spectral_type( void ) const { return spectral_type; }
	std::vector<double> get_spectral_param( void ) const { return (*spectrum).get_param(); }

	void set_domain( Domain* _domain ) { domain = _domain; }
	std::string get_domain_type( void ) { return (*domain).get_type(); }
	std::string get_domain_name( void ) { return (*domain).get_name(); }

	std::vector<Point> get_nodes( void ) { return (*domain).get_nodes(); }

	void move( std::vector<double> translation ){ (*domain).move(translation); }
	void rotate( double angle ){ (*domain).rotate(angle); }

	// TODO: think about a way to return the complete spectrum instead of the density at a given point
	// this could be done via a vector

	// overloaded functions for the density value for the wavelength w
	inline virtual double get_density( double w, int m = 0 ) { return (*spectrum)(w); }
	inline virtual double get_density( Point x, double w, int m = 0 ) { return ((*domain)(x) && on)? (*spectrum)(w): 0.0; } // questions the external switch (on)
	inline virtual double get_density( Point x, double t, double w, int m = 0 ) { return ((*domain)(x) && on && (*time)(t)) ? (*spectrum)(w) : 0.0; } // questions the internal time function

	// the arguments n and m are only needed for the LED lattice
	void set_intensity( double intensity );
	virtual void set_all_intensity( double _intensity, int m = -1 );
	virtual void set_intensity( double _intensity, Point x, int m = -1 );
	virtual void set_intensity( double _intensity, int n, int m = -1 );
	inline virtual double get_intensity( void ) { return on ? intensity : 0.0; }
	inline virtual double get_intensity( int n, int m = 0 ) { return on ? intensity : 0.0; }
	inline virtual double get_intensity( Point x, int m = 0 ) { return on && (*domain)(x) ? intensity : 0.0; }

	// for external switch on/off
	inline void switch_on( void ) { on = true; }
	inline void switch_off( void ) { on = false; }
	inline void toggle( void ) { on = !on; }

	// inverting the domain result
	bool invert = false;

	// light internal function to control on-off times; the bool on (if set false) overrides this like an external switch
	bool (*time)( double t ) = [] ( double t ){ return true; };

	// overloaded functions to check whether a point is inside, light is on, or both
	inline virtual bool check( void ) { return on; }
	inline virtual bool check( Point x ) { return ((invert ^ (*domain)(x)) && on); }
	inline virtual bool check( double t ) { return ((*time)(t) && on); }
	inline virtual bool check( Point x, double t ) { return ((invert ^ (*domain)(x)) && ((*time)(t) && on)); }
	// these functions are only needed in LED lattice
	inline virtual bool check( Point x, int m ) { return m==0 ? ((invert ^ (*domain)(x)) && on) : false; }
	inline virtual bool check( Point x, double t, int m ) { return m==0 ? ((invert ^ (*domain)(x)) && ((*time)(t) && on)) : false; };
};

struct Light_Type
{
	std::string type;
	std::vector<double> param;
	Light_Type() : param(2) {;}
	Light_Type( std::string _type, std::vector<double> _param ) : type{_type}, param{_param} {;}
};

class LED : public Light
{
private:
	Gauss gauss;
public:
	LED() { type = "LED"; spectrum = &gauss; }
	LED( std::string _type ) { type = _type; spectrum = &gauss; }
	LED( std::vector<double> _param ) : gauss(_param) { type = "LED"; spectrum = &gauss; }
	LED( std::string _type, std::vector<double> _param ) : gauss(_param) { type = _type; spectrum = &gauss; }
	LED( Light_Type definition ) : gauss(definition.param) { type = definition.type; spectrum = &gauss; }
	LED( const LED& led ); // copy constructor
};

class LASER : public Light
{
private:
	Delta delta;
public:
	LASER() { type = "LASER"; spectrum = &delta; }
	LASER( std::string _type ) { type = _type; spectrum = &delta; }
	LASER( std::vector<double> _param ) : delta(_param) { type = "LASER"; spectrum = &delta; }
	LASER( std::string _type, std::vector<double> _param ) : delta(_param) { type = _type; spectrum = &delta; }
	LASER( Light_Type definition ) : delta(definition.param) { type = definition.type; spectrum = &delta; }
	LASER( const LASER& laser ); // copy constructor
};

class BAND : public Light
{
private:
	Uniform uniform;
public:
	BAND() { type = "BAND"; spectrum = &uniform; }
	BAND( std::string _type ) { type = _type; spectrum = &uniform; }
	BAND( std::vector<double> _param ) : uniform(_param) { type = "BAND"; spectrum = &uniform; }
	BAND( std::string _type, std::vector<double> _param ) : uniform(_param) { type = _type; spectrum = &uniform; }
	BAND( Light_Type definition ) : uniform(definition.param) { type = definition.type; spectrum = &uniform; }
	BAND( const BAND& band ); // copy constructor
};

// combines several LEDs (e.g. different colors) and one (!) spatial illumination domain
class LEDs_Domain
{
private:
	std::string type = "LEDs_Domain";
	std::vector<std::string> led_types;
	std::map<std::string, int> led_index;
	int number_leds = 0;
public:
	LEDs_Domain() : leds(0), led_types(0), led_index() {;}
	void switch_on( int m = -1 ); // switch on LED m (starting with 0)
	void switch_off( int m = -1 ); // m=-1 results in switching ALL LEDs
	void toggle( int m = -1 );
	void set_intensity( double intensity, int m = -1 );
	double get_intensity( int m = 0 );
	double get_density( double w, int m = 0 );
	double get_density( Point x, double w, int m = 0 );
	double get_density( Point x, double t, double w, int m = 0 );
	std::string get_type( void ) const { return type; }

	void add_LED( LED& led );
	int get_number_leds( void ) { return number_leds; }
	std::string get_led_type( int m ) { return (m>=0&&m<number_leds) ? led_types[m] : "INDEX OUT OF RANGE"; }
	int get_led_index( std::string type ) { return led_index[type]; }

	std::vector<LED> leds;
	Rectangle domain;
};

// a lattice of LED domains
class LED_Lattice : public Light
{
private:
	bool is_there_light = false;
	bool set = false;
	std::vector<std::string> led_types;
	std::map<std::string, int> led_index;
	int number_leds = 0;
public:
	std::vector<LEDs_Domain> lattice;
	LED_Lattice();
	bool set_lattice( int n ); // creates a nxn lattice
	int get_lattice_size( void ) const { return lattice.size(); }

	void add_LED( std::string type, std::vector<double> param );
	void add_LED( Light_Type definition );
	int get_number_leds( void ) { return number_leds; }
	std::string get_led_type( int m ) { return (m>=0&&m<number_leds) ? led_types[m] : "INDEX OUT OF RANGE"; }
	int get_led_index( std::string type ) { return led_index[type]; }

	void switch_side_on( int n, int m = -1 ); // switch on LED m (starting with 0) on side n
	void switch_side_off( int n, int m = -1 ); // m=-1 results in switching ALL LEDs
	void switch_side_on( Point x, int m = -1 ); // switch on LED m (starting with 0) on site at point x
	void switch_side_off( Point x, int m = -1 ); // m=-1 results in switching ALL LEDs
	void toggle_side( int n, int m = -1 );
	void toggle_side( Point x, int m = -1 );
	void switch_all_sides_on( int m = -1 ); // switch ALL LEDs on
	void switch_all_sides_off( int m = -1 ); // switch ALL LEDs off

	void set_all_intensity( double intensity, int m = -1 ) override; // sets the intensity for all LEDs of type m
	void set_intensity( double intensity, int n, int m = -1  ) override; // sets the intensity at lattice side n
	void set_intensity( double intensity, Point x, int m = -1  ) override; // sets the intensity at point x
	double get_intensity( int n, int m = 0 ) override; // gets the intensity at lattice side n of LED type m
	double get_intensity( Point x, int m = 0 ) override; // sets the intensity at point x of LED type m

	// overloaded functions for the density value for the wavelength w
	double get_density( double w, int m = 0 ) override;
	double get_density( Point x, double w, int m = 0 ) override;
	double get_density( Point x, double t, double w, int m = 0 ) override;

	// find LED lattice sub-domain into which x falls and return the pointer
	// to the LEDs_Domain; returns NULL if x is outside the lattice
	LEDs_Domain* search_domain( Point x );

	// overloaded functions to check whether a point is inside, light is on, or both
	bool check( void ) override { return on; }
	bool check( Point x ) override;
	bool check( Point x, double t ) override;
	// functions specifc for the lattice to check whether there is light of LED type m
	bool check( Point x, int m ) override;
	bool check( Point x, double t, int m ) override;
};

#endif /* CUSTOM_MODULES_LIGHT_H_ */
