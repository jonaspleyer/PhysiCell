/*
 * spectrum.h
 *
 *  Created on: 17.01.2020
 *      Author: cfleck
 */

#ifndef CUSTOM_MODULES_SPECTRUM_H_
#define CUSTOM_MODULES_SPECTRUM_H_

#include "domain.h"

// For debugging
#ifndef _HI_
#define _HI_ std::cout << "HERE I AM" << std::endl;
#endif


// TODO: provide a general spectrum class in which the spectrum is discretely or by a function defined
// TODO: find a way to combine the spectra of light sources.

class Spectral_Density
{
protected:
	std::string type;
	std::vector<double> param;
public:
	Spectral_Density() : param(2) {};
	virtual ~Spectral_Density() {};
	virtual void set_param( std::vector<double> _param ) = 0;
	virtual double operator()( double w ) = 0;
	std::string get_type( void ) { return type; }
	std::vector<double> get_param( void ) { return param; }
};

class Delta : public Spectral_Density
{
private:
	double width = 0.05; // band width of the approximate delta distribution in nm
	double norm = 2.0 * width;
public:
	Delta() { type = "Delta"; }
	Delta( std::vector<double> _param ) { type = "Delta"; param = _param; }
	void set_param( std::vector<double> _param ) { param = _param; }
	inline double operator()( double w ) { return ((w >= param[0] - width) && (w <= param[0] + width)) ? 1.0/norm : 0.0; }
};

class Uniform : public Spectral_Density
{
private:
	double norm;
public:
	Uniform() { type = "Uniform"; norm = 0;};
	Uniform( std::vector<double> _param ) { type = "Uniform"; param = _param; norm = param[0] - param[1]; }
	void set_param( std::vector<double> _param ) { param = _param; norm = param[0] - param[1]; }
	inline double operator()( double w ) { return ((w >= param[0]) && (w <= param[1])) ? 1.0/norm : 0.0; }
};

class Gauss : public Spectral_Density
{
private:
	double norm;
public:
	Gauss() { type = "Gauss"; norm = 0; }
	Gauss( std::vector<double> _param ) { type = "Gauss"; param = _param; norm = sqrt(2.0*PhysiCell::PhysiCell_constants::pi)*param[1]; }
	void set_param( std::vector<double> _param ) { param = _param; norm = sqrt(2.0*PhysiCell::PhysiCell_constants::pi)*param[1]; }
	inline double operator()( double w ) { return exp( -pow(w - param[0],2.0)/(2.0*param[1]*param[1]) )/norm; }
};

#endif /* CUSTOM_MODULES_SPECTRUM_H_ */
