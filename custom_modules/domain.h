/*
 * domain.h
 *
 *  Created on: 17.01.2020
 *      Author: cfleck
 */

#ifndef CUSTOM_MODULES_DOMAIN_H_
#define CUSTOM_MODULES_DOMAIN_H_

#include "../core/PhysiCell.h"

#include <vector>
#include <algorithm>

// For debugging
#ifndef _HI_
#define _HI_ std::cout << "HERE I AM" << std::endl;
#endif

// For a 2D point
// ToDo: extend to 3D
struct Point
{
	double x;
	double y;
	void operator=(const Point& P){x = P.x; y = P.y; }
	void set(double _x, double _y){x = _x; y = _y; }
	void move( std::vector<double> translation ){ x += translation[0]; y += translation[1]; }
	std::string tuple( void ) { return "(" + std::to_string(x) + ", " + std::to_string(y) + ")"; }
};

class Domain
{
 protected:
	std::string type = "Basic";
	std::string name = "NoName";
	bool set = false;
	int number_nodes = 0;
	int min_number_nodes = 1;
	double area = 0;
	std::vector < Point > nodes;
	Point center; // for the center of mass

 public:
	Domain() : nodes(1) {;}
	Domain( std::vector < Point > _nodes, std::string _name );
	virtual ~Domain() {};
	Point get_center( void ) const { return center; }
	virtual bool point_inside( const Point& x ) = 0;
	virtual void move( std::vector<double> translation ) = 0;
	virtual void rotate( double angle ) = 0; // rotate around the center
	inline bool operator()( const Point& x ) { return point_inside(x); }
	std::vector < Point > get_nodes ( void ) const;
	std::string get_type( void ) const { return type; }
	void set_name( std::string _name ){ name = _name; }
	std::string get_name( void ) const { return name; }
	double get_area( void ) const { return area; }
};

// The rectangle is defined by giving two points and a length. The two points define the vector AB and the lenght the length of the
// orthogonal vector AD
class Rectangle : public Domain
{
protected:
	void set_nodes( std::vector < Point > _nodes, double length );

	// used for a first simple check
	double x_min = 0;
	double x_max = 0;
	double y_min = 0;
	double y_max = 0;

 public:
	// vectors needed for the point inside rectangle test
	// if the rectangle is moved, they need to be updated
	std::vector<double> AB;
	std::vector<double> AD;
	std::vector<double> AP;

	double AB_AB = 0;
	double AD_AD = 0;

	Rectangle();
	Rectangle( std::vector < Point > _nodes, double lenght, std::string _name = "NoName" );

	void set_up( std::vector < Point > _nodes, double lenght, std::string _name = "NoName" );

	void move( std::vector<double> translation );
	void rotate( double angle );
	bool point_inside( const Point& x );
};

class Default_Domain : public Rectangle
{
private:

public:
	Default_Domain();
	void set_up( void );
	bool point_inside( const Point& x ) override;
	void move( std::vector<double> translation ) override { ; }
	void rotate( double angle ) override { ; }
};


class Circle : public Domain
{
 private:
	void set_nodes( Point _center  );
	double radius;
	double r2;

	// vectors needed for the point inside rectangle test
	// if the rectangle is moved, they need to be updated
	std::vector<double> CP;

 public:
	Circle();
	Circle( Point _center, double radius, std::string _name = "NoName" );
	void set_up( Point _center, double _radius, std::string _name = "NoName" );
	void move( std::vector<double> translation );
	void rotate( double angle ){ ; }
	bool point_inside( const Point& x );
};

class Polygon : public Domain
{
 private:
	bool set = false;
	std::vector < Point > nodes;

	int isLeft( Point P0, Point P1, Point P2 );
	int cn_PnPoly( Point P, Point* V, int n );
	int wn_PnPoly( Point P, Point* V, int n );

 public:
	Polygon( );
	Polygon( std::vector < Point > _nodes, std::string _name = "NoName" );

	// to be implemented later
	void move( std::vector<double> translation ) {;}
	void rotate( double angle ) { ; }
//	bool point_inside( Point x );
};

class Domain_Union : public Domain
{
private:
	std::vector<Domain*> domains;
	std::vector<bool> invert;
public:
	Domain_Union( std::string _name = "NoName" ) : domains(0)  { name = _name; type = "DomainUnion"; }
	void add_domain( Domain* _domain, bool _invert = false );

	// to be removed later
	void set_area( double _area ){ area = _area; }

	// to be implemented later
	void remove_domain( int n ){;}

	// to be implemented later
	void move( std::vector<double> translation ) {;}
	void rotate( double angle ) {;}

	bool point_inside( const Point& x );
};

class Domain_Intersection : public Domain
{
private:
	std::vector<Domain*> domains;
	std::vector<bool> invert;
public:
	Domain_Intersection( std::string _name = "NoName" ) : domains(0) { name = _name; type = "DomainIntersection"; }
	void add_domain( Domain* _domain, bool _invert = false );

	// to be removed later
	void set_area( double _area ){ area = _area; }

	// to be implemented later
	void remove_domain( int n ) {;}

	// to be implemented later
	void move( std::vector<double> translation ) {;}
	void rotate( double angle ) {;}

	bool point_inside( const Point& x );
};

#endif /* CUSTOM_MODULES_DOMAIN_H_ */
