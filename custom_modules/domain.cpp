/*
 * domain.cpp
 *
 *  Created on: 17.01.2020
 *      Author: cfleck
 */

#include "domain.h"

std::vector < Point > Domain::get_nodes( ) const
{
	return nodes;
}

Default_Domain::Default_Domain()
{
	// Attention: the order of initialization matter; if default domain
	// is instantiated BEFORE the actual domain size is set from the XML file
	// on should call the set_up function again AFTER the XML file is read.
	type = "Rectangle";
	name = "DefaultDomain";
	min_number_nodes = 2;
	nodes.resize(4);
	set_up();
}

void Default_Domain::set_up( void )
{
	Point P0, P1;
	std::vector < Point > _nodes(2);
	number_nodes = 4;
	x_min = default_microenvironment_options.X_range[0];
	x_max = default_microenvironment_options.X_range[1];
	y_min = default_microenvironment_options.Y_range[0];
	y_max = default_microenvironment_options.Y_range[1];
	P0.set(x_min,y_max); P1.set(x_min,y_min);
	_nodes = {P0,P1};
	set_nodes( _nodes, x_max - x_min );
	set = true;
}

inline bool Default_Domain::point_inside( const Point& x )
{
	return true; // assumes that cells can only be inside the simulation domain
}

Rectangle::Rectangle() : AB(2), AD(2), AP(2)
{
	type = "Rectangle";
	min_number_nodes = 2;
	nodes.resize(4);
}

Rectangle::Rectangle( std::vector < Point > _nodes, double length, std::string _name ) : AB(2), AD(2), AP(2)
{
	type = "Rectangle";
	name = _name;
	min_number_nodes = 2;
	nodes.resize(4);
	number_nodes = _nodes.size();
	set_up( _nodes, length );
}

void Rectangle::set_up( std::vector < Point > _nodes, double length, std::string _name )
{
	name = _name;
	number_nodes = _nodes.size();
	if ( number_nodes < min_number_nodes )
	{
		set = false;
	}
	else
	{
		set_nodes( _nodes, length );
		set = true;
	}
	return;
}

void Rectangle::set_nodes( std::vector < Point > _nodes, double length )
{
	// construct the 2D vectors between the nodes
	AB = {_nodes[1].x - _nodes[0].x, _nodes[1].y - _nodes[0].y};
	AB_AB = std::inner_product(std::begin(AB),std::end(AB), std::begin(AB), 0.0);
	AD_AD = length*length;
	AD = {-AB[1]*sqrt(AD_AD/AB_AB), AB[0]*sqrt(AD_AD/AB_AB)};
//	std::cout << AB << std::endl;
//	std::cout << AD << std::endl;

	nodes[0] = _nodes[0]; nodes[1] = _nodes[1];
	nodes[2] = _nodes[1]; nodes[2].move(AD);
	nodes[3] = _nodes[0]; nodes[3].move(AD);

//	std::cout << nodes[0].tuple() << std::endl;
//	std::cout << nodes[1].tuple() << std::endl;
//	std::cout << nodes[2].tuple() << std::endl;
//	std::cout << nodes[3].tuple() << std::endl;

	// find the max and min values
	std::vector<double> xval;
	std::vector<double> yval;

	xval = {nodes[0].x, nodes[1].x, nodes[2].x, nodes[3].x};
	yval = {nodes[0].y, nodes[1].y, nodes[2].y, nodes[3].y};

	x_min = *std::min_element(std::begin(xval), std::end(xval));
	x_max = *std::max_element(std::begin(xval), std::end(xval));
	y_min = *std::min_element(std::begin(yval), std::end(yval));
	y_max = *std::max_element(std::begin(yval), std::end(yval));

	// calculate the area
	area = sqrt( AB_AB*AD_AD );

	// determine center of the rectangle
	center.set( nodes[0].x + 0.5*(AB[0]+AD[0]), nodes[0].y + 0.5*(AB[1]+AD[1]) );
}

// are the boundaries a problem?
inline bool Rectangle::point_inside( const Point& P )
{
	bool inside;

//	std::cout << "Checking whether point P is inside the rectangle" << std::endl;

	//check whether domain is defined
	if(!set) { inside = false; return inside; }

	// first simple test
	if( P.x < x_min || P.x > x_max || P.y < y_min || P.y > y_max )
	{
		inside = false;
	}
	else
	{
		// the point P lies inside the rectangle iff
		// (0<AP*AB < ||AB||^2) && (0<AP*AD < ||AD||^2)

		AP = {P.x - nodes[0].x, P.y - nodes[0].y};
	//	std::cout << std::to_string(AP[0]) + ", " + std::to_string(AP[1]) << std::endl;

	//	std::cout << "Calculating scalar products" << std::endl;
		double AP_AB = std::inner_product(std::begin(AP),std::end(AP), std::begin(AB), 0.0);
	//	std::cout << "AP*AB = " << AP_AB << std::endl;
		double AP_AD = std::inner_product(std::begin(AP),std::end(AP),std::begin(AD), 0.0);
	//	std::cout << "AP*AD = " << AP_AD << std::endl;

		if( ((0<AP_AB) && (AP_AB<AB_AB)) && ((0<AP_AD) && (AP_AD<AD_AD) ) )
		{
			inside = true;
		}
		else
		{
			inside = false;
		}
	}
	return inside;
}

void Rectangle::move( std::vector<double> translation )
{
	if( translation.size() > 1 )
	{
		for( auto iter = nodes.begin(); iter != nodes.end(); ++iter )
		{
			iter->x += translation[0];
			iter->y += translation[1];
		}
//		 find the max and min values
		std::vector<double> xval;
		std::vector<double> yval;

		xval = {nodes[0].x, nodes[1].x, nodes[2].x, nodes[3].x};
		yval = {nodes[0].y, nodes[1].y, nodes[2].y, nodes[3].y};

		x_min = *std::min_element(std::begin(xval), std::end(xval));
		x_max = *std::max_element(std::begin(xval), std::end(xval));
		y_min = *std::min_element(std::begin(yval), std::end(yval));
		y_max = *std::max_element(std::begin(yval), std::end(yval));

		// determine the new center of the rectangle
		center.set( nodes[0].x + 0.5*(AB[0]+AD[0]), nodes[0].y + 0.5*(AB[1]+AD[1]) );
	}
}

void Rectangle::rotate( double angle )
{
	if(angle == 0){ return; }

	angle = angle*2.0*PhysiCell::PhysiCell_constants::pi/360.0;

	std::vector < Point > _nodes(2);
	std::vector<double> v1(2), v2(2);

	v1 = {(nodes[0].x - center.x)*cos(angle) - (nodes[0].y - center.y)*sin(angle),
			(nodes[0].x - center.x)*sin(angle) + (nodes[0].y - center.y)*cos(angle)};
	v2 = {(nodes[1].x - center.x)*cos(angle) - (nodes[1].y - center.y)*sin(angle),
			(nodes[1].x - center.x)*sin(angle) + (nodes[1].y - center.y)*cos(angle)};

	_nodes[0].set(v1[0]+center.x, v1[1]+center.y);
	_nodes[1].set(v2[0]+center.x, v2[1]+center.y);
//	std::cout << _nodes[0].tuple() << std::endl;
//	std::cout << _nodes[1].tuple() << std::endl;
	set_nodes( _nodes, sqrt(AD_AD) );
}

Circle::Circle() : CP(2)
{
	type = "Circle";
	min_number_nodes = 1;
	nodes.resize(1);
	radius = 0; r2 =0;
}

Circle::Circle( Point _center, double _radius, std::string _name ) : CP(2)
{
	type = "Circle";
	name = _name;
	min_number_nodes = 1;
	nodes.resize(1);
	number_nodes = 1;;
	set_up( _center, radius );
}

void Circle::set_up( Point _center, double _radius, std::string _name )
{
	name = _name;
	radius = _radius;
	r2 = radius*radius;
	area =r2*PhysiCell::PhysiCell_constants::pi;
	set_nodes( _center );
	set = true;
}

void Circle::set_nodes( Point _center  )
{
	nodes[0] = _center;
	center = _center;
}

void Circle::move( std::vector<double> translation )
{
	if( translation.size() > 1 )
	{
		nodes[0].x += translation[0];
		nodes[0].y += translation[1];
		center = nodes[0];
	}
}

inline bool Circle::point_inside( const Point& P )
{
	bool inside;
//	Point p = P;
//	std::cout << "Checking whether point " << p.tuple() << " is inside the circle" << std::endl;
	//check whether domain is defined
	if(!set) {
		inside = false;
		return inside;
	}
	// vector from center of the circle to point P
	CP = {P.x - nodes[0].x, P.y - nodes[0].y};

//	std::cout << std::to_string(CP[0]) + ", " + std::to_string(CP[1]) << std::endl;

	if( std::inner_product(std::begin(CP),std::end(CP), std::begin(CP), 0.0) <= r2 )
	{
		inside = true;
	}
	else
	{
		inside = false;
	}

	return inside;
}

// Attention: the Polygon class is not finished yet. Not tested.
Polygon::Polygon( )
{
	min_number_nodes = 3;
	number_nodes = 0;
	set = false;
}

Polygon::Polygon( std::vector < Point > _nodes, std::string _name )
{
	name = _name;
	min_number_nodes = 3;
	number_nodes = _nodes.size();
	if ( number_nodes < 3 )
	{
		set = false;
	}
	else
	{
		nodes = _nodes;
		set = true;
	}
}

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"



inline int Polygon::isLeft( Point P0, Point P1, Point P2 )
{
    return ( (P1.x - P0.x) * (P2.y - P0.y)
            - (P2.x -  P0.x) * (P1.y - P0.y) );
}

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
int Polygon::cn_PnPoly( Point P, Point* V, int n )
{
    int    cn = 0;    // the  crossing number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {    // edge from V[i]  to V[i+1]
       if (((V[i].y <= P.y) && (V[i+1].y > P.y))     // an upward crossing
        || ((V[i].y > P.y) && (V[i+1].y <=  P.y))) { // a downward crossing
            // compute  the actual edge-ray intersect x-coordinate
            float vt = (float)(P.y  - V[i].y) / (V[i+1].y - V[i].y);
            if (P.x <  V[i].x + vt * (V[i+1].x - V[i].x)) // P.x < intersect
                 ++cn;   // a valid crossing of y=P.y right of P.x
        }
    }
    return (cn&1);    // 0 if even (out), and 1 if  odd (in)

}


// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
int Polygon::wn_PnPoly( Point P, Point* V, int n )
{
    int    wn = 0;    // the  winding number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {   // edge from V[i] to  V[i+1]
        if (V[i].y <= P.y) {          // start y <= P.y
            if (V[i+1].y  > P.y)      // an upward crossing
                 if (isLeft( V[i], V[i+1], P) > 0)  // P left of  edge
                     ++wn;            // have  a valid up intersect
        }
        else {                        // start y > P.y (no test needed)
            if (V[i+1].y  <= P.y)     // a downward crossing
                 if (isLeft( V[i], V[i+1], P) < 0)  // P right of  edge
                     --wn;            // have  a valid down intersect
        }
    }
    return wn;
}

void Domain_Union::add_domain( Domain* _domain, bool _invert )
{
	if( _domain == NULL ) { return; }
	domains.push_back(_domain);
	invert.push_back(_invert);
	// calculation of the total area needed
	return;
}

bool Domain_Union::point_inside( const Point& x )
{
	if( domains.size() == 0 ){ return false; }

	bool inside = false;
	auto in = invert.begin();
	for( auto iter = domains.begin(); iter < domains.end(); ++iter )
	{
		inside = (inside || ((*in) ^ (*iter)->point_inside(x)));
		++in;
	}
	return inside;
}

void Domain_Intersection::add_domain( Domain* _domain, bool _invert )
{
	if( _domain == NULL ) { return; }
	domains.push_back(_domain);
	invert.push_back(_invert);
	// calculation of the area needed
	return;
}

bool Domain_Intersection::point_inside( const Point& x )
{
	if( domains.size() == 0 ){ return false; }

	bool inside = true;
	auto in = invert.begin();
	for( auto iter = domains.begin(); iter < domains.end(); ++iter )
	{
		inside = (inside && ((*in) ^ (*iter)->point_inside(x)));
		++in;
	}
	return inside;
}
