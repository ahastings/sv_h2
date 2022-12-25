#pragma once
#include "icVector.H"
#include <vector>
#include <list>
#include <utility>


class POLYLINE {
public:
	std::list<icVector3> m_vertices;
	icVector3 m_rgb = icVector3(1.0, 0.0, 0.0);
	double m_weight = 1.0;

	bool isNeighbor(const POLYLINE& line);
	void merge(const POLYLINE& line);//was  bool merge(const POLYLINE& line);

private:

};

void display_polyline(std::vector<POLYLINE>& polylines);

class LineSegment
{
public:

	// fields
	icVector3 start, end;
	double len;

	// constructors

	LineSegment(icVector3 start_in, icVector3 end_in)
	{
		start = start_in;
		end = end_in;
		len = length(end - start);
	}

	LineSegment(double sx, double sy, double sz, double ex, double ey, double ez)
	{
		start = icVector3(sx, sy, sz);
		end = icVector3(ex, ey, ez);
		len = length(end - start);
	}

	// methods

	icVector3 midpoint()
	{
		icVector3 diff = end - start;
		return start + (0.5 * diff);
	}
};

// PolyLine is a list of connected line segments
typedef std::vector<LineSegment> PolyLine;