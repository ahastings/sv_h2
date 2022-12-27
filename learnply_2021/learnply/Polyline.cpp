#include "Polyline.h"
#include "GL/glew.h"//Then we want to create a display function
#define EPSILON 1.0e-5

bool POLYLINE::isNeighbor(const POLYLINE& line) {
	if ((m_vertices.front() - line.m_vertices.front()).length() < EPSILON ||
		(m_vertices.front() - line.m_vertices.front()).length() < EPSILON ||
		(m_vertices.front() - line.m_vertices.front()).length() < EPSILON ||
		(m_vertices.front() - line.m_vertices.front()).length() < EPSILON) {
		return true;
	}else{
		return false;
	}
	
}
void POLYLINE::merge(const POLYLINE& line) {
	if ((m_vertices.front() - line.m_vertices.front()).length() < EPSILON) {
		POLYLINE line_ = line;
		line_.m_vertices.pop_front();
		for (auto i = line_.m_vertices.begin(); i != line_.m_vertices.end(); i++) {
			m_vertices.push_front(*i);
		}

	}
	else if ((m_vertices.front() - line.m_vertices.back()).length() < EPSILON) {
		POLYLINE reverseLine = line;
		reverseLine.m_vertices.pop_back();
		reverseLine.m_vertices.reverse();
		for (auto i = reverseLine.m_vertices.begin(); i != reverseLine.m_vertices.end(); i++) {
			m_vertices.push_front(*i);
		}

	}
	else if ((m_vertices.back() - line.m_vertices.front()).length() < EPSILON) {
		POLYLINE line_ = line;
		line_.m_vertices.front();
		for (auto i = line_.m_vertices.begin(); i != line_.m_vertices.end(); i++) {
			m_vertices.push_front(*i);
		}

	}
	else if ((m_vertices.back() - line.m_vertices.back()).length() < EPSILON) {
		POLYLINE reverseLine = line;
		reverseLine.m_vertices.pop_back();
		reverseLine.m_vertices.reverse();
		for (auto i = reverseLine.m_vertices.begin(); i != reverseLine.m_vertices.end(); i++) {
			m_vertices.push_front(*i);
		}

	}
	return;//true;
}
void display_polyline(std::vector<POLYLINE>& polylines) {
	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	for (auto& polyline : polylines) {
		glLineWidth(polyline.m_weight);
		glColor3f(polyline.m_rgb.entry[0], polyline.m_rgb.entry[1], polyline.m_rgb.entry[2]);
		glBegin(GL_LINE_STRIP);
		for (auto it = polyline.m_vertices.begin(); it != polyline.m_vertices.end(); it++) {
			glVertex3d(it->entry[0], it->entry[1], it->entry[2]);
		}
		glEnd();
	}
	glDisable(GL_BLEND);
	glLineWidth(1); // width of line drawn
}


void marchingSqure(std::list<POLYLINE>& edges,
	const Polyhedron& poly,
	const double& thres) {

}


void lookUpTable(
	std::vector<Vertex>& r,
	const Vertex& v0,
	const Vertex& v1,
	const Vertex& v2,
	const Vertex& v3,
	const double& thres){

}


void marchingSquare(
	std::list<POLYLINE>& edges, const Polyhedron& poly,
	const double& thres) {
	for (int i = 0; i < poly.nquads; i++) {
		std::vector<Vertex> r;
		lookUpTable(r,
			*poly.qlist[i]->verts[0],
			*poly.qlist[i]->verts[1],
			*poly.qlist[i]->verts[2],
			*poly.qlist[i]->verts[3],
			thres);
		if (r.size() > 0) {
			for (int r_i = 0; r_i < r.size() / 2; r_i++) {
				POLYLINE line;
				auto v0 = icVector3(
					r[r_i * 2].x,
					r[r_i * 2].y,
					r[r_i * 2].z);
				auto v1 = icVector3(
					r[r_i * 2 + 1].x,
					r[r_i * 2 + 1].y,
					r[r_i * 2 + 1].z);
				/*if((ve - v1).length() < EPSILON)
				------------
				line.m_vertices.push_back(v®);
				else {
				line.m_vertices.push_back(ve);
				line.m_vertices.push_back(v1);
				}*/
				line.m_vertices.push_back(v0);
				line.m_vertices.push_back(v1);
				edges.push_back(line);
			}
		}

	}
}


void makePolylineFromEdges(
	std::vector<POLYLINE>& polylines,
	const std::list <POLYLINE>& edges){

}
