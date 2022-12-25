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
