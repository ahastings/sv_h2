#include "Polyline.h"
#include "GL/glew.h"//Then we want to create a display function
bool POLYLINE::isNeighbor(const POLYLINE& line) {
	//if ((m_vertices.front() - line.m_vertices.)
	return true;
}
bool POLYLINE::merge(const POLYLINE& line) {
	return true;
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
