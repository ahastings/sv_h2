#include "Polyline.h"
//#include "polyhedron.h"
//Polyhedron* poly;
#include "GL/glew.h"//Then we want to create a display function
#define EPSILON 1.0e-5
#include "iostream"
#include "gl/freeglut.h"
#include "gl/glew.h"
extern Polyhedron* poly;

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


/*void marchingSqure(std::list<POLYLINE>& edges,
	const Polyhedron& poly,
	const double& thres) {

}*/

Vertex linearInterpolateByScalar(const Vertex& v0, const Vertex& v1, const double& thres) {
	double f0_f1 = v0.scalar - v1.scalar;
	Vertex r(0.0, 0.0, 0.0);
	if (std::abs(f0_f1) < 1.0e-5) {
		r.x = (v0.x + v1.x) / 2;
		r.y = (v0.y + v1.y) / 2;
		r.z = (v0.z + v1.z) / 2;
	}
	else {
		double t = std::abs((v0.scalar - thres) / ((v0.scalar - thres) - (v1.scalar - thres)));
		r.x = v0.x + t * (v1.x - v0.x);
		r.y = v0.y + t * (v1.y - v0.y);
		r.z = v0.z + t * (v1.z - v0.z);
	}
	return r;
}


void lookUpTable(
	std::vector<Vertex> &r,
	const Vertex& v0,
	const Vertex& v1,
	const Vertex& v2,
	const Vertex& v3,
	const double& thres) {

	r.reserve(2);
	int id = 0;
	if (v0.scalar <= thres + EPSILON) {
		id += 1;
	}
	if (v1.scalar <= thres + EPSILON) {
		id += 2;
	}
	if (v2.scalar <= thres + EPSILON) {
		id += 4;
	}
	if (v3.scalar <= thres + EPSILON) {
		id += 8;
	}
	double center = 0;
	switch (id)
	{
	case 0:
		break;
	case 1:
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 2:
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		break;
	case 3:
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 4:
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		break;
	case 5:
		center = v0.scalar, v1.scalar, v2.scalar, v3.scalar;
		center /= 4;
		if (center <= thres) {
			r.push_back(linearInterpolateByScalar(v0, v1, thres));
			r.push_back(linearInterpolateByScalar(v1, v2, thres));
			r.push_back(linearInterpolateByScalar(v2, v3, thres));
			r.push_back(linearInterpolateByScalar(v0, v3, thres));
		}
		else {
			r.push_back(linearInterpolateByScalar(v0, v1, thres));
			r.push_back(linearInterpolateByScalar(v0, v3, thres));
			r.push_back(linearInterpolateByScalar(v1, v2, thres));
			r.push_back(linearInterpolateByScalar(v2, v3, thres));
		}
		break;
	case 6://same as 9
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		break;
	case 7:
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 8:
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 9: //same as 6
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
		break;
	case 10:
		center = v0.scalar + v1.scalar + v2.scalar + v3.scalar;
		center /= 4;
		if (center <= thres) {
			r.push_back(linearInterpolateByScalar(v0, v1, thres));
			r.push_back(linearInterpolateByScalar(v0, v3, thres));
			r.push_back(linearInterpolateByScalar(v1, v2, thres));
			r.push_back(linearInterpolateByScalar(v2, v3, thres));
		}
		else {
			r.push_back(linearInterpolateByScalar(v0, v1, thres));
			r.push_back(linearInterpolateByScalar(v1, v2, thres));
			r.push_back(linearInterpolateByScalar(v2, v3, thres));
			r.push_back(linearInterpolateByScalar(v0, v3, thres));
		}break;
	case 11:
		// v1 - v2 //v2- v3
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		r.push_back(linearInterpolateByScalar(v2, v3, thres));
	break; case 12:
		// v2 v2 //ve- v3
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
	break; case 13:
		// v0 - vi //v1- v2
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v1, v2, thres));
		break;
	case 14:
		// v0 - v1 //v0- v3
		r.push_back(linearInterpolateByScalar(v0, v1, thres));
		r.push_back(linearInterpolateByScalar(v0, v3, thres));
		break;
	case 15:
		break;
		

	}
}

void marchingSqure(
	std::list<POLYLINE>& edges,
	const Polyhedron& poly,
	const double& thres) {
	for (int i = 0; i < poly.nquads; i++) {
		std::vector<Vertex>r;
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
	polylines.reserve(edges.size());
	std::list<POLYLINE> edges_temp(edges);
	while (edges_temp.size() > 0) {
		polylines.push_back(edges_temp.front());
		edges_temp.erase(edges_temp.begin());
		int init_size = 0;
		while (init_size != edges_temp.size()) {
			init_size = edges_temp.size();
			for (auto i = edges_temp.begin(); i != edges_temp.end();) {
				if (polylines.back().isNeighbor(*i)) {
					polylines.back().merge(*i);
					i = edges_temp.erase(i);
				}
				else
					i++;
			}
		}
	}
}
//Polyhedron* poly1;
void findMm(double& M, double& m) {
	m = INFINITY;
	M = -m;
	for (int i = 0; i < poly->nverts; i++) {

		auto& vertex = poly->vlist[i];

		if (vertex->scalar < m) {
			m = vertex->scalar;
		}
		if (vertex->scalar > M) {
			M = vertex->scalar;
		}
	}

}
void project1() {
	std::cout << "Assignment color for a scalar value" << std::endl;
	double m = INFINITY;
	double M = -m;

	findMm(M, m);
	std::cout << "Min: " << m << std::endl;
	std::cout << "Max: " << M << std::endl;

	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		double s_v = vertex->scalar;
		double gray_ = (s_v - m) / (M - m);

		vertex->R = vertex->G = vertex->B = gray_; // Graying the colors r g b

	}
	glutPostRedisplay();
}
void project1b() {
	std::cout << "Assignment color for a rgb value" << std::endl;


	double m = INFINITY;
	double M = -m;
	findMm(M, m);
	double Mminusm = M - m;

	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];

		double s_v = vertex->scalar;
		double l = (s_v - m) / (M - m);
		double r = (M - s_v) / (M - m);

		icVector3 c1(1.0, 0.0, 0.0);
		icVector3 c2(0.0, 0.0, 1.0);
		//vertex->R = C1R * 1 + C2R * r;
		icVector3 c = c1 * 1 + c2 * r;

		vertex->R = c.x;
		vertex->G = c.y;
		vertex->B = c.z;

		//vertex->R = vertex->G = vertex->B = gray_; // Graying the colors r g b

	}
	glutPostRedisplay();
}


void HSVtoRGB(
	icVector3& rgb,
	const icVector3& hsv) {
	double H = hsv.x;
	double S = hsv.y;
	double v = hsv.z;
	double C = S * v;
	double X = C * (1 - abs(fmod(H / 60.0, 2) - 1));
	double m = v - C;
	double r, g, b;
	if (H >= 0 && H < 60) {
		r = C, g = X, b = 0;
	}
	else if (H >= 60 && H < 120) {
		r = X, g = C, b = 0;
	}
	else if (H >= 120 && H < 180) {
		r = 0, g = C, b = X;
	}
	else if (H >= 180 && H < 240) {
		r = 0, g = X, b = C;
	}
	else if (H >= 240 && H < 300) {
		r = X, g = 0, b = C;
	}
	else {
		r = C, g = 0, b = X;
	}
	rgb.x = (r + m);
	rgb.y = (g + m);
	rgb.z = (b + m);
}

void RGBtoHSV(
	icVector3& hsv,
	const icVector3& rgb) {
	double r = rgb.x;
	double g = rgb.y;
	double b = rgb.z;
	// h , s, v = hue, saturation, value
	double cmax = std::max(r, std::max(g, b)); // Max of r, g, b
	double cmin = std::min(r, std::min(g, b)); // Min of r, g, b
	double diff = cmax - cmin; // diff of cmax and min

	double& h = hsv.x;
	double& s = hsv.y;
	double& v = hsv.z;
	// if cmax equal to cmin then h = 0
	if (cmax == cmin)
		h = 0;
	// if cmax equal to r then compute h
	else if (cmax == r)
		h = fmod(60 * ((g - b) / diff) + 360, 360);
	// if cmax equal to g then compute h
	else if (cmax == g)
		h = fmod(60 * ((b - r) / diff) + 120, 360);
	// if cmax equal to b then compute h
	else if (cmax == b)
		h = fmod(60 * ((r - g) / diff) + 240, 360);

	//if max is zero -- edge/base case
	if (cmax == 0)
		s = 0;
	else
		s = (diff / cmax);

	//compute v
	v = cmax;
}

void project1c() {
	std::cout << "C" << std::endl;
	double m;
	double M;
	findMm(M, m);
	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		double s_v = vertex->scalar;
		icVector3 c1(1.0, 0.0, 0.0);
		icVector3 c2(0.0, 0.0, 1.0);
		icVector3 HSVc1, HSVc2;
		RGBtoHSV(HSVc1, c1);
		RGBtoHSV(HSVc2, c2);
		double l = (s_v - m) / (M - m);
		double r = (M - s_v) / (M - m);
		icVector3 HSVc = HSVc1 * l + HSVc2 * r;
		icVector3 RGBc;
		HSVtoRGB(RGBc, HSVc);
		vertex->R = RGBc.x;
		vertex->G = RGBc.y;
		vertex->B = RGBc.z;

	}
	glutPostRedisplay();
}


void project1_2() {
	std::cout << "Height field method" << std::endl;

	double m;
	double M;
	findMm(M, m);
	for (int i = 0; i < poly->nverts; i++) {
		auto& vertex = poly->vlist[i];
		double s_v = vertex->scalar;
		double l = (s_v - m) / (M - m);
		vertex->z = 10 * l;

	}
	glutPostRedisplay();
}

