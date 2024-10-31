// 31/10/2024	 	Fix bug in stl_edge::dump
// 25/02/2023	1.32	Fix bug in -scalemm.
// 06/07/2019	1.31	Correct file reading problems.
// 13/04/2019	1.30	Adaptation to visual studio 2015. -scale, -scalemm, -verbose added. -silent by default.
// 11/07/2011	1.22	Remove trace in nodebug mode.
// 03/05/2011	1.21	Correct bug that stopped reading of file when very small facet was encountered.
// 03/05/2011	1.20	Add option -oeps used for remove unuseful facets, Debug : Use smaller eps when reading the facets
// 31/01/2011	1.19	Add option -raw (conversion without any optimisation or edge calculation)
// 19/01/2011	1.18	Add option -out (output file name)
// 18/01/2011	1.17	Correct bug of Author and partname was wrong when -o option is present
// 15/12/2010	1.16	Suppress define DEBUG
// 20/11/2007	1.15	Allow spaces in description in header
// 21/02/2007	1.14	Conform to new header specs http://www.ldraw.org/article/398
// 23/04/2005	1.13	Add option -deps (max determinant) and test non coplanar quads
// 23/02/2003	1.12	Add option -o5 to force creation of all optional lines
// 18/07/2002	1.01	Add options -no1 to -no5 to exclude line type
// 14/08/2001	b10	Use two or three cylinders primitives when the plain primitive does not exist (use new 1-16cyli and 1-16edge)
// 12/11/2000	b9	Improve detection of cylinders by checking the height 
// 01/09/2000	b8	Improve detection of unuseful edges by counting the different normal of faces surrounding the edge
// 02/05/2000	b7	Remove unuseful edges (experiment), transformation matrix
// 13/03/2000	b6	Reads binary stl files. Set origin
// 13/03/2000	b5	Detects cylinder primitives by basing only on angle between facet
// 18/02/2000	b4	Scan edges at begining, topologic structure (adjency)
// 14/02/2000	b3	Type 2 lines corrected, merge3 called before calculate_edge
// 10/02/2000	b2	Creates now standard primitives n-4cyli.dat and n-4edge.dat
// 28/01/2000	b1	First release on internet
// stl2dat.cpp : Defines the entry point for the console application.
//
//#ifdef _DEBUG
//#define DEBUG
//#endif

//#include "stdafx.h"
#include "math.h"
#include "fstream"
#include "sstream"
#include "iomanip"
//#include "io.h"
#include "unistd.h"

#include "vect.h"
#include "ort.h"
#include "c_reel.h"
#include "c_sea.h"

#include "stl2dat.h"
#include <string.h>
#include <stdio.h>
#include <algorithm>

using namespace std;

#define version "1.32"

static double ag_lim;
static double ag_lim_q;
static double ag_lim_t;
static double determinant_eps;
static bool cyl_norm;
static bool no_plane;
static bool no_cylinder;
static bool no_tangent;
static bool opt5;
static int col_line1;
static int col_line2;
static int col_line3;
static int col_line4;
static int col_line5;
static bool no_line1;
static bool no_line2;
static bool no_line3;
static bool no_line4;
static bool no_line5;
static bool rand_col_surfaces;
static bool rand_col_adj_cyl;
static bool ldr_opt;
static bool ldr_out;
static bool bfc;
static bool silent;
static bool ldraw;
static bool stl_raw;
static double optim_eps;

static bool noprim;
static bool debug;
static bool nodebug;

//#define DEBUG

istream & operator >> (istream & is, string & string)
{
	char s[256];
	if (is)
	{
		is >> s;
		string = s;
	}
	return is;
}

double ftrunc(double v, int precision);

stl_v stl_edge::coord(int i, int face)
{
	if (!adjacent[face])
		face = 1 - face;
	if (!adjacent[face])
		return stl_v_zero;

	int usage = adjacent[face]->usage;
	if (usage == 0)
		usage = 4;
	return(adjacent[face]->coords[(coords[face] + i) % usage]);
}

stl_vertex * stl_edge::vertx(int i, int face)
{
	int usage = adjacent[face]->usage;
	if (usage == 0)
		usage = 4;
	return(adjacent[face]->vertex[(coords[face] + i) % usage]);
}

void stl_edge::set_coord(int i, int face, const stl_v & val)
{
	if (adjacent[face])
	{
		int usage = adjacent[face]->usage;
		if (usage == 0)
			usage = 4;
		adjacent[face]->coords[(coords[face] + i) % usage] = val;
	}
}

bool stl_edge::colinear(const stl_v & vd)
{
	return(coord(1) - coord(0)).paral(vd);
}

stl_facet::stl_facet()
{
	usage = 3;
	color = 16;
	surf = NULL;
	for (int i = 0; i < 4; i++)
	{
		edge[i] = 0;
		vertex[i] = 0;
		adjacent[i] = 0;
	}
	printed = false;
	m_surface = -1;
	precision = 2;
}

static long number_edge = 0;

stl_edge::stl_edge()
{
	linetype = 2;
	for (int i = 0; i < 2; i++)
	{
		coords[i] = 0;
		vertex[i] = 0;
		adjacent[i] = 0;
	}
	printed = false;
	deleted = false;
	color = 24;
	nbr_edge = number_edge++;
}


stl_vertex::stl_vertex(const stl_v & _coords)
{
	coords = _coords;
	edge_list = 0;
	deleted = false;
	to_remove = false;
	edge_to_remove = 0;
	precision = 2;
}

void stl_vertex::add_edge(stl_edge * edge)
{
	edge_list = new stl_edge_list(edge, edge_list);
}



stl_file::stl_file(const char * filename)
{
	is.open(filename, std::ios::in | std::ios::binary);
	ok = is.is_open();

	origin = stl_v_zero;
	transf = stl_l_id;
	name = filename;
	facet_list = 0;
	last_facet = 0;
	edge_list = 0;
	last_edge = 0;
	vertex_list = 0;
	prim_list = 0;
	last_prim = 0;
	surf_list = 0;
}

void stl_file::dump()
{
	cout << "Facets" << endl;
	stl_facet_list * cur = facet_list;
	while (cur)
	{
		cur->item->dump();
		cur = cur->next;
	}

	cout << "Edges" << endl;
	stl_edge_list * cure = edge_list;
	while (cure)
	{
		cure->item->dump();
		cure = cure->next;
	}
}

bool stl_file::token(const char * token)
{
	string s;
	if (is)
	{
		is >> s;
		if (s != token)
		{
			if (!silent)
				cout << "'" << s << "' found, '" << token << "' expected." << endl;
			return false;
		}
		return true;
	}
	return false;
}

bool stl_file::read_header()
{
	bool result = true;
	char header[80];
	memset(header, 0, 80);
	//is.setmode(ios_base::binary);
	is.read(header, 80);
	//	is >> header;

	int i = 0;
	bool cr = false;
	while (!cr && i < 80)
	{
		cr = (header[i] == '\r' || header[i] == '\n');
		i++;
	}
	if (cr)
	{
		stl_bin = false;
		is.seekg(0);
		//is.setmode(ios_base::text);
		result = token("solid");
		if (result)
		{
			is.get(header, 80);
			desc = header;
			desc.erase(0, desc.find_first_not_of(" \n\r\t"));
			desc.erase(desc.find_last_not_of(" \n\r\t")+1);
			//desc = desc. Trim();
		}
	}
	else
	{
		stl_bin = true;
		i = 79;
		while (i >= 0 && header[i] == ' ')
		{
			header[i] = 0;
			i--;
		}
		desc = string(header);
	}
	return(result);
}

bool stl_file::new_facet(const stl_v & normal, const stl_v vertex[])
{
	bool result = false;
	stl_facet * facet = new stl_facet;

	facet->read_normal = normal;

	bool identical = false;
	for (int i = 0; i < 3; i++)
	{
		//		facet->coords[i] = vertex[i] - origin;
		facet->coords[i] = transf * vertex[i];
		for (int j = 0; j < i; j++)
			//			if ((facet->coords[i] - facet->coords[j]).norm() < 0.001)
			if (facet->coords[i] == facet->coords[j])
				identical = true;
	}
	if (!identical)
	{
		facet->calculate_normal();
		//		facet->normal = ((facet->coords[1] - facet->coords[0]) * 
		//						(facet->coords[2] - facet->coords[0])).normaliz();
		if (facet->normal != facet->read_normal)
		{
			cout << "diff normal " << facet->normal << " ; " << facet->read_normal << " ; " << facet->normal - facet->read_normal << endl;
		}
		if (facet->normal.null())
			facet->normal = facet->read_normal;
		if (!facet->normal.null())
		{
			facet->nbr_facet = nbr_facet;
			stl_facet_list * facet_elem = new stl_facet_list;
			facet_elem->item = facet;
			facet_elem->next = 0;
			if (facet_list == 0)
				facet_list = facet_elem;
			if (last_facet != 0)
				last_facet->next = facet_elem;
			last_facet = facet_elem;
			result = true;
		}
		else
		{
			cout << "Error facet normal null !" << endl;
			cout << facet->coords[0] << " " << facet->coords[1] << " " << facet->coords[2] << " " << facet->read_normal << " " << facet->normal << endl;
			facet->dump();
		}
	}
	return(result);
}

int stl_file::read_facet_text()
{
	// Result :
	//	 0 : ok
	//	-1 : error 'facet' or 'endsolid' expected
	//	-2 : bad token found
	//	 1 : end of solid
	stl_v normal;
	stl_v vertex[3];

	string tok;
	is >> tok;
	if (tok == "endsolid")
		return 1;
	if (tok != "facet")
	{
		if (!silent)
			cout << "'" << tok << "' found, 'facet' or 'endsolid' expected." << endl;
		return -2;
	}
	bool result = token("normal");
	if (result) is >> normal;

	if (result) result = token("outer");
	if (result) result = token("loop");
	bool identical = false;
	for (int i = 0; i < 3 && result; i++)
	{
		if (result) result = token("vertex");
		if (result) is >> vertex[i];
	}
	if (result) result = token("endloop");
	if (result) result = token("endfacet");
	if (result) new_facet(normal, vertex);
	return(result ? 0 : -1);
}

stl_v read_v(ifstream & is)
{
	union {
		float f[3];
		char buf[12];
	} u;
	is.read(u.buf, 12);
	return(stl_v(u.f[0], u.f[1], u.f[2]));
}

bool stl_file::read_facet_bin()
{
	bool result;
	stl_v normal;
	stl_v vertex[3];

	normal = read_v(is);
	for (int i = 0; i < 3; i++)
		vertex[i] = read_v(is);

	char dummy[2];
	is.read(dummy, 2);
	result = !is.eof();
	if (result)
		result = new_facet(normal, vertex);
	return(result);
}

bool stl_file::read()
{
	bool result = true;
	nbr_facet = 0;
	if (stl_bin)
	{
		char header[80];
		is.read(header, 80);
		union {
			long l;
			char buf[4];
		} u;
		is.read(u.buf, 4);
		long face_number = u.l;
		while (nbr_facet < face_number && result)
		{
			if (result = read_facet_bin())
				nbr_facet++;
		}
		if (nbr_facet < face_number)
			cout << "face_number " << face_number << " face read " << nbr_facet << endl;
	}
	else
	{
		result = read_header();
		if (!result)
			return(result);
		int err = 0;
		bool done = false;
		while (is && !done)
		{
			err = read_facet_text();
			if (err == 0)
				nbr_facet++;
			else if (err < 0)
				cout << "error reading facet " << nbr_facet << endl;
			else if (err == 1)
				done = true;

		}
	}
	if (!silent && nbr_facet > 0)
	{
		cout << desc << endl;
		cout << nbr_facet << " facets readed" << endl;
	}
	return(result);
}


bool stl_facet::check()
{
	bool result = true;
	if (usage == 0)
		return(true);
	result = (usage == 3 || usage == 4);
	int max = usage > 4 ? 4 : usage;
	if (result && normal.null())
	{
		result = false;
		if (!silent)
			cout << "check facet " << nbr_facet << " normal = null" << endl;
	}
	for (int i = 0; i < max && result; i++)
	{
		if (adjacent[i] && adjacent[i]->usage == 0)
		{
			result = false;
			if (!silent)
				cout << "check facet " << nbr_facet << " adjacent[" << i << "] " << adjacent[i]->nbr_facet << " usage = 0" << endl;
		}
		if (!edge[i])
		{
			result = false;
			if (!silent)
				cout << "check facet " << nbr_facet << " edge[" << i << "] is null" << endl;
		}
		else if (edge[i]->deleted)
		{
			if (!silent)
				cout << "check facet " << nbr_facet << " edge[" << i << "] deleted" << endl;
		}

		if (!vertex[i])
		{
			if (!silent)
				cout << "check facet " << nbr_facet << " vertex[" << i << "] is null" << endl;
		}
		else if (vertex[i]->deleted)
		{
			if (!silent)
				cout << "check facet " << nbr_facet << " vertex[" << i << "] deleted" << endl;
		}
		else
		{
			if (coords[i] != vertex[i]->coords)
				if (!silent)
					cout << "check facet " << nbr_facet << " coords vertex[" << i << "]=" << vertex[i]->coords << " != " << coords[i] << endl;
		}
		if (result)
		{
			int a = edge[i]->side(this);
			if (edge[i]->adjacent[a] != this)
			{
				result = false;
				if (!silent)
					cout << "check facet " << nbr_facet << " no adjacent edge a=" << a << endl;
			}
			if (result && edge[i]->coords[a] != i)
			{
				result = false;
				if (!silent)
					cout << "check facet " << nbr_facet << " vertex " << edge[i]->coords[a] << " should be " << i << endl;
			}
		}

		for (int j = 0; j < max && result; j++)
		{
			if (j != i)
			{
				if ((coords[i] - coords[j]).norm() < optim_eps)
				{
					result = false;
					if (!silent)
						cout << "check facet " << nbr_facet << " same vertices " << i << " " << j << " "
						<< coords[i] << " " << coords[j] << " " << (coords[i] - coords[j]).norm() << endl;
				}
				else
				{
					stl_v vd = (coords[i] - coords[j]).normaliz();

					for (int k = 0; k < max && result; k++)
						if (k != j && k != i)
						{
							if (sapptdr(coords[k], coords[i], vd))
							{
								result = false;
								if (!silent)
									cout << "check facet " << nbr_facet << " colinear vertices " << i << " " << j << " " << k << endl;
								is_convex(!silent);
							}

						}
				}
			}
		}
	}
	return(result);
}

bool stl_edge::check()
{
	bool result = true;
	if (linetype == 0)
		return(true);
	result = (linetype == 2 || linetype == 5);
	if (result && adjacent[0] == 0)
	{
		result = false;
		if (!silent)
			cout << "check edge " << ident() << " : adjacent[0] == NULL" << endl;
	}
	if (result && adjacent[0] && adjacent[0] == adjacent[1])
	{
		result = false;
		if (!silent)
			cout << "check edge " << ident() << " : same adjacent facet" << endl;
	}
	if (result && adjacent[0] && adjacent[0]->usage == 0)
	{
		result = false;
		if (!silent)
			cout << "check edge " << ident() << " : adjacent facet[0] deleted" << endl;
	}
	if (result && adjacent[1] && adjacent[1]->usage == 0)
	{
		result = false;
		if (!silent)
			cout << "check edge " << ident() << " : adjacent facet[1] deleted" << endl;
	}

	if (result && adjacent[1] &&
		((coord(0) != coord(0, 1) && coord(0) != coord(1, 1)) ||
		(coord(1) != coord(1, 1) && coord(1) != coord(0, 1))))
	{
		result = false;
		if (!silent)
			cout << "check edge " << ident() << " : vertices do not correspond:"
			<< "[" << coords[0] << "]" << coord(0) << " != "
			<< "[" << coords[1] << "]" << coord(0, 1) << " & "
			<< "[" << (coords[0] + 1) % adjacent[0]->usage << "]" << coord(1) << " != "
			<< "[" << (coords[1] + 1) % adjacent[1]->usage << "]" << coord(1, 1)
			<< endl;
	}
	if (result && (coord(0) - coord(1)).norm() < optim_eps)
	{
		result = false;
		if (!silent)
			cout << "check edge " << ident() << " : identical vertices " << coord(0) << " " << coord(1) << " " << (coord(0) - coord(1)).norm() << endl;
	}
	if (result)
		for (int i = 0; i < 2; i++)
			if (adjacent[i] && adjacent[i]->adjacent[coords[i]] != adjacent[1 - i])
			{
				result = false;
				if (!silent)
				{
					cout << "check edge : adjacent[" << 1 - i << "] ";
					if (adjacent[1 - i])
						cout << adjacent[1 - i]->nbr_facet << " ";
					cout << " != facet ";
					if (adjacent[i] && adjacent[i]->adjacent[coords[i]])
						cout << adjacent[i]->adjacent[coords[i]]->nbr_facet;
					cout << endl;
				}
			}
	return(result);
}

bool stl_vertex::check()
{
	bool result = true;
	if (deleted)
		if (!silent)
			cout << "check vertex : deleted" << endl;
	stl_edge_list * cur = edge_list;
	while (cur)
	{
		if (cur->item->coord(0) != coords && cur->item->coord(1) != coords)
		{
			if (!silent)
				cout << "check vertex : coords do not correspond " << cur->item->ident() << endl;
			result = false;
		}
		cur = cur->next;
	}
	return(result);
}


bool stl_file::check()
{
	bool result = true;
	stl_facet_list * cur = facet_list;
	while (cur)
	{
		if (!cur->item->check())
			result = false;
		cur = cur->next;
	}
	stl_edge_list * cure = edge_list;
	while (cure)
	{
		if (!cure->item->deleted)
		{
			if (!cure->item->check())
				result = false;
			if (cure->item->adjacent[0] == 0)
			{
				if (!silent)
					cout << "check : " << cure->item->ident() << " edge adjacent[0] = NULL" << endl;
				result = false;
			}
			else if (cure->item->adjacent[1] == 0)
			{
				if (!silent)
					cout << "check : " << cure->item->ident() << " edge adjacent[1] = NULL" << endl;
				result = false;
			}
			else
			{
				stl_edge_list * curf = cure->next;
				while (curf)
				{
					if (!curf->item->deleted)
					{
						if (cure->item->linetype >= 0 && curf->item->linetype >= 0)
							if (curf->item->adjacent[0] == 0)
							{
								if (!silent)
									cout << "check : " << curf->item->ident() << " edge adjacent[0] = NULL" << endl;
								result = false;
							}
							else if (curf->item->adjacent[1] == 0)
							{
								if (!silent)
									cout << "check : " << curf->item->ident() << " edge adjacent[1] = NULL" << endl;
								result = false;
							}
							else if ((cure->item->coord(0) == curf->item->coord(0) && cure->item->coord(1) == curf->item->coord(1)) ||
								(cure->item->coord(0) == curf->item->coord(1) && cure->item->coord(1) == curf->item->coord(0)))
							{
								if (!silent)
									cout << "check : same edge " << cure->item->ident() << " and " << curf->item->ident() << endl;
								result = false;
							}
					}
					curf = curf->next;
				}
			}
		}
		cure = cure->next;
	}
	stl_vertex_list * curv = vertex_list;
	while (cur)
	{
		if (!curv->item->check())
			result = false;
		curv = curv->next;
	}
	return(result);
}



static int cur_rand_col = 0;

int random_color()
{
	cur_rand_col = (cur_rand_col + 1) % 16;
	return(cur_rand_col);
}

bool stl_facet::recalculate_precision()
{
	bool result = false;
	if (usage == 4 && !no_line4)
	{
		// On force la précision max pour tous les vertex
		for (int i = 0; i < 4; i++)
			precision = max(vertex[i]->precision, precision);
		for (int i = 0; i < 4; i++)
		{
			if (vertex[i]->precision < precision)
			{
				vertex[i]->precision = precision;
				result = true;
			}
		}
		int face = 0;
		if (!adjacent[face])
			face = 1 - face;
		if (adjacent[face])
		{
			double det = fabs(adjacent[face]->determinant());
			if (det > determinant_eps && precision < 8)
			{
				precision++;
				result = recalculate_precision();
			}
		}
	}
	return(result);
}

void stl_facet::write(ostream & os, int col3, int col4, int edgecolor, bool printface, bool printedge)
{
	if (!printed)
	{
		printed = true;
		switch (usage)
		{
			/*
					case 2:
						os << setprecision(5);
						if (printedge)
							os << "2 " << edgecolor << " " << coords[0] << " " << coords[1] << endl;
						break;
			*/
		case 3:
			if (!no_line3)
			{
				if (debug)
					os << "0 facet " << nbr_facet << " " << normal /*<< " " << read_normal*/ << " " << surface3pt(coords[0], coords[1], coords[2]) << endl;
				os << setprecision(7);
				if (stl_raw)
				{
					if (printface)
						os << "3 " << col3 << " " << coords[0] << " " << coords[1] << " " << coords[2] << endl;
					if (!printedge)
					{
						os << "2 " << edgecolor << " " << coords[0] << " " << coords[1] << endl;
						os << "2 " << edgecolor << " " << coords[1] << " " << coords[2] << endl;
						os << "2 " << edgecolor << " " << coords[2] << " " << coords[0] << endl;
					}
				}
				else
				{
					if (printface)
						os << "3 " << col3 << " " << vertex[0] << " " << vertex[1] << " " << vertex[2] << endl;
					if (!printedge)
					{
						os << "2 " << edgecolor << " " << vertex[0] << " " << vertex[1] << endl;
						os << "2 " << edgecolor << " " << vertex[1] << " " << vertex[2] << endl;
						os << "2 " << edgecolor << " " << vertex[2] << " " << vertex[0] << endl;
					}
				}
			}
			break;
		case 4:
			if (!no_line4)
			{
				if (debug)
					os << "0 facet " << nbr_facet << " " << normal /*<< " " << read_normal*/ << endl;
				os << setprecision(7);
				cout << setprecision(7);
				if (printface)
				{
					//					int tmpcol = col4;
					//					double d = determinant(-1);
					//if (fabs(d) > determinant_eps)
					//					if (d != 0)
					//					{
					//	tmpcol = 4;
					//cout << "4 " << 4 << " " << vertex[0] << " " << vertex[1] << " " << vertex[2] << " " << vertex[3] << " det=" << d << endl;
					//					}
					os << "4 " << col4 << " " << vertex[0] << " " << vertex[1] << " " << vertex[2] << " " << vertex[3] << endl;
					//					for (int i=0;i<4;vertex[i++]->precision=4);os << "4 " << col4 << " " << vertex[0] << " " << vertex[1] << " " << vertex[2] << " " << vertex[3] << endl << "0 // (3)" << determinant(3) << " (4)" << determinant(4) << " (5)" << determinant(5) << endl;
				}
				if (!printedge)
				{
					os << "2 " << edgecolor << " " << vertex[0] << " " << vertex[1] << endl;
					os << "2 " << edgecolor << " " << vertex[1] << " " << vertex[2] << endl;
					os << "2 " << edgecolor << " " << vertex[2] << " " << vertex[3] << endl;
					os << "2 " << edgecolor << " " << vertex[3] << " " << vertex[0] << endl;
				}
			}
			break;
			/*
					case 5:
						os << setprecision(5);
						if (printedge)
							os << "5 " << edgecolor << " " << coords[0] << " " << coords[1] << " " << coords[2] << " " << coords[3] << endl;
						break;
			*/
		}
		for (int i = 0; i < usage; i++)
			if (edge[i])
				edge[i]->write(os, col_line2, col_line5);
	}
}

string stl_facet::ident()
{
	char buf[16];
	sprintf(buf, "%d", nbr_facet);
	return(string(buf));
}

void stl_facet::calculate_normal()
{
	normal = ((coords[1] - coords[0]) *
		(coords[2] - coords[0])).normaliz();
}

void stl_facet::dump()
{
	cout << ident() << " : ";
	if (usage == 0)
	{
		cout << "deleted" << endl;
		return;
	}
	for (int i = 0; i < usage; i++)
		if (adjacent[i])
			cout << adjacent[i]->ident() << " ";
		else
			cout << ". ";
	cout << "  \t[ ";
	for (int i = 0; i < usage; i++){
            if (edge[i])
		cout << edge[i]->nbr_edge << " ";
            else
                cout << "NULL ";
        }
	cout << "]\t";
	for (int i = 0; i < usage; i++)
		cout << "<" << coords[i] << "> ";
	cout << endl;
}

void stl_edge::dump()
{
	cout << ident() << " : ";
	if (deleted)
	{
		cout << "deleted" << endl;
		return;
	}
	for (int i = 0; i < 2; i++)
		cout << coords[i] << " ";
	cout << "  ";
	for (int i = 0; i < 2; i++)
		cout << "<" << vertex[i]->coords << ">  ";
	cout << endl;
}

void stl_vertex::dump()
{
	cout << coords << " : ";
	if (deleted)
	{
		cout << "deleted" << endl;
		return;
	}
	stl_edge_list * cur = edge_list;
	while (cur)
	{
		cout << cur->item->ident() << " ";
		cur = cur->next;
	}
	cout << endl;
}

double small0(double v);

ostream & operator << (
	ostream & s,
	const	stl_vertex * v
	)
{
	if (v)
	{
		return (s <<
			small0(ftrunc(v->coords.x_coord(), -v->precision)) << " " <<
			small0(ftrunc(v->coords.y_coord(), -v->precision)) << " " <<
			small0(ftrunc(v->coords.z_coord(), -v->precision)));
	}
	else
	{
		return (s << "0 0 0");
	}
}

void stl_edge::write(ostream & os, int col2, int col5) {
	if (!printed)
	{
		switch (linetype)
		{
		case 2:
			if (!no_line2)
			{
				os << setprecision(7);
				//			os << "2 " << col2 << " " << coord(0) << " " << coord(1) << endl;
				os << "2 " << color << " " << vertx(0) << " " << vertx(1) << endl;
			}
			break;
		case 5:
			if (!no_line5)
			{
				os << setprecision(7);
				//			os << "5 " << col5 << " " << coord(0) << " " << coord(1) << " " << coord(2) << " " << coord(2,1) << endl;
				os << "5 " << color << " " << vertx(0) << " " << vertx(1) << " " << vertx(2) << " " << vertx(2, 1) << endl;
			}
			break;
		}
		printed = true;
	}
}

string stl_edge::ident()
{
	char buf[16];
	sprintf(buf, "%d", nbr_edge);
	string a0 = "-";
	string a1 = "-";
	if (adjacent[0])
		a0 = adjacent[0]->ident();
	if (adjacent[1])
		a1 = adjacent[1]->ident();
	return(string("[") + string(buf) + ":" + a0 + "," + a1 + "]");
}


void stl_prim::write(ostream & os)
{
	if (!no_line1)
	{
		if (bfc && !cw)
			os << "0 BFC INVERTNEXT" << endl;
		os << setprecision(5);
		os << "1 " << color << " " << transfo << " " << name << ".dat" << endl;
	}
}

std::string str_tolower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), 
                // static_cast<int(*)(int)>(std::tolower)         // wrong
                // [](int c){ return std::tolower(c); }           // wrong
                // [](char c){ return std::tolower(c); }          // wrong
                   [](unsigned char c){ return std::tolower(c); } // correct
                  );
    return s;
}

void stl_file::write_dat(bool print_geom)
{
	eps_empile(0.0001);
	string dat_name;
	if (ldr_out)
		dat_name = out_filename;
	else
	{
		dat_name = name;
		string extension;
		if (ldr_opt)
			extension = ".ldr";
		else
			extension = ".dat";
		dat_name = str_tolower(dat_name);
		size_t pos = dat_name.find(".stl");
		if (pos != std::string::npos)
			dat_name = dat_name.replace(pos, 4, extension);
		else
			dat_name += extension;
	}
	ofstream fout(dat_name);
	fout << "0 " << partname << endl;
	int f = max(dat_name.rfind("/"), dat_name.rfind("\\"));
	if (f > 0)
		dat_name = dat_name.substr(f + 1);
	fout << "0 Name: " << dat_name << endl;
	fout << "0 Author: " << author << endl;
	if (ldraw)
	{
		if (name.substr(0, 2) == "s\\")
			fout << "0 !LDRAW_ORG Unofficial_Subpart" << endl;
		else if (name.substr(0, 2) == "p\\")
			fout << "0 !LDRAW_ORG Unofficial_Primitive" << endl;
		else if (name.substr(0, 5) == "p\\48\\")
			fout << "0 !LDRAW_ORG Unofficial_48_Primitive" << endl;
		else
			fout << "0 !LDRAW_ORG Unofficial_Part" << endl;
		fout << "0 !LICENSE Redistributable under CCAL version 2.0 : see CAreadme.txt" << endl;
	}
	fout << endl;
	if (bfc)
		fout << "0 BFC CERTIFY CCW" << endl;
	fout << endl;
	fout << "0 // Created with stl2dat conversion tool v" << version << endl;
	stl_prim_list * curp = prim_list;
	int prim_cpt = 0;
	while (curp)
	{
		curp->item->write(fout);
		curp = curp->next;
		prim_cpt++;
	}
	if (!silent)
		cout << prim_cpt << " primitives detected" << endl;

	bool precision_ok = false;
	while (!precision_ok)
	{
		precision_ok = true;
		stl_facet_list * curf = facet_list;
		while (curf)
		{
			precision_ok = precision_ok && !curf->item->recalculate_precision();
			curf = curf->next;
		}
	}

	int surf_cpt = 0;
	stl_surf_list * curs = surf_list;
	while (curs)
	{
		stl_facet_list * cur = facet_list;
		stl_surf * surf = curs->item;
		if (!surf->prim)
		{
			fout << endl;
			//		eps_empile( 0.01);
			if (print_geom)
				curs->item->write(fout);
			//		eps_depile();
			int col3 = col_line3;
			int col4 = col_line4;
			if (rand_col_surfaces || (rand_col_adj_cyl && surf->adj_cyl))
				col3 = col4 = random_color();
			while (cur)
			{
				if (cur->item->surf == curs->item)
					cur->item->write(fout, col3, col4);
				cur = cur->next;
			}
		}
		curs = curs->next;
		surf_cpt++;
	}
	if (!silent)
		cout << surf_cpt << " geometric surfaces detected" << endl;

	fout << endl;
	stl_facet_list * cur = facet_list;
	while (cur)
	{
		int col3 = col_line3;
		int col4 = col_line4;
		if (rand_col_surfaces)
			col3 = col4 = random_color();
		if (cur->item->color != 16)
			col3 = col4 = cur->item->color;
		cur->item->write(fout, col3, col4);
		cur = cur->next;
	}
	stl_edge_list * cure = edge_list;
	while (cure)
	{
		if (cure->item->linetype > 0)
			cure->item->write(fout);
		cure = cure->next;
	}

	eps_depile();
}

void stl_file::add_edge(int usage, const stl_v& v0, const stl_v& v1, const stl_v& v2, const stl_v& v3)
{
	/*!!
	stl_facet edge;
	edge.usage = usage;
	edge.coords[0] = v0;
	edge.coords[1] = v1;
	edge.coords[2] = v2;
	edge.coords[3] = v3;
	stl_facet_list * edge_elem = new stl_facet_list;
	edge_elem->item = edge;
	edge_elem->next = 0;
	if (edge_list == 0)
		edge_list = edge_elem;
	if (last_edge != 0)
		last_edge->next = edge_elem;
	last_edge = edge_elem;
	*/
}

// create vertex if necessary, link to facet and edge
void create_vertex(stl_vertex_list * & vertex_list, stl_edge * edge, stl_facet * fct1, stl_facet * fct2, int i1, int i2)
{
	bool reaffect = false;
	stl_vertex * result = fct1->vertex[i1];
	stl_vertex * vertex2 = 0;
	if (!result && fct2)
		result = fct2->vertex[i2];
	else
		if (fct2 && fct2->vertex[i2] && !fct2->vertex[i2]->deleted && fct1->vertex[i1] != fct2->vertex[i2])
		{
			// merge 2 edge lists
			vertex2 = fct2->vertex[i2];
			stl_edge_list * prev = result->edge_list;
			while (prev->next)
				prev = prev->next;
			prev->next = vertex2->edge_list;
			vertex2->deleted = true;
			reaffect = true;
		}
	if (!result)
	{
		result = new stl_vertex(fct1->coords[i1]);
		vertex_list = new stl_vertex_list(result, vertex_list);
	}
	result->add_edge(edge);
	if (reaffect)
	{
		stl_edge_list * cur = vertex2->edge_list;
		while (cur)
		{
			for (int i = 0; i < 2; i++)
			{
				if (cur->item->vertex[i] == vertex2)
					cur->item->vertex[i] = result;
				stl_facet * adj = cur->item->adjacent[i];
				if (adj)
				{
					for (int v = 0; v < adj->usage; v++)
						if (adj->vertex[v] == vertex2)
							adj->vertex[v] = result;
				}
			}
			cur = cur->next;
		}
	}
	fct1->vertex[i1] = result;
	if (fct2)
		fct2->vertex[i2] = result;
	edge->vertex[edge->coord(0) == result->coords ? 0 : 1] = result;
}

void add_edge2(
	stl_edge_list * & edge_list,
	stl_vertex_list * & vertex_list,
	stl_facet * fct1,
	stl_facet * fct2,
	int i1, int i2, int j2
)
{
	//cout << "add_edge2( [" << fct1->ident() << "], [" << fct2->ident() << "], " << i1 << ", " << i2 << ", " << j2 << ")" << endl;
	stl_edge * edge = new stl_edge();
	fct1->edge[i1] = edge;
	fct1->adjacent[i1] = fct2;
	edge->adjacent[0] = fct1;
	edge->adjacent[1] = fct2;
	edge->coords[0] = i1;	// number of vertex in adjacent[0] (minimum)
	edge->coords[1] = i2;	// number of vertex in adjacent[1] (minimum)

	int k1 = (i1 + 1) % fct1->usage;
	int k2 = 0;
	if (fct2)
		k2 = (i2 == j2 ? (i2 + 1) % fct2->usage : i2);
	create_vertex(vertex_list, edge, fct1, fct2, i1, j2);
	create_vertex(vertex_list, edge, fct1, fct2, k1, k2);
	if (fct2)
	{
		fct2->edge[i2] = edge;
		fct2->adjacent[i2] = fct1;
	}
	edge_list = new stl_edge_list(edge, edge_list);
}

// Calculate edges between 2 facets
int stl_file::create_edge(stl_facet * fct1, stl_facet * fct2)
{
	int result = 0;
	for (int i1 = 0; i1 < fct1->usage; i1++)
		if (!fct1->edge[i1])
			for (int i2 = 0; i2 < fct2->usage; i2++)
				if (!fct2->edge[i2])
				{
					int j1 = (i1 + 1) % fct1->usage;
					int j2 = (i2 + 1) % fct2->usage;
					if (fct1->coords[i1] == fct2->coords[i2] && fct1->coords[j1] == fct2->coords[j2])
					{
						result++;
						add_edge2(edge_list, vertex_list, fct1, fct2, i1, i2, i2);
					}
					else if (fct1->coords[i1] == fct2->coords[j2] && fct1->coords[j1] == fct2->coords[i2])
					{
						result++;
						add_edge2(edge_list, vertex_list, fct1, fct2, i1, i2, j2);
					}
				}
	return(result);
}

// Create links between topological elements : vertex, edges, facets
void stl_file::topo()
{
	int nb_edge = 0;
	stl_facet_list * cur = facet_list;
	stl_facet * fct1;
	while (cur)
	{
		fct1 = cur->item;
		if (fct1->usage)
		{
			stl_facet_list * cur2 = cur->next;
			while (cur2)
			{
				stl_facet * fct2 = cur2->item;
				if (fct2->usage)
					nb_edge += create_edge(fct1, fct2);
				cur2 = cur2->next;
			}
			for (int i1 = 0; i1 < fct1->usage; i1++)
				if (!fct1->edge[i1])
				{
					add_edge2(edge_list, vertex_list, fct1, NULL, i1, -1, -1);
				}
			nb_edge++;
		}
		cur = cur->next;
	}
	if (!silent)
		cout << nb_edge << " edges created" << endl;
#ifdef DEBUG
	/*
	if (!nodebug)
		{
			stl_edge_list * cure = edge_list;
			while (cure)
			{
				if (cure->item->adjacent[1])
					cout << "edge adj(" << cure->item->adjacent[0]->nbr_facet << ","
						 << cure->item->adjacent[1]->nbr_facet << ")" << endl;
				else
					cout << "edge end(" << cure->item->adjacent[0]->nbr_facet << ")" << endl;
				cure = cure->next;
			}
			stl_facet_list * curf = facet_list;
			while (curf)
			{
				cout << "facet(" << curf->item->nbr_facet << ") " << endl;
				for (int i = 0 ; i < curf->item->usage ; i++)
					if (curf->item->adjacent[i])
						cout << "adj[" << i << "]=(" << curf->item->adjacent[i]->nbr_facet << ")" << endl;
				curf = curf->next;
			}
		}
	*/
#endif
	//check();
}

int other(int i, int j)
{
	int result = 0;
	if (i != 1 && j != 1)
		result = 1;
	if (i != 2 && j != 2)
		result = 2;
	return(result);
}

// Merge 2 triangles to a triangle
bool merge3(stl_facet * fct1, stl_facet * fct2)
{
	bool result = false;
	if (fct1->normal == fct2->normal && fct1->usage == 3 && fct2->usage == 3)
	{
		bool result = false;
		for (int i1 = 0; i1 < 3 && !result; i1++)
			for (int i2 = 0; i2 < 3 && !result; i2++)
				if (fct1->coords[i1] == fct2->coords[i2])
				{
					for (int j1 = 0; j1 < 3 && !result; j1++)
						for (int j2 = 0; j2 < 3 && !result; j2++)
							if (j1 != i1 && j2 != i2 && fct1->coords[j1] == fct2->coords[j2])
							{
								int k1 = other(i1, j1);
								int k2 = other(i2, j2);
								if (fct1->coords[i1].colinear(fct1->coords[k1], fct2->coords[k2]))
								{
									fct1->coords[i1] = fct2->coords[k2];
									fct2->usage = 0;
									return(true);
								}
								else if (fct1->coords[j1].colinear(fct1->coords[k1], fct2->coords[k2]))
								{
									fct1->coords[j1] = fct2->coords[k2];
									fct2->usage = 0;
									return(true);
								}
							}

				}
	}
	return(result);
}

// Merge 2 consecutives triangles to a triangle
void stl_file::optim3()
{
	int opt3 = 0;
	stl_facet_list * cur = facet_list;
	stl_facet * fct1;
	stl_facet * fct2;
	fct1 = cur->item;
	cur = cur->next;
	while (cur)
	{
		fct2 = cur->item;
		if (merge3(fct1, fct2))
			opt3++;
		stl_facet_list * cur2 = cur->next;
		for (int c = 0; c < 200 && cur2; c++)
		{
			if (merge3(fct1, cur2->item))
				opt3++;
			else
				cur2 = cur2->next;
		}
		cur = cur->next;
		fct1 = fct2;
	}
	if (!silent)
		cout << opt3 << " triangle facets merged" << endl;
}

void stl_facet::calcul_angles(double & totag, double & minag, double & maxag)
{
	totag = 0;
	minag = 1000;
	maxag = -1000;
	for (int a = 0; a < 4; a++)
	{
		int b = (a + 1) % 4;
		int c = (a + 2) % 4;
		double ag = (coords[a] - coords[b]).angle(coords[c] - coords[b], normal);
		totag += ag;
		if (ag < minag)
			minag = ag;
		if (ag > maxag)
			maxag = ag;
	}
}

bool stl_facet::is_convex(bool trace)
{
	if (normal.null())
		return(false);
	bool result = true;
	double totag;
	double minag;
	double maxag;
	calcul_angles(totag, minag, maxag);
	if (trace)
		cout << "(" << nbr_facet << ") tot=" << degre(totag) << "° min=" << degre(minag) << "° max=" << degre(maxag) << "°" << endl;
	if (equal(totag, 4 * pi))
	{
		stl_v tmp = coords[0];
		coords[0] = coords[1];
		coords[1] = coords[2];
		coords[2] = tmp;
		calcul_angles(totag, minag, maxag);
		if (trace)
			cout << "(" << nbr_facet << ") tot=" << degre(totag) << "° min=" << degre(minag) << "° max=" << degre(maxag) << "°" << endl;
	}
	result = (equal(totag, 2 * pi) && maxag < pi - ag_lim_q) || (equal(totag, 6 * pi) && minag > pi + ag_lim_q);
	return(result);
}

double stl_facet::surface()
{
	if (m_surface < 0)
	{
		m_surface = 0;
		if (usage >= 3)
		{
			m_surface = surface3pt(coords[0], coords[1], coords[2]);
			if (usage == 4)
				m_surface += surface3pt(coords[1], coords[2], coords[3]);
		}
	}
	return(m_surface);
}

double stl_facet::surface(const stl_v & v_ori, const stl_v & v_repl)
{
	double result = 0;
	if (usage >= 3)
	{
		stl_v c[4];
		for (int i = 0; i < usage; i++)
			c[i] = coords[i] == v_ori ? v_repl : coords[i];

		result = surface3pt(c[0], c[1], c[2]);
		if (usage == 4)
			result += surface3pt(c[1], c[2], c[3]);
	}
	return(result);
}

double stl_facet::determinant(int precision)
{
	double result = 0;
	if (usage >= 4)
	{
		stl_v v[4];
		for (int i = 0; i < 4; i++)
		{
			if (precision == -1)
				v[i] = vertex[i]->coords.ftrunc(-vertex[i]->precision);
			else
				v[i] = vertex[i]->coords.ftrunc(-max(precision, -vertex[i]->precision));
		}
		result = determinant4(v);
	}
	return(result);
}


// Merge 2 triangles to a rectangle if possible
bool merge4(stl_facet * fct1, stl_facet * fct2)
{
	return false;
}

int stl_edge::side(stl_facet * m_facet)
{
	if (adjacent[0] == m_facet)
		return(0);
	else
		return(1);
}

// Reaffect pointer to new facet
void stl_edge::reaffect(stl_facet * m_facet, stl_facet * repl_facet, int ve)
{
	int s = side(m_facet);
	adjacent[s] = repl_facet;
	coords[s] = ve;
	if (coords[1 - s] >= 0 && coords[1 - s] <= 3 && adjacent[1 - s])
		adjacent[1 - s]->adjacent[coords[1 - s]] = repl_facet;
}

// Reaffect border edges to new merged facet
void stl_facet::merge_facet(stl_edge * m_edge, stl_facet * m_facet, int ve0, int ve1)
{
	for (int i = 3; i > ve0 + 1; i--)
	{
		edge[i] = edge[i - 1];
		adjacent[i] = adjacent[i - 1];
		edge[i]->coords[edge[i]->side(this)] = i;
	}
	edge[ve0] = m_facet->edge[(ve1 + 1) % 3];
	edge[ve0 + 1] = m_facet->edge[(ve1 + 2) % 3];
	adjacent[ve0] = m_facet->adjacent[(ve1 + 1) % 3];
	adjacent[ve0 + 1] = m_facet->adjacent[(ve1 + 2) % 3];
	edge[ve0]->reaffect(m_facet, this, ve0);
	edge[ve0 + 1]->reaffect(m_facet, this, ve0 + 1);
}

// Merge 2 adjacent triangles facet to 1 rectangle facet
bool stl_edge::merge_adj()
{
	//	cout << this << " merge_adj" << endl;
	bool result = false;
	if (adjacent[0] && adjacent[1] &&
		(adjacent[0]->normal == adjacent[1]->normal || adjacent[0]->normal == -adjacent[1]->normal) &&
		adjacent[0]->usage == 3 && adjacent[1]->usage == 3)
	{
		// Insert new point to rectangle facet
		int j = coords[0] + 1;
		for (int i = 3; i > j; i--)
		{
			adjacent[0]->coords[i] = adjacent[0]->coords[i - 1];
			adjacent[0]->vertex[i] = adjacent[0]->vertex[i - 1];
		}
		adjacent[0]->coords[j] = coord(2, 1);
		adjacent[0]->vertex[j] = vertx(2, 1);
		// verify convex quadrilater
		adjacent[0]->usage = 4;
		int p = 2;
		double det = fabs(adjacent[0]->determinant(p));
		if (det > determinant_eps)
		{
			p++;
			while (p < 8 && det > determinant_eps)
			{
				det = fabs(adjacent[0]->determinant(p));
				if (det > determinant_eps)
					p++;
			}
			adjacent[0]->precision = max(adjacent[0]->precision, p);
			for (int i = 0; i < 4; i++)
			{
				adjacent[0]->vertex[i]->precision = max(adjacent[0]->vertex[i]->precision, p);
			}
			//if (det > determinant_eps)
			//{
			//cout << setprecision(7);
			//cout << "4 " << 1 << " " << adjacent[0]->vertex[0] << " " << adjacent[0]->vertex[1] << " " << adjacent[0]->vertex[2] << " " << adjacent[0]->vertex[3] << " p=" << p << " det=" << det << endl;
			//}
						//cout << "merge4 " << det << endl;
						//adjacent[0]->color = 4;
		}
		if (adjacent[0]->is_convex() && det <= determinant_eps)
		{
#ifdef DEBUG
			if (!nodebug)
				cout << "merge4 " << adjacent[0]->nbr_facet << " " << adjacent[1]->nbr_facet << endl;
#endif
			result = true;
			adjacent[0]->merge_facet(this, adjacent[1], coords[0], coords[1]);
			adjacent[0]->usage = 4;
			adjacent[1]->usage = 0;
			linetype = -1;
			//cout << "merge " << adjacent[0]->ident() << " with " << adjacent[1]->ident() << endl;

			//adjacent[0]->write(fout, 2);
		}
		else
		{
			//adjacent[0]->usage = 4;
			//adjacent[0]->write(fout, 1);
			for (int i = j; i < 3; i++)
			{
				adjacent[0]->coords[i] = adjacent[0]->coords[i + 1];
				adjacent[0]->vertex[i] = adjacent[0]->vertex[i + 1];
			}
			adjacent[0]->usage = 3;
		}
	}
	return(result);
}

// next edge
stl_edge * stl_edge::next(int & vx, stl_facet * & s_fct)
{
	/*
	if (!adjacent[1])
		return(0);
	int n_side = adjacent[0] == s_fct ? 1 : 0;
	stl_facet * n_fct = adjacent[n_side];
	stl_v s_v = coord( vx, 1 - n_side);
	bool ok = false;
	for (int vi = 0; !ok && vi < n_fct->usage; vi++)
	{
		stl_v c = coord( vi, n_side);
		if (s_v == c)
		{
			vx = vi;
			ok = true;
		}
	}
	if (!ok)
		return(0);

	s_fct = n_fct;
	if (n_fct)
	{
		if ( n_fct->edge[vx] != this)
			return( n_fct->edge[vx]);
		else
			return( n_fct->edge[(vx+1) % n_fct->usage]);
	}
	else
	*/
	return(0);
}

void stl_facet::replace_vertex(stl_vertex * ve_ori, stl_vertex * ve_repl)
{
	for (int i = 0; i < usage; i++)
		if (vertex[i] == ve_ori)
		{
			vertex[i] = ve_repl;
			coords[i] = ve_repl->coords;
			int j = edge[i]->side(this);
			edge[i]->vertex[j] = ve_repl;
			/*
						for (int j = 0; j < 2; j++)
							if (edge[i]->vertex[j] == ve_ori)
							{
								edge[i]->coords[j] = i;
							}
			*/
		}
	//		ed->vertex[i] = ve_repl;
}

void stl_facet::replace_adjacent(stl_vertex * ve_ori, stl_edge * ed, bool del_edge)
{
	if (usage == 0)
		return;
	if (del_edge)
		for (int i = 0; i < usage; i++)
			if (edge[i]->vertex[0] == ve_ori || edge[i]->vertex[1] == ve_ori)
				edge[i]->deleted = true;
	int s = side_edge(ed);
	stl_edge * ed0 = edge[(s + 1) % usage];
	stl_edge * ed1 = edge[(s + usage - 1) % usage];
	stl_facet * fct0 = adjacent[(s + 1) % usage];
	stl_facet * fct1 = adjacent[(s + usage - 1) % usage];
	int ad0;

	if (fct0)
	{
		// replace adjencent pointers for facets
		ad0 = fct0->side_edge(ed0);
		fct0->adjacent[ad0] = fct1;

		for (int i = 0; i < fct0->usage; i++)
			if (fct0->edge[i] == ed1)
				fct0->edge[i] = ed0;
		fct0->calculate_normal();
	}

	// replace adjacent pointers for edges
	ad0 = ed0->adjacent[0] == this ? 0 : 1;
	ed0->adjacent[ad0] = fct1;

	if (fct1)
	{
		// replace adjencent pointers for facets
		int ad1 = fct1->side_edge(ed1);
		fct1->adjacent[ad1] = fct0;

		// replace coords in modified edges
		for (int i = 0; i < fct1->usage; i++)
		{
			if (ed0->vertex[ad0] == fct1->vertex[i])
				ed0->coords[ad0] = i;
			if (fct1->edge[i] == ed1)
				fct1->edge[i] = ed0;
		}
		fct1->calculate_normal();
	}

	ed1->deleted = true;
	//	if (ed1->vertex[ad1] == fct0->vertex[i])
	//		ed1->coords[ad1] = i;

	usage = 0;
}

/*
	// replace adjency pointers for facets
	int ad0 = fct0->side_edge( ed0);
	int ad1 = fct1->side_edge( ed1);
	fct0->adjacent[ad0] = fct1;
	fct1->adjacent[ad1] = fct0;

	// replace adjacent pointers for edges
	ad0 = ed0->adjacent[0] == this ? 0 : 1;
	//ad1 = ed1->adjacent[0] == this ? 0 : 1;
	ed0->adjacent[ad0] = fct1;
	//ed1->adjacent[ad1] = fct0;

	// replace coords in modified edges
	for (int i = 0; i < fct1->usage; i++)
	{
		if (ed0->vertex[ad0] == fct1->vertex[i])
			ed0->coords[ad0] = i;
		if (fct1->edge[i] == ed1)
			fct1->edge[i] = ed0;
	}
	for (i = 0; i < fct0->usage; i++)
		if (fct0->edge[i] == ed1)
			fct0->edge[i] = ed0;
	ed1->deleted = true;
	//	if (ed1->vertex[ad1] == fct0->vertex[i])
	//		ed1->coords[ad1] = i;
	fct0->calculate_normal();
	fct1->calculate_normal();
*/

// remove a vertex
void stl_vertex::remove_vertex(stl_edge * ed_supr)
{
	stl_facet * fct0 = ed_supr->adjacent[0];
	stl_facet * fct1 = ed_supr->adjacent[1];
	stl_vertex * ve_repl = ed_supr->vertex[ed_supr->vertex[0] == this ? 1 : 0];
	if (!ve_repl)
		return;
	//cout << "remove_vertex(" << coords << ") repl=" << ve_repl->coords << endl;

	stl_edge_list * cur = edge_list;
	while (cur)
	{
		stl_edge * ed = cur->item;
		if (ed != ed_supr)
		{
			stl_facet * fct = ed->adjacent[0];
			if (fct/* && fct != fct0 && fct != fct1*/)
				fct->replace_vertex(this, ve_repl);

			fct = ed->adjacent[1];
			if (fct/* && fct != fct0 && fct != fct1*/)
				fct->replace_vertex(this, ve_repl);
		}
		cur = cur->next;
	}
	if (fct0)
		fct0->replace_adjacent(this, ed_supr, false);
	if (fct1)
		fct1->replace_adjacent(this, ed_supr, false);

	deleted = true;
	ed_supr->deleted = true;
}

// scan all edges around a vertex, search edge to remove : compare angles
/*
void stl_vertex::scan_edges()
{
	cout << "scan_edges(" << coords << ")" << endl;
	stl_edge_list * cur = edge_list;
	ed_min = 0;
	stl_edge * ed_found = 0;
	double ag_min = 99999;
	while (cur)
	{
		stl_edge * ed = cur->item;
		if (ed->adjacent[1])
		{
			double ag = ed->adjacent[0]->angle_between( ed->adjacent[1]);
			if (ag < ag_min)
			{
				ag_min = ag;
				ed_found = ed;
			}
			cout << "  fct(" << ed->ident() << ") ag=" << degre(ag) << "°" << " ag_min=" << degre(ag_min) << "°" << endl;
		}
		cur = cur->next;
	}
	if (ed_found && ag_min < ag_lim_t && ag_min > 0.01)
	{
		ed_min = ed_found;
		ed_min->color = 4;
		to_remove = true;
	}
}
*/

/*
// scan all edges around a vertex, search edge to remove : check edge length
void stl_vertex::scan_edges()
{
	//cout << "scan_edges(" << coords << ")" << endl;
	stl_edge_list * cur = edge_list;
	stl_edge * ed_min = 0;
	stl_edge * ed_found = 0;
	double ln_min = 99999;
	while (cur)
	{
		stl_edge * ed = cur->item;
		double ln = (ed->coord(0) - ed->coord(1)).norm();
		if (ln < ln_min)
		{
			ln_min = ln;
			ed_found = ed;
		}
		//cout << "  fct(" << ed->ident() << ") ln=" << ln << "°" << " ln_min=" << ln_min << "°" << endl;
		cur = cur->next;
	}
//	if (ed_found && ln_min < 1 && ln_min > 0.01)
	if (ed_found && ln_min < 0.005)
	{
		ed_min = ed_found;
		ed_min->color = 4;
		edge_to_remove = ed_min;
		to_remove = true;
	}
}
*/
// scan all edges around a vertex, search edge to remove : check edge length
/*
void stl_vertex::scan_edges()
{
	//cout << "scan_edges(" << coords << ")" << endl;
	stl_edge_list * cur = edge_list;
	stl_edge * ed_found = 0;
	while (cur && !ed_found)
	{
		stl_edge * ed = cur->item;
int nbr_f;
if (ed->adjacent[0]) nbr_f = ed->adjacent[0]->nbr_facet;
if (ed->adjacent[1]) nbr_f = ed->adjacent[1]->nbr_facet;
		if (ed->coord(0) == ed->coord(1))
		{
			ed_found = ed;
		}
		//cout << "  fct(" << ed->ident() << ") ln=" << ln << "°" << " ln_min=" << ln_min << "°" << endl;
		cur = cur->next;
	}
	if (ed_found)
	{
		edge_to_remove = ed_found;
		to_remove = true;
	}
}
*/
void stl_vertex::scan_edges(int nb_face)
{
	//cout << "scan_edges(" << coords << ")" << endl;
	stl_v vn[64];
	int n = 0;
	stl_edge_list * cur = edge_list;
	while (cur)
	{
		stl_edge * ed = cur->item;
		for (int a = 0; a < 2 && ed->adjacent[a]; a++)
		{
			bool found = false;
			for (int i = 0; i < n; i++)
				if (vn[i] == ed->adjacent[a]->normal)
					found = true;
			if (!found && n < 64)
				vn[n++] = ed->adjacent[a]->normal;
		}

		cur = cur->next;
	}
	if (n <= 2)
	{
		//cout << coords << " n=" << n << endl;
		cur = edge_list;
		while (cur)
		{
			stl_edge * ed = cur->item;
			//ed->color = n;
			//for (int a = 0; a < 2 && ed->adjacent[a]; a++)
			//	if (ed->adjacent[a]->color > n)
			//		ed->adjacent[a]->color = n;
			cur = cur->next;
		}
	}
	if (n == nb_face)
	{
		// If the surface, find an edge that can be removed 
		// without violating the boudaries of the surface.
		// This is done by checking the surface of the surface
		// The surface must be kept the same after a vertex has been removed.

		double init_surface = 0;
		cur = edge_list;
		while (cur)
		{
			for (int a = 0; a < 2 && cur->item->adjacent[a]; a++)
				init_surface += cur->item->adjacent[a]->surface();
			cur = cur->next;
		}
		init_surface /= 2;
		//cout << "init_surface=" << init_surface << endl;
		cur = edge_list;
		bool found = false;
		while (cur && !found)
		{
			double cur_surface = 0;
			stl_v repl = (coords == cur->item->coord(1) ? cur->item->coord(0) : cur->item->coord(1));
			stl_edge_list * cur2 = edge_list;
			while (cur2)
			{
				for (int a = 0; a < 2 && cur2->item->adjacent[a]; a++)
					cur_surface += cur2->item->adjacent[a]->surface(coords, repl);
				cur2 = cur2->next;
			}
			cur_surface /= 2;
			//cout << "cur_surface=" << cur_surface << endl;
			found = is_small(cur_surface - init_surface);
			if (!found)
				cur = cur->next;
		}

		if (cur)
		{
			remove_vertex(cur->item);
			edge_to_remove = cur->item;
			to_remove = true;
		}
		else
		{
			cur = edge_list;
			while (cur)
			{
				//cur->item->color = 4;	// Change color of facet near deleted edges
				cur = cur->next;
			}
		}
	}

}


// scan all vertices, if it is unuseful, then remove it.
void stl_file::optim_facets(int nb_face)
{
	//check();
	eps_empile(optim_eps);
	int opt_v = 0;
	stl_vertex_list * cur = vertex_list;
	while (cur)
	{
		//cur->item->dump();
		if (!cur->item->deleted)
			cur->item->scan_edges(nb_face);
		cur = cur->next;
	}

	cur = vertex_list;
	while (cur)
	{
		if (cur->item->deleted)
			opt_v++;
		/*
		//cur->item->dump();
		if (!cur->item->deleted)
		if (cur->item->to_remove)// && opt_v < 2
		//if (cur->item->coords.z_coord() == 25)//!!
		{
			cur->item->remove_vertex(cur->item->edge_to_remove);
//			cur->item->remove_vertex(cur->item->edge_list->item);
			opt_v++;
		}
		*/
		cur = cur->next;
	}
	if (!silent)
		cout << opt_v << " unuseful vertices removed" << endl;
	eps_depile();
#ifdef _DEBUG
	if (!nodebug)
		check();
#endif
}

// scan all triangle facets, if possible, merge two facets to one quadrangle
void stl_file::optim4()
{
	int opt4 = 0;
	stl_edge_list * cur = edge_list;
	while (cur)
	{
		if (cur->item->merge_adj())
		{
			opt4++;
			//check();
		}
		cur = cur->next;
	}
	if (!silent)
		cout << opt4 << " rectangle facets created" << endl;
	//check();
}

// return plane corresponding to facet
void stl_facet::facetplan(stl_v & op, stl_v & vp)
{
	//long orient;
	//scalplan(coords[0],coords[1],coords[2], op, vp, orient);
	op = stl_v_zero;
	for (int i = 0; i < usage; i++)
		op += coords[i];
	op /= usage;
	vp = normal;

}

bool stl_facet::is_border(const stl_v & pt)
{
	bool result = false;
	for (int i = 0; !result && i < usage; i++)
		result = (pt == coords[i]);
	return(result);
}

bool stl_facet::is_adjacent(stl_facet * fct1)
{
	bool result = false;
	for (int i = 0; !result && i < usage; i++)
		for (int j = 0; !result && j < fct1->usage; j++)
		{
			int i1 = (i + 1) % usage;
			int j1 = (j + 1) % fct1->usage;
			result = (coords[i] == fct1->coords[j] && coords[i1] == fct1->coords[j1]) ||
				(coords[i] == fct1->coords[j1] && coords[i1] == fct1->coords[j]);
		}
#ifdef DEBUG
	if (!nodebug)
		cout << "adjacent(" << nbr_facet << "," << fct1->nbr_facet << ") = " << result << endl;
#endif
	return(result);
}

void stl_file::add_prim(stl_prim * prim)
{
	stl_prim_list * prim_elem = new stl_prim_list;
	prim_elem->item = prim;
	prim_elem->next = 0;
	if (prim_list == 0)
		prim_list = prim_elem;
	if (last_prim != 0)
		last_prim->next = prim_elem;
	last_prim = prim_elem;
}

const stl_lin std_lin = stl_lin(stl_v(1.0, 0.0, 0.0),
	stl_v(0.0, 0.0, 1.0),
	stl_v(0.0, 1.0, 0.0),
	stl_v(0.0, 0.0, 0.0));

const int MAX_AXE = 128;

/*
Line type 1's format is:

Line Format:
1 colour x y z a b c d e f g h i part.dat
Fields a through i are orientation & scaling parameters, which can be used in 'standard' 3D transformation matrices. Fields x, y and z also fit into this matrix:

	| a d g 0 |
	| b e h 0 |
	| c f i 0 |
	| x y z 1 |

so that every point (x, y ,z) gets transformed to (x', y', z') :

	x' = a*x + b*y + c*z + x
	y' = d*x + e*y + f*z + y
	z' = g*x + h*y + i*z + z

or, in matrix-math style:
									| a d g 0 |
	| X' Y' Z' 1 | = | X Y Z 1 | x  | b e h 0 |
									| c f i 0 |
									| x y z 1 |
*/



stl_surf::stl_surf()
{
	count = 1;
	prim = 0;
	adj_cyl = false;
}

void stl_surf::write(ostream & os)
{
	if (adj_cyl)
		os << " adjacent to cylinder" << endl;
}

void stl_surf::set_facet(stl_facet * fct)
{
}

stl_cylinder::stl_cylinder(const stl_v & _od, const stl_v & _vd, double _radius, double _height) : stl_surf()
{
	od = _od;
	vd = _vd;
	sidrpl(od, vd, stl_v_zero, vd, od);
	radius = _radius;
	height = _height;
	cw = false;
}

void stl_cylinder::write(ostream & os)
{
	os << setprecision(5);
	os << "0 cylinder " << od << ", " << vd << ", " << Mfround(radius, 001) << ", " << Mfround(height, 001);
	if (adj_cyl)
		os << " adjacent to cylinder";
	os << endl;
}


bool stl_cylinder::compatible(stl_facet * fct, stl_facet * adj)
{
	// Test if the two facets have same axis, radius and height
	bool result = false;
	stl_cylinder * cylinder = fct->new_cylinder(adj);
	if (cylinder)
	{
		result = ((vd == cylinder->vd || vd == -cylinder->vd) &&
			sapptdr(cylinder->od, od, vd) &&
			equal(radius, cylinder->radius) &&
			equal(height, cylinder->height));
		delete cylinder;
		if (result)
		{
			// test CW of cylinder
			stl_v p1 = fct->vertex[0]->coords + fct->normal;
			stl_v p2 = fct->vertex[0]->coords - fct->normal;
			double d1 = spadr(p1, od, vd);
			double d2 = spadr(p2, od, vd);
			if (d1 > d2)
				cw = true;
		}
	}
	return(result);
}

stl_plane::stl_plane(const stl_v & _op, const stl_v & _vp) : stl_surf()
{
	op = _op;
	vp = _vp;
	sidrpl(stl_v_zero, vp, op, vp, op);
}

void stl_plane::write(ostream & os)
{
	os << setprecision(5);
	os << "0 plane " << op << ", " << vp;
	if (adj_cyl)
		os << " adjacent to cylinder";
	os << endl;
}

bool stl_plane::compatible(stl_facet * fct, stl_facet *)
{
	bool result = false;
	stl_plane * plane = fct->new_plane();
	if (plane)
	{
#ifdef DEBUG
		if (!nodebug)
		{
			cout << "compatible [" << op << "," << vp
				<< "] and [" << plane->op << "," << plane->vp << "] "
				<< (vp - plane->vp).norm() << " "
				<< spaplan(plane->op, op, vp) << " = "
				<< (vp == plane->vp && sapptpl(plane->op, op, vp)) << endl;
		}
#endif
		result = (vp == plane->vp && sapptpl(plane->op, op, vp));
		delete plane;
	}
	return(result);
}

stl_tangent::stl_tangent(stl_facet * fct) : stl_surf()
{
	facet = fct;
}

void stl_tangent::set_facet(stl_facet * fct)
{
	facet = fct;
}

bool stl_tangent::compatible(stl_facet * fct, stl_facet *)
{
	bool result = false;
	bool adj = false;
	for (int i = 0; i < facet->usage && !adj; i++)
		if (facet->adjacent[i] == fct)
			adj = true;
	if (adj && facet->angle_between(fct) < ag_lim)
		result = true;
	return(result);
}

void indent(int level)
{
	for (int i = 0; i < level; i++)
		cout << "  ";
}

stl_plane * stl_facet::new_plane()
{
	stl_plane * result = 0;
	if (!cyl_norm || normal == read_normal)
	{
		stl_v f_op;
		stl_v f_vp;
		facetplan(f_op, f_vp);
		result = new stl_plane(f_op, f_vp);
	}
	return(result);
}


stl_cylinder * stl_facet::new_cylinder(stl_facet * fct2, const stl_v & op, const stl_v & vp, const double angle_cyl)
{
	stl_cylinder * result = 0;
	if (usage != 4 || fct2->usage != 4)
		return(result);

	stl_v v_ref = (normal * fct2->normal).normaliz();
	if (!v_ref.null())
	{
		double ag = normal.angle(fct2->normal, v_ref);
#ifdef DEBUG
		if (!nodebug)
			cout << "angle(" << nbr_facet << "," << fct2->nbr_facet << ") = " << degre(ag) << "°" << endl;
#endif
		if (equal(ag, angle_cyl))
		{
			stl_v op2;
			stl_v vp2;
			fct2->facetplan(op2, vp2);

			stl_v od1;
			stl_v vd1;
			double dist;
			stl_v p1;
			stl_v p2;
			int res = s_dist_drdr(op, vp, op2, vp2, dist, p1, p2);
			od1 = (p1 + p2) / 2;
			vd1 = (vp * vp2).normaliz();
			double rad = 0;
			for (int i = 0; i < fct2->usage; i++)
			{
				rad += spadr(fct2->coords[i], od1, vd1);
				//cout << "dist v[" << i << "]=" << fct2->coords[i] << " " << od1 << " " << " " << vd1 << " rad=" << rad << endl;
			}
			rad /= fct2->usage;

			double max = -1e5;
			double min = 1e5;
			for (int k = 0; k < usage; k++)
			{
				double dist = spaplan(coords[k], od1, vd1);
				if (dist < min)
					min = dist;
				if (dist > max)
					max = dist;
			}
			for (int k = 0; k < fct2->usage; k++)
			{
				double dist = spaplan(fct2->coords[k], od1, vd1);
				if (dist < min)
					min = dist;
				if (dist > max)
					max = dist;
			}
			bool ok = true;

			for (int k = 0; k < usage && ok; k++)
			{
				double dist = spaplan(coords[k], od1, vd1);
				ok = equal(dist, min) || equal(dist, max);
			}
			for (int k = 0; k < fct2->usage && ok; k++)
			{
				double dist = spaplan(fct2->coords[k], od1, vd1);
				ok = equal(dist, min) || equal(dist, max);
			}

			if (ok)
				if (angle_cyl == pi / 8)
					result = new stl_cylinder(od1, vd1, rad, max - min);
				else
					result = new stl_cylinder48(od1, vd1, rad, max - min);
			else
				result = 0;
		}
	}
	return(result);
}

stl_cylinder * stl_facet::new_cylinder(stl_facet * fct2, const stl_v & op, const stl_v & vp)
{
	return new_cylinder(fct2, op, vp, pi / 8);
}

stl_cylinder * stl_facet::new_cylinder48(stl_facet * fct2, const stl_v & op, const stl_v & vp)
{
	return new_cylinder(fct2, op, vp, pi / 24);
}

stl_cylinder * stl_facet::new_cylinder(stl_facet *	adj)
{
	stl_cylinder * result = 0;
	if (!cyl_norm || normal != read_normal)
	{
		stl_v op;
		stl_v vp;
		facetplan(op, vp);
		//		vp = read_normal;

		if (adj)
			result = new_cylinder(adj, op, vp);
		else
			for (int i = 0; i < usage && !result; i++)
				if (adjacent[i]/* && !adjacent[i]->surf*/)
				{
					result = new_cylinder(adjacent[i], op, vp);
				}
		if (!result)
		{
			if (adj)
				result = new_cylinder48(adj, op, vp);
			else
				for (int i = 0; i < usage && !result; i++)
					if (adjacent[i]/* && !adjacent[i]->surf*/)
					{
						result = new_cylinder48(adjacent[i], op, vp);
					}
		}
	}
	return(result);
}

stl_tangent * stl_facet::new_tangent()
{
	stl_tangent * result = 0;
	for (int i = 0; i < usage && !result; i++)
		if (adjacent[i])
		{
			if (angle_between(adjacent[i]) < ag_lim)
				result = new stl_tangent(this);
		}
	return(result);
}

void stl_facet::calc_surf(stl_facet * fct, stl_surf_list * & surf_list)
{
	stl_surf * f_surf = 0;
	if (fct)
		f_surf = fct->surf;
	//indent(level);
#ifdef DEBUG
	if (!nodebug)
		cout << "calc_surf(" << nbr_facet << ")" << endl;
#endif
	if (usage == 0 || surf)
		return;

	bool first_facet = false;
	if (f_surf)
	{
		if (f_surf->compatible(this, fct))
		{
			surf = f_surf;
			surf->count++;
		}
	}
	else
	{
		// Look for primitive only if not testing compatibility with other primitive
		if (!no_cylinder)
			surf = new_cylinder(0);
		if (!no_plane && !surf)
			surf = new_plane();
		if (!no_tangent && !surf)
			surf = new_tangent();
		if (surf)
			surf_list = new stl_surf_list(surf, surf_list);
		first_facet = true;
	}
	if (surf)
	{
		for (int i = 0; i < usage; i++)
			if (adjacent[i])
			{
				surf->set_facet(this);
				//cout << "test " << nbr_facet << " " << adjacent[i]->nbr_facet << endl;
				adjacent[i]->calc_surf(this, surf_list);
				//edge[i]->color = 4;
			}
		if (first_facet && surf->count == 1)
		{
			if (surf_list->item == surf)
				surf_list = surf_list->next;
			surf = 0;
		}
	}
}

int stl_facet::side_edge(stl_facet * fct)
{
	int result = 0;
	bool found = false;
	while (result < usage && !found)
	{
		if (adjacent[result] == fct)
			found = true;
		else
			result++;
	}
	if (!found)
		return(-1);
	return(result);
}

int stl_facet::side_edge(stl_edge * edg)
{
	int result = 0;
	bool found = false;
	while (result < usage && !found)
	{
		if (edge[result] == edg)
			found = true;
		else
			result++;
	}
	if (!found)
		return(-1);
	return(result);
}

stl_facet * stl_facet::next_surf(stl_facet * fct)
{
	stl_facet * result = 0;
	for (int i = 0; i < usage && !result; i++)
		if (adjacent[i] && adjacent[i] != fct && adjacent[i]->surf == surf)
			result = adjacent[i];
	return(result);
}

void stl_surf::calc_prim(stl_facet * facet, stl_file & stl)
{
}

int nb_zero(const stl_v & v)
{
	int result = 0;
	for (int i = 0; i < 3; i++)
		result += is_small(v.coord(i)) ? 1 : 0;
	return(result);
}

void stl_cylinder::calc_extrem_partial(
	stl_facet * facet,
	stl_facet * & fct1,
	stl_facet * & fct2,
	stl_facet * & fct3,
	stl_facet * & fct4
)
{
	// Search extremity of partial cylinder
	fct1 = 0;
	fct2 = facet;
	stl_facet * sfct1 = fct1;
	stl_facet * sfct2 = fct2;
	int i = 0;
	int tmp_count = count * nb_arete() / 16;
	while (fct2 && (fct2 != facet || i == 0) && i < count)
	{
		//		cout << "(1) fct2=" << fct2->nbr_facet << endl;
		sfct2 = fct2;
		sfct1 = fct1;
		fct2 = fct2->next_surf(fct1);
		fct1 = sfct2;
		i++;
	}
	if (!fct2)
	{
		fct2 = sfct2;
		fct1 = sfct1;
	}
	// fct4 is the extremity
	fct3 = fct1;
	fct4 = fct2;

#ifdef DEBUG
	if (!nodebug)
		cout << "(2)"
		<< " fct1=" << (fct1 ? fct1->nbr_facet : 0)
		<< " fct2=" << (fct2 ? fct2->nbr_facet : 0)
		<< " fct3=" << (fct3 ? fct3->nbr_facet : 0)
		<< " fct4=" << (fct4 ? fct4->nbr_facet : 0)
		<< endl;
#endif
	// Search next orthogonal facet of cylinder
	sfct1 = fct1;
	sfct2 = fct2;
#ifdef DEBUG
	if (!nodebug)
		cout << "(3)"
		<< " fct1=" << (fct1 ? fct1->nbr_facet : 0)
		<< " fct2=" << (fct2 ? fct2->nbr_facet : 0)
		<< endl;
#endif
	int nb_fac_orth = 3 * nb_arete() / 16;
	for (i = 0; i < nb_fac_orth && fct1; i++)
	{
		sfct1 = fct1;
		sfct2 = fct2;
		fct1 = fct1->next_surf(fct2);
		fct2 = sfct1;
#ifdef DEBUG
		if (!nodebug)
			cout << "(3)"
			<< " fct1=" << (fct1 ? fct1->nbr_facet : 0)
			<< " fct2=" << (fct2 ? fct2->nbr_facet : 0)
			<< endl;
#endif
	}
	if (!fct1)
	{
		fct2 = sfct2;
		fct1 = sfct1;
	}
}

void stl_cylinder::calc_extrem_full(
	stl_facet * facet,
	stl_facet * & fct1,
	stl_facet * & fct2,
	stl_facet * & fct3,
	stl_facet * & fct4
)
{
	//Find best first facet of cylinder
	fct1 = 0;
	fct2 = facet;
	fct3 = facet;
	fct4 = 0;
	stl_facet * sfct1 = fct1;
	stl_facet * sfct2 = fct2;
	stl_facet * sfct3 = fct3;
	stl_facet * sfct4 = fct4;
	stl_v max_p;
	stl_v min_p;
	stl_v max_o;
	stl_v min_o;
	bool best = false;
	int best_nb = 0;

	fct3 = fct3->next_surf(fct4);
	fct4 = sfct3;
	int i = 0;
	do
	{
		calc_center(fct3, fct4, max_o, min_o, max_p, min_p);
#ifdef DEBUG
		if (!nodebug)
			cout << "(2)"
			<< " fct3=" << (fct3 ? fct3->nbr_facet : 0)
			<< " fct4=" << (fct4 ? fct4->nbr_facet : 0)
			<< " nb_zero=" << nb_zero(min_p - min_o) << endl;
#endif
		if (nb_zero(min_p - min_o) > best_nb)
		{
			sfct3 = fct3;
			sfct4 = fct4;
			best_nb = nb_zero(min_p - min_o);
#ifdef DEBUG
			if (!nodebug)
				cout << "(2)"
				<< " fct3=" << (fct3 ? fct3->nbr_facet : 0)
				<< " fct4=" << (fct4 ? fct4->nbr_facet : 0)
				<< " best_nb=" << best_nb << endl;
#endif
		}

		sfct1 = fct4;
		sfct2 = fct3;
		if (fct3)
			fct3 = fct3->next_surf(fct4);
		fct4 = sfct2;
		i++;
	} while (best_nb < 2 && i < 4 && fct3);

	fct3 = sfct3;
	fct4 = sfct4;
	fct1 = fct3;
	fct2 = fct4;

	// Search next orthogonal facet of cylinder
	for (i = 0; i < 6 && fct2; i++)
	{
#ifdef DEBUG
		if (!nodebug)
			cout << "(3)"
			<< " fct1=" << (fct1 ? fct1->nbr_facet : 0)
			<< " fct2=" << (fct2 ? fct2->nbr_facet : 0)
			<< endl;
#endif
		sfct1 = fct1;
		sfct2 = fct2;
		fct2 = fct2->next_surf(fct1);
		fct1 = sfct2;
	}
}
void stl_cylinder::calc_center(
	stl_facet * fct3,
	stl_facet * fct4,
	stl_v & max_o,
	stl_v & min_o,
	stl_v & max_p,
	stl_v & min_p
)

{
	//cout << "calc_center( fct3=" << (fct3 ? fct3->nbr_facet : 0) << " fct4=" << (fct4 ? fct4->nbr_facet : 0) << ")" << endl;
	// Calculate projection of facet on cylinder axis
	// to get cylinder height

	double max = -1e5;
	double min = 1e5;
	//	int i = fct3->side_edge(fct3->next_surf(fct4));
	int i = fct4->side_edge(fct3);
	for (int j = 2; j < 4; j++)
	{
		int k = (i + j) % fct4->usage;
		double dist = spaplan(fct4->coords[k], od, vd);
		if (dist > max)
		{
			max = dist;
			max_p = fct4->coords[k];
			sidrpl(od, vd, fct4->coords[k], vd, max_o);
		}
		if (dist < min)
		{
			min = dist;
			min_p = fct4->coords[k];
			sidrpl(od, vd, fct4->coords[k], vd, min_o);
		}
	}
	// Calculate the orthogonal vector
	/*
	stl_v min_p2 = (max_o - min_o) * (max_p - min_o);
	min_p2.normalize();
	min_p2 *= min_o.distance(max_p);
	if (min_p2.ps(min_p - min_o) < 0)
		min_p = -(min_o + min_p2);
	else
		min_p = (min_o + min_p2);
		*/
}

void stl_cylinder::calc_prim(stl_facet * facet, stl_file & stl)
{
	int nb_cyl = 0;
	int primitives = 0;
	stl_v max_p;
	stl_v min_p;
	stl_v max_o;
	stl_v min_o;
	stl_v min_z;
	stl_v max_z;
	//	stl_facet * fct0;
	int tmp_count = count * 16 / nb_arete();

	if (tmp_count >= 4)
	{
		// search first facets of cylinder in following order :
		// fct1, fct2, fct3, fct4
		stl_facet * fct1 = 0;
		stl_facet * fct2 = 0;
		stl_facet * fct3 = 0;
		stl_facet * fct4 = 0;
		int i = 0;
		if (tmp_count < 16)
			calc_extrem_partial(facet, fct1, fct2, fct3, fct4);
		else
			calc_extrem_full(facet, fct1, fct2, fct3, fct4);
#ifdef DEBUG
		if (!nodebug)
		{
			cout << "(4)"
				<< " fct1=" << (fct1 ? fct1->nbr_facet : 0)
				<< " fct2=" << (fct2 ? fct2->nbr_facet : 0)
				<< " fct3=" << (fct3 ? fct3->nbr_facet : 0)
				<< " fct4=" << (fct4 ? fct4->nbr_facet : 0)
				<< endl;
			fct1->dump();
			fct2->dump();
			fct3->dump();
			fct4->dump();
		}
#endif
		if (!fct2)
			return;

		calc_center(fct3, fct4, max_o, min_o, max_p, min_p);

		// Calculate projection of facet on cylinder axis
		// to get cylinder width
		double max = -1e5;
		double min = 1e5;
		i = fct1->side_edge(fct2);
		for (int j = 2; j < 4; j++)
		{
			int k = (i + j) % fct1->usage;
			double dist = spaplan(fct1->coords[k], od, vd);
			//cout << "k=" << k << " coords=" << fct1->coords[k] << " dist=" << dist << endl;
			if (dist < min)
			{
				min = dist;
				min_z = fct1->coords[k];
				//cout << "min=" << min << " min_z=" << min_z << endl;
			}
			if (dist > max)
				max = dist;
		}
		//cout << "min=" << min <<  " max=" << max << endl;

				// Verify that cylinder is not truncated
		stl_facet * lst = facet;
		stl_facet * cur = lst->next_surf(0);
		do
		{
			double lmax = -1e5;
			double lmin = 1e5;
			for (int j = 0; j < 4; j++)
			{
				double dist = spaplan(cur->coords[j], od, vd);
				if (dist < lmin)
					lmin = dist;
				if (dist > lmax)
					lmax = dist;
			}
			//cout << "lmin=" << lmin << " lmax=" << lmax << endl;
			if (!equal(min, lmin) || !equal(max, lmax))
				return;
			stl_facet * nxt = cur->next_surf(lst);
			lst = cur;
			cur = nxt;
		} while (cur && cur != facet);

		if (rand_col_adj_cyl)
		{
			// Set adjacent surface
			lst = facet;
			cur = lst->next_surf(0);
			do
			{
				for (int j = 0; j < cur->usage; j++)
				{
					if (cur->adjacent[j] && cur->adjacent[j]->surf && cur->adjacent[j]->surf != this)
						cur->adjacent[j]->surf->adj_cyl = true;
				}
				stl_facet * nxt = cur->next_surf(lst);
				lst = cur;
				cur = nxt;
			} while (cur && cur != facet);
		}

		// Calculate transformation matrix
		stl_lin std_lin = stl_lin(stl_v(1.0, 0.0, 0.0),
			stl_v(0.0, 1.0, 0.0),
			stl_v(0.0, 0.0, 1.0),
			stl_v(0.0, 0.0, 0.0));
		stl_lin lin_dec = stl_lin(stl_v(1.0, 0.0, 0.0),
			stl_v(0.0, 1.0, 0.0),
			stl_v(0.0, 0.0, 1.0),
			stl_v(0.0, 0.0, 0.0));

		stl_v min_z2 = (max_o - min_o) * (min_p - min_o);
		min_z2.normalize();
		min_z2 *= min_o.distance(min_p);
		if (min_z2.ps(min_z - min_o) < 0)
			min_z = (min_o - min_z2);
		else
			min_z = (min_o + min_z2);

		stl_lin lin_to(max_o - min_o, min_z - min_o, min_p - min_o, min_o);
		int dec_count_tot = 0;

		while (tmp_count > 0)
		{
			string prefix("");
			int dec_count = 0;
			//if (tmp_count ==  4 || tmp_count ==  8 || tmp_count ==  12 || tmp_count ==  16)
			if (tmp_count >= 4 && tmp_count != 6 && tmp_count != 7)
			{
				prefix = string(1, (char)('0' + tmp_count / 4)) + string("-4");
				dec_count = int(tmp_count / 4) * 4;
			}
			else if (tmp_count >= 2)
			{
				prefix = string(1, (char)('0' + tmp_count / 2)) + string("-8");
				dec_count = int(tmp_count / 2) * 2;
			}
			else if (tmp_count == 1)
			{
				prefix = "1-16";
				dec_count = 1;
			}
			int col = col_line1;
			if (rand_col_surfaces)
				col = random_color();
			if (nb_arete() == 48)
				prefix = "48\\" + prefix;
			stl_prim * prim_cyl = new stl_prim(prefix + "cyli", col, cw);
			stl_prim * prim_ed1 = new stl_prim(prefix + "edge", 16);
			stl_prim * prim_ed2 = new stl_prim(prefix + "edge", 16);
			prim_cyl->transfo = stl_lin(std_lin, stl_lin(min_p - min_o, max_o - min_o, min_z - min_o, min_o)) * lin_dec;
			prim_ed1->transfo = stl_lin(std_lin, stl_lin(min_p - min_o, (max_o - min_o).normaliz(), min_z - min_o, min_o)) * lin_dec;
			prim_ed2->transfo = stl_lin(std_lin, stl_lin(max_p - max_o, (max_o - min_o).normaliz(), min_z - min_o, max_o)) * lin_dec;
			stl.add_prim(prim_cyl);
			stl.add_prim(prim_ed1);
			stl.add_prim(prim_ed2);
			prim = prim_cyl;
			tmp_count -= dec_count;
			dec_count_tot += dec_count;
			lin_dec = stl_lin(stl_v(cos(dec_count_tot * pi / 8), 0.0, sin(dec_count_tot * pi / 8)),
				stl_v(0.0, 1.0, 0.0),
				stl_v(cos((dec_count_tot + 4) * pi / 8), 0.0, sin((dec_count_tot + 4) * pi / 8)),
				stl_v(0.0, 0.0, 0.0));
		}
		primitives++;

		//cout << prim_cyl.transfo << endl;
	}
#ifdef DEBUG
	if (!nodebug)
	{
		cout << "stl_cylinder::calc_prim max_o=" << max_o - min_o
			<< " min_z=" << min_z - min_o
			<< " min_p=" << min_p - min_o
			<< " min_o=" << min_o
			<< endl;
		cout << "2 4 " << min_p << " " << min_o << endl;
		cout << "2 2 " << max_o << " " << min_o << endl;
		cout << "2 1 " << min_z << " " << min_o << endl;
	}
#endif
}

void stl_surf::calc_curv(stl_file & stl)
{
}

void stl_cylinder::calc_curv(stl_file & stl)
{
}

void stl_plane::calc_curv(stl_file & stl)
{
}

void stl_file::primitives()
{
	//	eps_empile(0.01);
		// Bind facet together to faces
	stl_facet_list * cur = facet_list;
	while (cur)
	{
		stl_facet * fct = cur->item;
		if (fct->usage && !fct->surf)
		{
			fct->calc_surf(0, surf_list);
		}
		cur = cur->next;
	}
#ifdef OPTIM_CURVE
	// Search curves
	stl_surf_list * sur_ptr = surf_list;
	while (sur_ptr)
	{
		stl_surf * surf = sur_ptr->item;
		surf->calc_curv(*this);
		sur_ptr = sur_ptr->next;
	}
	// Search circle primitives
#endif

	if (!noprim)
	{
		// Search primitives
		cur = facet_list;
		while (cur)
		{
			stl_facet * fct = cur->item;
			if (fct->usage && fct->surf && !fct->surf->prim)
			{
				fct->surf->calc_prim(fct, *this);
			}
			cur = cur->next;
		}
		cur = facet_list;
		while (cur)
		{
			stl_facet * fct = cur->item;
			if (fct->usage && fct->surf && fct->surf->prim)
				fct->surf->clean(fct);
			cur = cur->next;
		}
	}
	//	eps_depile();
}

void stl_facet::clean()
{
	for (int i = 0; i < usage; i++)
		edge[i]->printed = true;
	printed = true;
}

void stl_surf::clean(stl_facet * fct)
{
	for (int i = 0; i < fct->usage; i++)
		fct->edge[i]->printed = true;
	fct->printed = true;
}

void stl_cylinder::clean(stl_facet * fct)
{
	for (int i = 0; i < fct->usage; i++)
	{
		if ((fct->adjacent[i] && fct->adjacent[i]->surf == fct->surf) ||
			!fct->edge[i]->colinear(vd))
			fct->edge[i]->printed = true;
	}
	fct->printed = true;
}

double stl_facet::angle_between(stl_facet * fct)
{
	if (normal == fct->normal)
		return(0);
	stl_v v_ref = (normal * fct->normal).normaliz();
	if (v_ref.null())
		return(pi);
	else
	{
		double ag = (normal).angle(fct->normal, v_ref);
		return(ag);
	}
}

bool stl_edge::calc_edge()
{
	bool result = false;
	if (!adjacent[0] || !adjacent[1])
	{
		// border edge
		linetype = 0;//2
	}
	else if (adjacent[0]->normal == adjacent[1]->normal)
	{
		// same plane
		linetype = 0;
		result = true;
	}
	else
	{
		double angle = adjacent[0]->angle_between(adjacent[1]);
		if (angle > ag_lim)
			linetype = 2;
		else
		{
			// Check if the optional edge is external or internal
			if (opt5)
			{
				stl_v dir_e = (coord(1) - coord(2)).normaliz();
				int side = dir_e.ps(adjacent[0]->normal + adjacent[1]->normal) > 0;
				if (side)
					linetype = 5;
				else
					linetype = 0;
			}
			else
			{
				linetype = 5;
			}
		}
	}
	return(result);
}

// Calculate edges between 2 adjacent facets
void stl_file::calc_edge()
{
	//	eps_empile(0.0001);
	int nb_edge = 0;
	stl_edge_list * cur = edge_list;
	while (cur)
	{
		if (cur->item->linetype >= 0 && !cur->item->calc_edge())
			nb_edge++;
		cur = cur->next;
	}
	//	eps_depile();
	if (!silent)
		cout << nb_edge << " edges cleared" << endl;
}
// Calculate edges between 2 adjacent facets
	/*!!
void stl_file::calc_edge()
{
	/*!!
	int nb_edge = 0;
	stl_facet_list * cur = facet_list;
	stl_facet * fct1;
	while (cur)
	{
		fct1 = &cur->item;
		if (fct1->usage)
		{
			stl_facet_list * cur2 = cur->next;
			while (cur2)
			{
				stl_facet * fct2 = &cur2->item;
				if (fct2->usage)
					nb_edge += create_edge(fct1,fct2);
				cur2 = cur2->next;
			}
		}
		cur = cur->next;
	}
	cout << nb_edge << " edges created" << endl;
}
	*/


double get_arg(char* argv[], int idx)
{
	double result = 0;
	stringstream s;
	s << argv[idx];
	s >> result;
	//sscanf(argv[idx],"%g",&result);
	return(result);
}

int iget_arg(char* argv[], int idx)
{
	int result = 0;
	stringstream s;
	s << argv[idx];
	s >> result;
	//sscanf(argv[idx],"%g",&result);
	return(result);
}

int usage()
{
	cout << "stl2dat v" << version << endl;
	cout << "Usage : stl2dat filename [] [name [author]]" << endl;
	cout << "        -out filename : name of the output file" << endl;
	cout << "        -scalemm : units of input file are in mm instead of ldraw units" << endl;
	cout << "        -scale S : S is scale factor" << endl;
	cout << "        -o X Y Z : origin point" << endl;
	cout << "        -m X Y Z A B C D E F G H I : transformation matrix" << endl;
	cout << "        -a angle : angle is limit for optional edges" << endl;
	cout << "        -aq angle : angle is limit for quadrangles" << endl;
	cout << "        -at angle : angle is limit for remove unuseful facets" << endl;
	cout << "        -eps value" << endl;
	cout << "        -teps value : minimum distance for coincidence" << endl;
	cout << "        -deps value : maximum determinant for coplanar quads" << endl;
	cout << "        -oeps value : epsilon used for remove unuseful facets" << endl;
	//		cout << "        -o3 : merge triangles" << endl;
	cout << "        -o4 : no merge to quadrangles" << endl;
	cout << "        -o5 : remove of internal (concave) optional edges" << endl;
	cout << "        -op : no primitives calculation" << endl;
	cout << "        -of : no remove of unuseful facets" << endl;
	cout << "        -pp : print geometric surfaces" << endl;
	cout << "        -cn : detect cylinder with normal difference" << endl;
	cout << "        -oe : no creation of edges" << endl;
	cout << "        -np : no plane geometric surfaces" << endl;
	cout << "        -nc : no cylinder geometric surfaces" << endl;
	cout << "        -nt : no tangent geometric surfaces" << endl;
	cout << "        -c1 color : color of primitives" << endl;
	cout << "        -c2 color : color of edges" << endl;
	cout << "        -c3 color : color of triangles" << endl;
	cout << "        -c4 color : color of quadrangles" << endl;
	cout << "        -c5 color : color of optional edges" << endl;
	cout << "        -no1 : no primitives" << endl;
	cout << "        -no2 : no edges" << endl;
	cout << "        -no3 : no triangles" << endl;
	cout << "        -no4 : no quadrangles" << endl;
	cout << "        -no5 : no optional edges" << endl;
	cout << "        -cs : random color for surfaces" << endl;
	cout << "        -cc : random color for adjacent to cylinder" << endl;
	cout << "        -dat : produces a .dat file" << endl;
	cout << "        -ldr : produces a .ldr file" << endl;
	cout << "        -nobfc : no BFC instructions" << endl;
	cout << "        -silent : does not display information (DEFAULT)" << endl;
	cout << "        -verbose : display information" << endl;
	cout << "        -ldraw : does not add ldraw standard header" << endl;
	cout << "        -raw : conversion without any optimisation or edge calculation" << endl;
	cout << "        -h, --help : show this text and exit" << endl;
	return 0;
}

int main(int argc, char* argv[])
{
	if (argc < 2 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))
	{
		return usage();
	}
#ifdef DEBUG
	if (!nodebug)
		cout << "File:" << argv[1] << endl;
#endif

	bool opt3 = false;
	bool opt4 = true;
	opt5 = false;
	bool opt_prim = true;
	bool opt_edge = true;
	bool opt_facets = true;

	ag_lim = radian(22.51);
	ag_lim_q = pi / 20;
	ag_lim_t = ag_lim / 2;
	bool b_ag_lim_t = false;
	bool auth_done = false;
	bool name_done = false;
	bool eps_done = false;
	bool print_geom = false;
	no_plane = false;
	no_cylinder = false;
	no_tangent = false;
	cyl_norm = false;
	col_line1 = 16;
	col_line2 = 24;
	col_line3 = 16;
	col_line4 = 16;
	col_line5 = 24;
	no_line1 = false;
	no_line2 = false;
	no_line3 = false;
	no_line4 = false;
	no_line5 = false;
	rand_col_surfaces = false;
	ldr_opt = false;
	ldr_out = false;
	bfc = true;
	silent = true;
	ldraw = true;
	stl_raw = false;
	determinant_eps = 0.1;

	noprim = false;
	debug = false;
	nodebug = false;

	eps_empile(0.01);
	double topo_eps = 0.001;
	optim_eps = 0.001;

	if (access(argv[1], 4) == -1)
	{
		cout << "File " << argv[1] << " does not exist." << endl;
		exit(-1);
	}

	stl_file stl(argv[1]);

	string attr;
	string arg_partname;
	int idx = 2;
	while (idx < argc)
	{
		attr = argv[idx];
		if (attr == "-h" || attr == "--help"){
			return usage();
		}
		else if (attr == "-o")
		{
			double c[3];
			for (int i = 0; i < 3; i++)
			{
				idx++;
				c[i] = 0;
				if (idx < argc)
					c[i] = get_arg(argv, idx);
			}
			stl.origin = stl_v(c[0], c[1], c[2]);
			stl.transf = stl_lin(
				stl_v(1.0, 0.0, 0.0),
				stl_v(0.0, 1.0, 0.0),
				stl_v(0.0, 0.0, 1.0),
				stl.origin);
		}
		else if (attr == "-m")
		{
			double c[12];
			for (int i = 0; i < 12; i++)
			{
				idx++;
				c[i] = 0;
				if (idx < argc)
					c[i] = get_arg(argv, idx);
			}
			stl.transf = stl_lin(
				stl_v(c[3], c[4], c[5]),
				stl_v(c[6], c[7], c[8]),
				stl_v(c[9], c[10], c[11]),
				stl_v(c[0], c[1], c[2]));
			stl.origin = stl_v(c[0], c[1], c[2]);
		}
		else if (attr == "-scale")
		{
			double s;
			idx++;
			if (idx < argc)
				s = get_arg(argv, idx);
			stl.transf = stl_lin(
				stl_v(s, 0.0, 0.0),
				stl_v(0.0, s, 0.0),
				stl_v(0.0, 0.0, s),
				stl.origin);
		}
		else if (attr == "-scalemm")
		{
			double s = 20.0 / 8;	// mm to ldraw units
			stl.transf = stl_lin(
				stl_v(s, 0.0, 0.0),
				stl_v(0.0, s, 0.0),
				stl_v(0.0, 0.0, s),
				stl.origin);
		}
		else if (attr == "-a")
		{
			idx++;
			if (idx < argc)
			{
				ag_lim = fabs(radian(get_arg(argv, idx)));
				if (is_small(ag_lim))
					ag_lim = pi * 0.01;
				if (!b_ag_lim_t)
					ag_lim_t = ag_lim / 2;
			}
		}
		else if (attr == "-aq")
		{
			idx++;
			if (idx < argc)
			{
				ag_lim_q = fabs(radian(get_arg(argv, idx)));
				if (is_small(ag_lim_q))
					ag_lim_q = pi * 0.01;
			}

		}
		else if (attr == "-at")
		{
			idx++;
			if (idx < argc)
			{
				ag_lim_t = fabs(radian(get_arg(argv, idx)));
				if (is_small(ag_lim_t))
					ag_lim_t = pi * 0.01;
				b_ag_lim_t = true;
			}

		}
		else if (attr == "-eps")
		{
			idx++;
			if (idx < argc)
			{
				double eps = get_arg(argv, idx);
				if (!eps_done)
				{
					eps_empile(eps);
					eps_done = true;
				}
			}
		}
		else if (attr == "-teps")
		{
			idx++;
			if (idx < argc)
			{
				topo_eps = get_arg(argv, idx);
			}
		}
		else if (attr == "-deps")
		{
			idx++;
			if (idx < argc)
			{
				determinant_eps = get_arg(argv, idx);
			}
		}
		else if (attr == "-oeps")
		{
			idx++;
			if (idx < argc)
			{
				optim_eps = get_arg(argv, idx);
			}
		}
		else if (attr == "-c1")
		{
			idx++;
			if (idx < argc)
				col_line1 = iget_arg(argv, idx);
		}
		else if (attr == "-c2")
		{
			idx++;
			if (idx < argc)
				col_line2 = iget_arg(argv, idx);
		}
		else if (attr == "-c3")
		{
			idx++;
			if (idx < argc)
				col_line3 = iget_arg(argv, idx);
		}
		else if (attr == "-c4")
		{
			idx++;
			if (idx < argc)
				col_line4 = iget_arg(argv, idx);
		}
		else if (attr == "-c5")
		{
			idx++;
			if (idx < argc)
				col_line5 = iget_arg(argv, idx);
		}
		else if (attr == "-no1")
			no_line1 = true;
		else if (attr == "-no2")
			no_line2 = true;
		else if (attr == "-no3")
			no_line3 = true;
		else if (attr == "-no4")
			no_line4 = true;
		else if (attr == "-no5")
			no_line5 = true;
		//		else if (attr == "-o3")
		//			opt3 = false;
		else if (attr == "-o4")
			opt4 = false;
		else if (attr == "-o5")
			opt5 = true;
		else if (attr == "-op")
			opt_prim = false;
		else if (attr == "-oe")
			opt_edge = false;
		else if (attr == "-of")
			opt_facets = false;
		else if (attr == "-pp")
			print_geom = true;
		else if (attr == "-nc")
			no_cylinder = true;
		else if (attr == "-np")
			no_plane = true;
		else if (attr == "-nt")
			no_tangent = true;
		else if (attr == "-cn")
			cyl_norm = true;
		else if (attr == "-cs")
			rand_col_surfaces = true;
		else if (attr == "-cc")
			rand_col_adj_cyl = true;
		else if (attr == "-silent")
			silent = true;
		else if (attr == "-verbose")
			silent = false;
		else if (attr == "-ldraw")
			ldraw = false;
		else if (attr == "-raw")
			stl_raw = true;
		else if (attr == "-noprim")
			noprim = true;
		else if (attr == "-ldr")
			ldr_opt = true;
		else if (attr == "-dat")
			ldr_opt = false;
		else if (attr == "-nobfc")
			bfc = false;
		else if (attr == "-debug")
			debug = true;
		else if (attr == "-nodebug")
			nodebug = true;
		else if (attr == "-out")
		{
			idx++;
			if (idx < argc)
			{
				stl.out_filename = argv[idx];
				ldr_out = true;
			}
		}
		else if (!name_done && argv[idx][0] != '-')
		{
			arg_partname = argv[idx];
			name_done = true;
		}
		else if (!auth_done && argv[idx][0] != '-')
		{
			stl.author = argv[idx];
			auth_done = true;
		}
		idx++;
	}
	if (!auth_done && idx < argc)
		stl.author = argv[idx];

	if (!silent)
		cout << "stl2dat v" << version << " (c) Marc Klein" << endl;

	eps_empile(optim_eps);
	stl.read();
	eps_depile();

	if (stl.nbr_facet <= 0)
	{
		if (!silent)
			cout << "File " << argv[1] << " is empty." << endl;
		exit(0);
	}

	if (name_done)
		stl.partname = arg_partname;
	else
		stl.partname = stl.get_desc();

	if (!stl_raw)
	{
		eps_empile(topo_eps);
		stl.topo();
		eps_depile();

		/**/
//stl.dump();
		//stl.check();

		if (opt_facets)
		{
			//cout << "optim_facets ag_lim_t=" << degre(ag_lim_t) << endl;
			stl.optim_facets(2);
			stl.optim_facets(1);
		}
		//stl.dump();
		//stl.check();
				/**/

		if (opt4)
			stl.optim4();
		//stl.check();
		eps_empile(optim_eps);
		if (opt_prim)
			stl.primitives();
		eps_depile();
		if (opt_edge)
			stl.calc_edge();
		//		if (opt3)
		//			stl.optim3();
	}

	stl.write_dat(print_geom);

	return 0;
}

