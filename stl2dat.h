#include <string>
#include <iostream>
#include "ort.h"

#define LIST_OF(cl) \
class cl##_list \
{ \
public: \
	cl * item; \
	cl##_list * next; \
	cl##_list(); \
	cl##_list( cl * it, cl##_list * ne ) { item = it; next = ne; }; \
	void list_append( cl * it, cl##_list * & list, cl##_list * & last ); \
}; \
cl##_list::cl##_list() \
{ \
	next = 0; \
} \
void cl##_list::list_append( cl * it, cl##_list * & list, cl##_list * & last ) \
{ \
	cl##_list * elem = new cl##_list; \
	elem->item = it; \
	elem->next = 0; \
	if (list == 0) \
		list = elem; \
	if (last != 0) \
		last->next = elem; \
}


/*
The standard ldraw primitives
*/
class stl_prim
{
public:
	stl_lin transfo;
	std::string name;
	int color;
	bool cw;	// for bfc

	stl_prim() {};
	stl_prim(int col)
	{
		color = col;
	};
	stl_prim(const std::string & nm, int col, bool _cw = true)
	{
		name = nm;
		color = col;
		cw = _cw;
	};
	void write(std::ostream & os);
};

LIST_OF(stl_prim)

class stl_file;
class stl_facet;
class stl_facet_list;
class stl_curv;
class stl_curv_list;

enum surf_type
{
	t_surf,
	t_plane,
	t_cylinder,
	t_tangent
};

/*
The class for geometric primitives
*/
class stl_surf
{
public:
	int count;
	stl_prim * prim;
	bool adj_cyl;

	stl_surf();
	virtual surf_type type() { return(t_surf); };

	virtual bool compatible( stl_facet * fct, stl_facet * adj) = 0;
	virtual void write( std::ostream & os);
	virtual void set_facet( stl_facet * fct);
	virtual void calc_prim( stl_facet * facet, stl_file & stl);
	virtual void clean(stl_facet * fct);
	virtual void calc_curv( stl_file & stl);
};

class stl_cylinder : public stl_surf
{
public:
	stl_v od;
	stl_v vd;
	double radius;
	double height;
	stl_v normal;
	bool cw;	// for bfc

	stl_cylinder( const stl_v & od, const stl_v & vd, double radius, double height);
	virtual surf_type type() { return(t_cylinder); };

	virtual void write( std::ostream & os);
	virtual bool compatible( stl_facet * fct, stl_facet * adj);
	virtual void calc_prim( stl_facet * facet, stl_file & stl);
	virtual void clean(stl_facet * fct);

	void calc_extrem_partial(
			stl_facet * facet,
			stl_facet * & fct1,
			stl_facet * & fct2,
			stl_facet * & fct3,
			stl_facet * & fct4
			);
	void calc_extrem_full(
			stl_facet * facet,
			stl_facet * & fct1,
			stl_facet * & fct2,
			stl_facet * & fct3,
			stl_facet * & fct4
			);
	void calc_center(
			stl_facet * fcta,
			stl_facet * fctb,
			stl_v & max_o,
			stl_v & min_o,
			stl_v & max_p,
			stl_v & min_p
			);
	virtual void calc_curv( stl_file & stl);
	virtual const int nb_arete() { return 16; }
	virtual const double angle_cyl() { return pi/8; }
};

class stl_cylinder48 : public stl_cylinder
{
public:
	stl_cylinder48( const stl_v & od, const stl_v & vd, double radius, double height) : 
	  stl_cylinder( od, vd,  radius, height) {};
	virtual const int nb_arete() { return 48; }
	virtual const double angle_cyl() { return pi/24; }
};

// Plane
class stl_plane : public stl_surf
{
public:
	stl_v op;	// Origin of plane
	stl_v vp;	// Direction vector of plane

	stl_plane( const stl_v & op, const stl_v & vp);
	virtual surf_type type() { return(t_plane); };

	virtual bool compatible( stl_facet * fct, stl_facet *);
	virtual void write( std::ostream & os);
	virtual void calc_curv( stl_file & stl);
};

// Tangent facet
class stl_tangent : public stl_surf
{
public:
	stl_facet * facet;	// Tangent facet

	stl_tangent( stl_facet * facet);
	virtual surf_type type() { return(t_tangent); };

	virtual bool compatible( stl_facet * fct, stl_facet *);
	virtual void set_facet( stl_facet * fct);
};

LIST_OF(stl_surf)

class stl_edge;
class stl_vertex;

// Curves : list of 
class stl_curv
{
public:
	int count;
	stl_prim * prim;
	bool adj_cyl;

	stl_curv();

	virtual bool compatible( stl_edge * edg, stl_edge * adj) = 0;
	virtual void write( std::ostream & os);
	virtual void set_facet( stl_edge * edg);
	virtual void calc_prim( stl_edge * facet, stl_file & stl);
	virtual void clean(stl_edge * edg);
};

LIST_OF(stl_curv)

/*
Store information for a facet with 3 or 4 points
*/
class stl_facet
{
public:
	int usage;	// 0 not used
				// 3 edges
				// 4 edges
	int nbr_facet;		// number of facet in the solid
	int color;			// color
	stl_v normal;		// normal calculated
	stl_v read_normal;	// normal read in stl file
	stl_v coords[4];		// coordinates of vertices
	stl_edge * edge[4];			// the edges defining the facet
	stl_facet * adjacent[4];	// the facet adjacent on edge[n]
	stl_vertex * vertex[4];		// the vertices
	stl_surf * surf;			// primitive surface
	bool printed;
	double m_surface;
	int precision;

	std::string ident();
	void dump();

	stl_facet();
	void write(std::ostream & os, int col3 = 16, int col4 = 16, int edgecolor = 24, bool printface = true, bool printedge = true);
	void facetplan( stl_v & oq, stl_v & wq);
	void calculate_normal();
	void clean();
	bool is_convex(bool trace = false);
	void calcul_angles( double & totag, double & minag, double & maxag);
	double angle_between(stl_facet * fct);
	bool is_border(const stl_v & pt);
	bool is_adjacent(stl_facet * fct1);
	void merge_facet(stl_edge * m_edge, stl_facet * m_facet, int ve1, int ve2);
	void calc_surf(stl_facet * fct, stl_surf_list * & surf_list);
	void replace_vertex(stl_vertex * ve_ori, stl_vertex * ve_repl);
	void replace_adjacent(stl_vertex * ve_ori, stl_edge * ed, bool del_edge = true);
	bool recalculate_precision();

	stl_facet * next_surf( stl_facet * fct);
	int side_edge( stl_facet * fct);
	int side_edge( stl_edge * edge);
	stl_plane * new_plane();
	stl_cylinder * new_cylinder( stl_facet * adj);
	stl_cylinder * new_cylinder( stl_facet * fct2, const stl_v & op, const stl_v & vp, const double angle_cyl);
	stl_cylinder * new_cylinder( stl_facet * fct2, const stl_v & op, const stl_v & vp);
	stl_cylinder * new_cylinder48( stl_facet * fct2, const stl_v & op, const stl_v & vp);
	stl_tangent * new_tangent();
	bool check();
	double surface();
	double surface(const stl_v & v_ori, const stl_v & v_repl);
	double determinant(int precision = -1);
};

LIST_OF(stl_facet)

/*
Store information for an edge between 2 facets
*/
class stl_edge
{
public:
	stl_facet * adjacent[2];	// The 2 facets adjacent
	int linetype;				// ldraw line type 2 or 5
	int coords[2];				// number of vertex in facet adjacent[i]
	stl_vertex * vertex[2];		// the vertices
	bool printed;
	bool deleted;
	int color;
	int nbr_edge;

	stl_edge();

	stl_v coord(int i, int face = 0);
	stl_vertex * vertx(int i, int face = 0);
	void set_coord(int i, int face, const stl_v & val);
	void write(std::ostream & os, int col2 = 24, int col5 = 24);
	std::string ident();
	bool merge_adj();
	int side(stl_facet * m_facet);
	void reaffect(stl_facet * m_facet, stl_facet * repl_facet, int ve);
	bool calc_edge();
	bool check();
	void dump();
	bool colinear( const stl_v & vd);
	stl_edge * next(int & vx, stl_facet * & s_fct);
};

LIST_OF(stl_edge)

/*
Store information for an vertex
*/
class stl_vertex
{
public:
	stl_v coords;
	stl_edge_list * edge_list;
	bool deleted;
	bool to_remove;
	int precision;
	stl_edge * edge_to_remove;

	stl_vertex(const stl_v & _coords);
	void add_edge( stl_edge * edge);
	bool check();
	void dump();
	void scan_edges(int nb_face);
	void remove_vertex(stl_edge * ed_supr);

    friend std::ostream & operator << ( 
		std::ostream & s, 
		const	stl_vertex * v 
	);
    /*-
     ! Use:	Write a given point/vector in an ostream.		
     ! remark:	Format is < x1 , x2 , x3 >		
     */
};

LIST_OF(stl_vertex)

class stl_file
{
	std::ifstream is;
	std::string desc;
	std::string name;
public:
	std::string partname;
	std::string author;
	std::string out_filename;
	int ok;
	stl_v origin;
	stl_lin transf;
	bool swapyz = false;
	int nbr_facet;
private:
	stl_facet_list * facet_list;
	stl_facet_list * last_facet;
	stl_edge_list * edge_list;
	stl_edge_list * last_edge;
	stl_vertex_list * vertex_list;
	stl_prim_list * prim_list;
	stl_prim_list * last_prim;
	stl_surf_list * surf_list;
	stl_curv_list * curv_list;
	bool stl_bin;

public:
	stl_file(const char * filename);
	void dump();

	bool read();
	void write_dat(bool print_geom);
	bool check();
	void topo();
	void optim3();
	void optim4();
	void optim_facets(int nb_face);
	void primitives();
	void calc_edge();
	std::string get_desc() 
	{
		return(desc);
	};
	void add_prim( stl_prim * prim);

private:
	bool token(const char * token);
	bool read_header();
	int read_facet_text();
	bool read_facet_bin();
	bool new_facet(const stl_v & normal, const stl_v vertex[]);
	int create_edge(stl_facet * fct1, stl_facet * fct2);
	void add_edge(int usage, const stl_v& v0, const stl_v& v1, const stl_v& v2 = stl_v(0,0,0), const stl_v& v3 = stl_v(0,0,0));

};

