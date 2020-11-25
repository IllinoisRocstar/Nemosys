#include "libsupermesh.h"

#define TRI_BUF_SIZE 22
#define TET_BUF_SIZE 81

#ifdef LIBSUPERMESH_DOUBLE_PRECISION
typedef double c_kind;
#else
typedef float c_kind;
#endif

#ifdef __cplusplus
#include <cmath>
extern "C" {
#else
#include <math.h>
#endif
void libsupermesh_sort_intersection_finder_set_input(long* nnodes_a, int* dim_a, long* nelements_a, int* loc_a, long* nnodes_b, int* dim_b, long* nelements_b, int* loc_b, c_kind* positions_a, long* enlist_a, c_kind* positions_b, long* enlist_b);
void libsupermesh_sort_intersection_finder_query_output(long* nindices);
void libsupermesh_sort_intersection_finder_get_output(long* nelements, long* nindices, long* indices, long* ind_ptr);
void libsupermesh_tree_intersection_finder_set_input(long* nnodes_a, int* dim_a, long* nelements_a, int* loc_a, long* nnodes_b, int* dim_b, long* nelements_b, int* loc_b, c_kind* positions_a, long* enlist_a, c_kind* positions_b, long* enlist_b);
void libsupermesh_tree_intersection_finder_query_output(long* nindices);
void libsupermesh_tree_intersection_finder_get_output(long* nelements, long* nindices, long* indices, long* ind_ptr);
void libsupermesh_intersect_tris_real(c_kind* tri_a, c_kind* tri_b, c_kind* tris_c, int* n_tris_c);
void libsupermesh_intersect_tets_real(c_kind* tet_a, c_kind* tet_b, c_kind* tets_c, int* n_tets_c);
void libsupermesh_triangle_area(c_kind* tri, c_kind* area);
void libsupermesh_tetrahedron_volume(c_kind* tet, c_kind* volume);
void libsupermesh_interval_size(c_kind* interval, c_kind* size);
void libsupermesh_intersect_intervals(c_kind* interval_a, c_kind* interval_b, c_kind* intervals_c, int* n_intervals_c);
void libsupermesh_intersect_tet_hex(c_kind* tet_a, c_kind* hex_b, c_kind* tets_c, int* n_tets_c);
void libsupermesh_intersect_tri_quad(c_kind* tri_a, c_kind* quad_b, c_kind* tris_c, int* n_tris_c);

#ifdef __cplusplus
}
#endif

int libsupermesh_max_n_tris_c(int n_lines_b, int n_lines_a){
   int out = n_lines_a*pow(2, n_lines_b) - 2;
   return out;
}

int libsupermesh_max_n_tets_c(int n_planes_b){
   int out = pow(3, n_planes_b);
   return out;
}

