#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "triangle.h"

#define TWOPI 6.28318530718

//Command Prompt switches:

// -p Triangulates a Planar Straight Line Graph (.poly file).
// -r Refines a previously generated mesh.
// -q Quality mesh generation with no angles smaller than 20 degrees. An
//      alternate minimum angle may be specified after the `q'.
// -a Imposes a maximum triangle area constraint. A fixed area constraint (that
//      applies to every triangle) may be specified after the `a', or varying area constraints may be read from a .poly file or .area file.
// -u Imposes a user-defined constraint on triangle size.
// -A Assigns a regional attribute to each triangle that identifies what
//      segment-bounded region it belongs to.
// -c Encloses the convex hull with segments.
// -D Conforming Delaunay: use this switch if you want all triangles in the
//      mesh to be Delaunay, and not just constrained Delaunay; or if you want
//      to ensure that all Voronoi vertices lie within the triangulation.
// -j Jettisons vertices that are not part of the final triangulation from the
//      output .node file (including duplicate input vertices and vertices
//      ``eaten'' by holes).
// -e Outputs (to an .edge file) a list of edges of the triangulation.
// -v Outputs the Voronoi diagram associated with the triangulation. Does not
//      attempt to detect degeneracies, so some Voronoi vertices may be
//      duplicated.
// -n Outputs (to a .neigh file) a list of triangles neighboring each triangle.
// -g Outputs the mesh to an Object File Format (.off) file, suitable for
//      viewing with the Geometry Center's Geomview package.
// -B Suppresses boundary markers in the output .node, .poly, and .edge output
//      files.
// -P Suppresses the output .poly file. Saves disk space, but you lose the
//      ability to maintain constraining segments on later refinements of the
//      mesh.
// -N Suppresses the output .node file.
// -E Suppresses the output .ele file.
// -I Suppresses mesh iteration numbers.
// -O Suppresses holes: ignores the holes in the .poly file.
// -X Suppresses exact arithmetic.
// -z Numbers all items starting from zero (rather than one). Note that this
//      switch is normally overrided by the value used to number the first
//      vertex of the input .node or .poly file. However, this switch is useful
//      when calling Triangle from another program.
// -o2 Generates second-order subparametric elements with six nodes each.
// -Y Prohibits the insertion of Steiner points on the mesh boundary. If
//      specified twice (-YY), it prohibits the insertion of Steiner points on
//      any segment, including internal segments.
// -S Specifies the maximum number of added Steiner points.
// -i Uses the incremental algorithm for Delaunay triangulation, rather than
//      the divide-and-conquer algorithm.
// -F Uses Steven Fortune's sweepline algorithm for Delaunay triangulation,
//      rather than the divide-and-conquer algorithm.
// -l Uses only vertical cuts in the divide-and-conquer algorithm. By default,
//      Triangle uses alternating vertical and horizontal cuts, which usually
//      improve the speed except with vertex sets that are small or short and
//      wide. This switch is primarily of theoretical interest.
// -s Specifies that segments should be forced into the triangulation by
//      recursively splitting them at their midpoints, rather than by
//      generating a constrained Delaunay triangulation. Segment splitting is
//      true to Ruppert's original algorithm, but can create needlessly small
//      triangles. This switch is primarily of theoretical interest.
// -C Check the consistency of the final mesh. Uses exact arithmetic for
//      checking, even if the -X switch is used. Useful if you suspect Triangle
//      is buggy.
// -Q Quiet: Suppresses all explanation of what Triangle is doing, unless an
//      error occurs.
// -V Verbose: Gives detailed information about what Triangle is doing. Add
//      more `V's for increasing amount of detail. `-V' gives information on
//    algorithmic progress and detailed statistics.
// -h Help: Displays complete instructions

//Structure for a 2-D vector

typedef struct{
  double x;
  double y;
} vec2;

//Fill points on a circle

double circle(int i, int vertices, double r){
  if (i%2 == 0){
    return (r * cos((double)(i/2) / vertices * 2 * 3.1415926));
  } else{
    return (r * sin((double)(i/2) / vertices * 2 * 3.1415926));
  }
}

//Initiate the arrays that will be used for input into triangles and that then
//  can be read afterwards

void init_triangleio( struct triangulateio *in){
  in->numberofpoints = 0;
  in->numberofpointattributes = 0;
  in->pointlist = (double *) NULL;
  in->pointmarkerlist = (int *) NULL;
  in->pointattributelist = (double *) NULL;

  in->numberofsegments = 0;
  in->segmentlist = (int *) NULL;
  in->segmentmarkerlist = (int *) NULL;

  in->numberoftriangles = 0;
  in->numberofcorners = 0;
  in->numberoftriangleattributes = 0;

  in->trianglelist = (int *) NULL;
  in->triangleattributelist = (double *) NULL;
  in->trianglearealist = (double *) NULL;

  in->numberofholes = 0;
  in->holelist = (double *) NULL;

  in->numberofregions = 0;
  in->regionlist = (double *) NULL;

  in->numberofedges = 0;
  in->edgelist = (int *) NULL;
  in->edgemarkerlist = (int *) NULL;
  in->normlist = (double *) NULL;

}

//Fill  the arrays that will go into triangle

void fill_circle( struct triangulateio *in, int vert, double r){
  int i;

  in->numberofpoints = vert;
  in->pointlist = (double *) malloc( vert*2*sizeof(double));
  in->pointmarkerlist = (int *) malloc( vert*sizeof(int));

  for( i=0; i<vert; i++){
    in->pointlist[2*i+0] = r*cos( ((double) i)/ vert * TWOPI);
    in->pointlist[2*i+1] = r*sin( ((double) i)/ vert * TWOPI);
    in->pointmarkerlist[i] = 1;
  }

  in->numberofsegments = vert;
  in->segmentlist = (int *) malloc( vert*2*sizeof(int));
  in->segmentmarkerlist = (int *) malloc( vert*sizeof(int));

  in->segmentlist[2*0+0] = vert-1;
  in->segmentlist[2*0+1] = 0;
  in->segmentmarkerlist[0] = 1;
  for( i=1; i<vert; i++){
    in->segmentlist[2*i+0] = i-1;
    in->segmentlist[2*i+1] = i;
    in->segmentmarkerlist[i] = 1;
  }
}

//Free the arryas that contain integrals if not free already

void free_ifnoti(int * list){
  if(list != NULL){
    free(list);
  }
}

//Free the arrays that contain double if not freed already

void free_ifnotd(double * list){
  if(list != NULL){
    free(list);
  }
}

//Actually free the arrays

void free_circle(struct triangulateio *in){
  free_ifnotd(in->pointlist);
  free_ifnoti(in->pointmarkerlist);
  free_ifnoti(in->segmentlist);
  free_ifnoti(in->segmentmarkerlist);
  free_ifnotd(in->pointattributelist);
  free_ifnoti(in->trianglelist);
  free_ifnotd(in->trianglearealist);
  free_ifnotd(in->triangleattributelist);
  free_ifnotd(in->holelist);
  free_ifnotd(in->regionlist);
  free_ifnoti(in->edgelist);
  free_ifnoti(in->edgemarkerlist);
  free_ifnotd(in->normlist);
}

//Calculate the jacobian. abs() did not work, I believe this is because it was
//too small.

double jacobian(vec2 a, vec2 b, vec2 c){
  double m = b.x - a.x;
  double n = b.y - a.y;
  double o = c.x - a.x;
  double p = c.y - a.y;
  double ret = ((m*o) - (n*p));
  if(ret < 0){
    ret = -ret;
  }

  printf("Jacobian: %f\n", ret);

  return ret;
}

//Array structure
//Not used unless Transformation() is being used

typedef struct{
  vec2 a;
  vec2 b;
} arr22;

//Transformation from x to eta space
//Currently not being used, This addition changed the value from ~451 to ~158
//Actual: 56 for hemisphere

arr22 Transformation(vec2 a, vec2 b, vec2 c){
  arr22 ret;
  double xab = b.x - a.x;
  double xac = c.x - a.x;
  double yab = b.y - a.y;
  double yac = c.y - a.y;

  double constant = 1.0/((xab * yac) - (xac * yab));

  ret.a.x = yac * constant;
  ret.b.x = -xac * constant;
  ret.a.y = -yab * constant;
  ret.b.y = xab * constant;

  return ret;
}

//Sum the function times the weight for each node in an element

double gaussianIntSingle(vec2 a, vec2 b, vec2 c, double r){
  double xab = b.x - a.x;
  double xac = c.x - a.x;
  double yab = b.y - a.y;
  double yac = c.y - a.y;

  double node = 0.0;
  double gaussian;
  double hemisphere;

  double weight[4] = {-27.0/48.0, 25.0/48.0, 25.0/48.0, 25.0/48.0};
  double Xi[4] = {1.0/3.0, 0.2, 0.2, 0.6};
  double Eta[4] = {1.0/3.0, 0.2, 0.6, 0.2};

  for(int i = 0; i < 4; i++){
    double x = (xab * Xi[i]) + (xac * Eta[i]) + a.x;
    double y = (yab * Xi[i]) + (yac * Eta[i]) + a.y;

    gaussian = exp(-(pow(x,2.0) + pow(y,2.0))) * weight[i];
    hemisphere = r - sqrt(pow(x,2) + pow(y,2)) * weight[i];

    node += hemisphere;
  }

  return node * jacobian(a,b,c);
}

//Structure for which nodes are in a triangle

typedef struct{
  int a;
  int b;
  int c;
} tripoint;

//Finding the x and y location for each node of a triangle and then finding the
//  value for that element that will eventually be summed.

double passPoints(tripoint triangle, double * pointlist, double r){
  vec2 a = {pointlist[2 * triangle.a], pointlist[2 * triangle.a + 1]};
  vec2 b = {pointlist[2 * triangle.b], pointlist[2 * triangle.b + 1]};
  vec2 c = {pointlist[2 * triangle.c], pointlist[2 * triangle.c + 1]};

  double element = gaussianIntSingle(a,b,c,r);
  return element;
}

int main(){
//Standard start to triangle

  struct triangulateio in;
  struct triangulateio out;
  struct triangulateio *vorout = (struct triangulateio *) NULL;
  char Switches[] = "pq30a0.1z";

//Initialize the the arrays that will be used to pass data to triangle and to
//take data from triangle

  init_triangleio( &in);
  init_triangleio( &out);

//fill the circle for points

  fill_circle( &in, 50, 3.0);

//Execute triangulate

  triangulate( Switches, &in, &out, vorout);

//Initialize the variable for the the final summation

  tripoint triangle;
  double element;

//Print he triangle, the nodes associated with that triangle and the moving sum
//  of the elements

  printf("Number of Triangles: %i\n", out.numberoftriangles);
  for(int i = 0; i < out.numberoftriangles; i++){
    printf("Triangle %i: ",i);

    for(int j = 0; j < out.numberofcorners; j++){
      printf("\t %d", out.trianglelist[i*out.numberofcorners + j]);

      if(j % 3 == 0){
        triangle.a = out.trianglelist[i*out.numberofcorners + j];
      } else if(j % 3 == 1){
        triangle.b = out.trianglelist[i*out.numberofcorners + j];
      } else {
        triangle.c = out.trianglelist[i*out.numberofcorners + j];
      }
    }

    printf("\n");
    printf("\t %d %d %d \n", triangle.a, triangle.b, triangle.c);
    element += passPoints(triangle, out.pointlist, 3.0);
    printf("Element: %f", element);
    printf("\n");
  }

  printf("The integration over this mesh is %f \n", element);

//Freeing the arrays

  free_circle( &in);
  free_circle( &out);
//  trifree( &out);

  return 0;
}
