#include <iostream>
#include <cmath>
#include <triangle.h>

using namespace std;

double boundaryx(double a){

}

double boundaryy(double a){

}

class vec2{
  double x;
  double y;

public:
  vec2 operator+(vec2 other);
  vec2 operator-(void);
  vec2 operator-(vec2 other);
  vec2 operator*(double a);
  vec2 operator/(double a);
  double dot(vec2 other);
  double cross(vec2 other);
  double det(vec2 other);
}

vec2::operator+(vec2 other){
  vec2 ret;
  ret.x = x + other.x;
  ret.y = y + other.y;
  return ret;
}

vec2::operator-(void){
  vec2 ret;
  ret.x = -x;
  ret.y = -y;
  return ret;
}

vec2::operator-(vec2 other){
  vec2 ret;
  ret.x = x - other.x;
  ret.y = y - other.y;
  return ret;
}

vec2::operator*(double a){
  vec2 ret;
  ret.x = a * x;
  ret.y = a * y;
  return ret;
}

vec2::operator/(double a){
  vec2 ret;
  ret.x = x / a;
  ret.y = y / a;
  return ret;
}

vec2::dot(vec2 other){
  return (x*other.x) + (y*other.y);
}

vec2::cross(vec2 other){
  return (y*other.x) + (x*other.y);
}

vec2::det(vec2 other){
  return fabs()(y*other.x) - (x*other.y));
}

class triangle{

public:
  void init(struct triangulateio *in);
  void fill(struct triangulateio *in);
  void freeIfNoti(int *list);
  void freeIfNotd(double *list);
  void free(struct triangulateio *in);
}

triangle::init(struct triangulateio *in){
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

triangle::fill(struct triangulateio *in){
  int i;

  in->numberofpoints = vert;
  in->pointlist = (double *) malloc( vert*2*sizeof(double));
  in->pointmarkerlist = (int *) malloc( vert*sizeof(int));

  for( i=0; i<vert; i++){
    in->pointlist[2*i+0] = boundaryx((double) i/2.0);
    in->pointlist[2*i+1] = boundaryy((double) i/2.0);
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

triangle::freeIfNoti(int *list){
  if (list != NULL){
    free(list);
  }
}

triangle::freeIfNotd(double *list){
  if (list != NULL){
    free(list);
  }
}

triangle::free(struct triangulateio *in){
  freeIfNotd(in->pointlist);
  freeIfNoti(in->pointmarkerlist);
  freeIfNoti(in->segmentlist);
  freeIfNoti(in->segmentmarkerlist);
  freeIfNotd(in->pointattributelist);
  freeIfNoti(in->trianglelist);
  freeIfNotd(in->trianglearealist);
  freeIfNotd(in->triangleattributelist);
  freeIfNotd(in->holelist);
  freeIfNotd(in->regionlist);
  freeIfNoti(in->edgelist);
  freeIfNoti(in->edgemarkerlist);
  freeIfNotd(in->normlist);
}

class integrate{
  int a;
  int b;
  int c;

  vec2 x;
  vec2 y;
  vec2 z;

public:
  void updateTri(int i, int j, int k);
  void passPoints(struct triangulateio *out);
  double Jacobian(void);
  double sumNode(void);
  double sumElements(struct triangulateio *out);
}

integrate::updateTri(i, j, k){
  a = i;
  b = j;
  c = k;
}

integrate::passPoints(struct triangulateio *out){
  x.x = out->pointlist[2 * a];
  x.y = out->pointlist[2 * a + 1];

  y.x = out->pointlist[2 * b];
  y.y = out->pointlist[2 * b + 1];

  z.x = out->pointlist[2 * c];
  z.y = out->pointlist[2 * c + 1];
}

integrate::Jacobian(void){
  vec2 m = y - x;
  vec2 n = z - x;

  return m.det(n);
}

integrate::sumNode(void){
  double weight[4] = {-27.0/48.0, 25.0/48.0, 25.0/48.0, 25.0/48.0};
  double Xi[4] = {1.0/3.0, 0.2, 0.2, 0.6};
  double Eta[4] = {1.0/3.0, 0.2, 0.6, 0.2};

  vec2 m = y-x;
  vec2 n = z-x;

  int i;
  double nodeSum = 0.0;

  for(i = 0; i < 4; i++){
    double xcomp = (m.x * Xi[i]) + (n.x * Eta[i]) + x.x;
    double ycomp = (m.y * Xi[i]) + (n.y * Eta[i]) + x.y;

    nodeSum += exp(-(xcomp*xcomp + ycomp*ycomp)) * weight[i];
  }

  return nodeSum;
}

integrate::sumElements(struct triangulateio *out){
  numElements = out->numberoftriangles;

  int i;
  int j;
  double elementSum = 0.0;

  for(i = 0; i < numElements; i++){
    for(j = 0; j < 3; j++){
      if(j % 3 == 0){
        a = out->trianglelist[i*out->numberofcorners + j];
      } else if(j % 3 == 1){
        b = out->trianglelist[i*out->numberofcorners + j];
      } else {
        c = out->trianglelist[i*out->numberofcorners + j];
      }
    }

    updateTri(a,b,c);
    passPoints(out);

    elementSum += sumNode() * Jacobian();
  }

  return elementSum;
}
