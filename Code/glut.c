/***********************************************************
 * A Template for building OpenGL applications using GLUT
 *
 * Author: Perspective @ cprogramming.com
 * Date  : Jan, 2005
 *
 * Description: 
 *   This code initializes an OpenGL ready window
 * using GLUT.  Some of the most common callbacks
 * are registered to empty or minimal functions.
 *
 * This code is intended to be a quick starting point
 * when building GLUT applications.
 *
 ***********************************************************/
 
//#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "GL/openglut.h"
//#include <GLUT/glut.h>
//#include <OpenGL/gl.h>
#include <vector>
#include "matrix_t.h"

void stiffnessMatrix (double x1, double y1, double x2, double y2, double ea , double k[4][4])
{
  double t = atan2(y2-y1,x2-x1);
  double c = cos(t);
  double s = sin(t);
  double l = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  double f = ea/l;
  k[0][0] = k[2][2] = c*c*f;
  k[1][1] = k[3][3] = s*s*f;
  k[0][2] = k[2][0] = -c*c*f;
  k[3][1] = k[1][3] = -s*s*f;
  k[0][1] = k[1][0] = k[3][2] = k[2][3] = s*c*f;
  k[3][0] = k[0][3] = k[2][1] = k[1][2] = -s*c*f;
}

void postPro (int B, double *xy, int *ind, 
	     double *ea, double *x, double *n) {
  double stiff [4][4];
  for (int i=0;i<B;i++){
    int n1 = ind[2*i];   
    int n2 = ind[2*i+1];
    int indx[4] = {2*n1,2*n1+1,2*n2,2*n2+1};
    stiffnessMatrix ( xy[indx[0]], xy[indx[1]], 
                      xy[indx[2]], xy[indx[3]], ea[i] , stiff );
    double f [4] = {0,0,0,0};
    for (int j=0;j<4;j++) {
      for (int k=0;k<4;k++) {
	f[j] += stiff[j][k] * x[indx[k]];  
      }
    }
    double t = atan2( xy[indx[1]]- xy[indx[3]], xy[indx[0]]- xy[indx[2]]);
    n[i] = f[0] * cos(t) + f[1] * sin(t);
    double l2 = sqrt((xy[indx[1]]+x[indx[1]]- xy[indx[3]]-x[indx[3]])*
		     (xy[indx[1]]+x[indx[1]]- xy[indx[3]]-x[indx[3]])+
		     (xy[indx[0]]+x[indx[0]]- xy[indx[2]]-x[indx[2]])*
		     (xy[indx[0]]+x[indx[0]]- xy[indx[2]]-x[indx[2]]));
    double l1 = sqrt((xy[indx[1]]- xy[indx[3]])*(xy[indx[1]]- xy[indx[3]])+
		     (xy[indx[0]]- xy[indx[2]])*(xy[indx[0]]- xy[indx[2]]));
    double dl = l2 - l1;
    double eps = dl / l1;
    double N = eps*ea[i];
    
    printf("n[%d] = %g N = %g\n",i,n[i],N);
  }
}

int trussSolver (int N, int B, int M, double *xy, int *ind, 
                 double *f, int *nc, double *vc, double *ea,  
                 double *x, double *r) {
  double stiff [4][4];
  matrix_t *K = create_matrix (2*N+M, 2*N+M);
  double *F = (double*) malloc ((2*N+M)*sizeof(double));
  double *X = (double*) malloc ((2*N+M)*sizeof(double));
  for (int i=0;i<B;i++){
    int n1 = ind[2*i];   
    int n2 = ind[2*i+1];
    int indx[4] = {2*n1,2*n1+1,2*n2,2*n2+1};
    stiffnessMatrix ( xy[indx[0]], xy[indx[1]], 
                      xy[indx[2]], xy[indx[3]], ea[i] , stiff );
    for (int j=0;j<4;j++) for (int k=0;k<4;k++)
      add_to_matrix(K,indx[j],indx[k],stiff[j][k]);
  }
  for (int i=0;i<2*N;i++) F[i] = f[i];
  for (int i=0;i<M;i++) {
    int indx[2] = {2*nc[i], 2*nc[i] +1};
    double v[3] = {vc[3*i], vc[3*i+1], vc[3*i+2]};
    double norm = sqrt(v[0]*v[0] + v[1]*v[1]);
    add_to_matrix(K,2*N+i,indx[0], v[0]/norm);
    add_to_matrix(K,2*N+i,indx[1], v[1]/norm);
    add_to_matrix(K,indx[0],2*N+i, v[0]/norm);
    add_to_matrix(K,indx[1],2*N+i, v[1]/norm);
    F[2*N+i] = v[2];
  }
  bool result = solve_linear_system (K, F, X);
  if (result == false){
    printf("ERROR : the truss is not stable\n");
  }
  for (int i=0;i<2*N;i++) x[i] = X[i];
  for (int i=0;i<M;i++) {r[i] = X[2*N+i];printf("r(%d) = %g\n",i,r[i]);}

  delete_matrix(K);
  free(F);
  free(X);
  return result;
}

void draw() ;

struct frame_s { 
  double _min[2],_max[2];
  std::vector<double> _pos;
  std::vector<double> _dx;
  std::vector<double> _f;
  std::vector<double> _ea;
  std::vector<double> _r;
  std::vector<double> _vc;
  std::vector<int> _topo;
  std::vector<int> _nc;

  frame_s () {
    _min[0] = _min[1] = -1;  
    _max[0] = _max[1] = 1;  
  }

  bool solve () {
    _dx = _pos;
    trussSolver (_pos.size()/2, 
		 _topo.size()/2, 
		 _nc.size()/3, 
		 _pos.data(), 
		 _topo.data(), 
                 _f.data(),
		 _nc.data(),
		 _vc.data(), 
		 _ea.data(),
		 _dx.data(),
		 _r.data());    
  }

  void addJoint (double x, double y) {
    _min[0] = std::min(x,_min[0]);
    _min[1] = std::min(y,_min[1]);
    _max[0] = std::max(x,_max[0]);
    _max[1] = std::max(y,_max[1]);    
    _pos.push_back(x);
    _pos.push_back(y);
  }
  void addBar (int a, int b) {
    if (a != b){
      for (int i=0;i<_topo.size();i+=2){
	int i1 = _topo[i];
	int i2 = _topo[i+1];
	if (i1 == a && i2 == b)return;
	if (i1 == b && i2 == a)return;
      }
      _topo.push_back(a);
      _topo.push_back(b);
    }
  }
};

static frame_s FRAME;
static int START_POINT;
static double XPOS,YPOS;

/* process menu option 'op' */
void menu(int op) {
 
  switch(op) {
  case 'Q':
  case 'q':
    exit(0);
  }
}
 
void getCoords (double &xmin, double &xmax, double &ymin, double &ymax){
  double sx = glutGet(GLUT_WINDOW_WIDTH);
  double sy = glutGet(GLUT_WINDOW_HEIGHT);
  
  const double MAX = std::max( FRAME._max[0]-FRAME._min[0] , 
			       FRAME._max[1]-FRAME._min[1]);
  const double MIDX = 0.5*(FRAME._max[0]+FRAME._min[0]);
  const double MIDY = 0.5*(FRAME._max[1]+FRAME._min[1]);
  if (sx > sy){
    xmin=MIDX-sx*.7*MAX/sy;
    xmax=MIDX+sx*.7*MAX/sy;
    ymin=MIDY-.7*MAX;
    ymax=MIDY+.7*MAX;
  }
  else{
    xmin=MIDX-.7*MAX;
    xmax=MIDX+.7*MAX;
    ymin=MIDY-sy*.7*MAX/sx;
    ymax=MIDY+sy*.7*MAX/sx;
  }
}
/* executed when a regular key is pressed */
void keyboardDown(unsigned char key, int x, int y) {
 
  switch(key) {
  case 'p':
    {
      double xmin, xmax, ymin, ymax;
      getCoords (xmin, xmax, ymin, ymax);
      double sx = glutGet(GLUT_WINDOW_WIDTH);
      double sy = glutGet(GLUT_WINDOW_HEIGHT);
      double xi  = 1+(x - sx) / sx;
      double eta = -(y - sy) / sy;
      printf("xi = %g eta = %g\n",xi,eta);
      FRAME.addJoint(xmin * (1.-xi) + xmax * xi,ymin * (1.-eta) + ymax * eta);
      draw();
    }
    break;
  case 'b':
    break;
  case 'Q':
  case 'q':
  case  27:   // ESC
    exit(0);
  }
}
 
/* executed when a regular key is released */
void keyboardUp(unsigned char key, int x, int y) {
 
}
 
/* executed when a special key is pressed */
void keyboardSpecialDown(int k, int x, int y) {
 
}
 
/* executed when a special key is released */
void keyboardSpecialUp(int k, int x, int y) {
 
}
 
/* reshaped window */
void reshape(int width, int height) {
 
  GLfloat fieldOfView = 90.0f;
  glViewport (0, 0, (GLsizei) width, (GLsizei) height);
 
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(fieldOfView, (GLfloat) width/(GLfloat) height, 0.1, 500.0);
 
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}
 
/* executed when button 'button' is put into state 'state' at screen position ('x', 'y') */
void mouseClick(int button, int state, int x, int y) {
  if (FRAME._pos.size() < 4)return;
  double xmin, xmax, ymin, ymax;
  getCoords (xmin, xmax, ymin, ymax);
  double sx = glutGet(GLUT_WINDOW_WIDTH);
  double sy = glutGet(GLUT_WINDOW_HEIGHT);
  double xi  = 1+(x - sx) / sx;
  double eta = -(y - sy) / sy;
  double X = xmin * (1.-xi) + xmax * xi;
  double Y = ymin * (1.-eta) + ymax * eta;
  double closestDistance = 1.e22;
  if (state == 0)  {
    for (int i=0;i<FRAME._pos.size();i+=2){
      double xx = FRAME._pos[i];
      double yy = FRAME._pos[i+1];
      double d = sqrt ((X-xx)*(X-xx)+(Y-yy)*(Y-yy));
      if (d < closestDistance){
	closestDistance = d;
	START_POINT = i/2;
      }
    }
    XPOS=X;
    YPOS=Y;
  }
  else if (state == 1)  {
    int END_POINT;
    for (int i=0;i<FRAME._pos.size();i+=2){
      double xx = FRAME._pos[i];
      double yy = FRAME._pos[i+1];
      double d = sqrt ((X-xx)*(X-xx)+(Y-yy)*(Y-yy));
      if (d < closestDistance){
	closestDistance = d;
	END_POINT = i/2;
      }
    }
    FRAME.addBar (START_POINT,END_POINT);
    START_POINT = -1;
  }
  draw();
}
 
/* executed when the mouse moves to position ('x', 'y') */
void mouseMotion(int x, int y) {
  if (START_POINT != -1){
    double xmin, xmax, ymin, ymax;
    getCoords (xmin, xmax, ymin, ymax);
    double sx = glutGet(GLUT_WINDOW_WIDTH);
    double sy = glutGet(GLUT_WINDOW_HEIGHT);
    double xi  = 1+(x - sx) / sx;
    double eta = -(y - sy) / sy;
    XPOS = xmin * (1.-xi) + xmax * xi;
    YPOS = ymin * (1.-eta) + ymax * eta;
    draw();
  }
}
 
/* render the scene */
void draw() {
 
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity();  
  
  double sx = glutGet(GLUT_WINDOW_WIDTH);
  double sy = glutGet(GLUT_WINDOW_HEIGHT);

  const double MAX = std::max( FRAME._max[0]-FRAME._min[0] , 
			       FRAME._max[1]-FRAME._min[1]);
  const double MIDX = 0.5*(FRAME._max[0]+FRAME._min[0]);
  const double MIDY = 0.5*(FRAME._max[1]+FRAME._min[1]);

  if (sx > sy){
    glOrtho (MIDX-sx*MAX*.7/sy,MIDX+sx*MAX*.7/sy,MIDY-MAX*.7,MIDY+MAX*.7,-1,1);
  }
  else{
    glOrtho (MIDX-MAX*.7,MIDX+MAX*.7,MIDY-sy*MAX*.7/sx,MIDY+sy*MAX*.7/sx,-1,1);
  }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
 
  /* render the scene here */

  glBegin(GL_LINES);
  glColor3d(1.,0.,0.);
  for (unsigned int i=0;i<FRAME._topo.size();i+=2){
    int i1 = FRAME._topo[i];
    int i2 = FRAME._topo[i+1];
    double x1 = FRAME._pos[2*i1];
    double y1 = FRAME._pos[2*i1+1];
    double x2 = FRAME._pos[2*i2];
    double y2 = FRAME._pos[2*i2+1];
    glVertex2d(x1,y1);
    glVertex2d(x2,y2);
  }
  glEnd();


  glPointSize(10.0f);
  glColor3d(0.,1.,0.);
  glEnable( GL_POINT_SMOOTH );
  glEnable( GL_BLEND );
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
  glBegin(GL_POINTS);
  for (unsigned int i=0;i<FRAME._pos.size();i+=2){
    double x1 = FRAME._pos[i];
    double y1 = FRAME._pos[i+1];
    glVertex2d(x1,y1);
  }
  glEnd();

  if (START_POINT != -1){
    glBegin(GL_LINES);
    glColor3d(0.,0.,1.);
    double x1 = FRAME._pos[2*START_POINT];
    double y1 = FRAME._pos[2*START_POINT+1];
    glVertex2d(x1,y1);
    glVertex2d(XPOS,YPOS);
    glEnd();
  }


  glFlush();
  glutSwapBuffers();
}
 
/* executed when program is idle */
void idle() { 
 
}
 
/* initialize OpenGL settings */
void initGL(int width, int height) {
 
  reshape(width, height);
 
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glClearDepth(1.0f);
 
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
}
 
/* initialize GLUT settings, register callbacks, enter main loop */
int main(int argc, char** argv) {   
  START_POINT = -1;
  
  glutInit(&argc, argv);
 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(800, 600);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("Truss Solver");
 
  // register glut call backs
  glutKeyboardFunc(keyboardDown);
  glutKeyboardUpFunc(keyboardUp);
  glutSpecialFunc(keyboardSpecialDown);
  glutSpecialUpFunc(keyboardSpecialUp);
  glutMouseFunc(mouseClick);
  glutMotionFunc(mouseMotion);
  glutReshapeFunc(reshape);
  glutDisplayFunc(draw);  
  glutIdleFunc(idle);
  glutIgnoreKeyRepeat(true); // ignore keys held down
 
  // create a sub menu 
  int subMenu = glutCreateMenu(menu);
  glutAddMenuEntry("Do nothing", 0);
  glutAddMenuEntry("Really Quit", 'q');
 
  // create main "right click" menu
  glutCreateMenu(menu);
  glutAddSubMenu("Sub Menu", subMenu);
  glutAddMenuEntry("Quit", 'q');
  glutAttachMenu(GLUT_RIGHT_BUTTON);
 
  initGL(800, 600);
 
  glutMainLoop();
  return 0;
}
