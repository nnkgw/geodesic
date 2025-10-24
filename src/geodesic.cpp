// geodesic.cpp
// Minimal .obj triangle mesh viewer + graph shortest path (Dijkstra on 1-skeleton)
// Libraries: freeglut, OpenGL (fixed pipeline), glm
// Usage: ./app model.obj
// Keys: 1/2 switch method, Space overlay toggle, R re-sample endpoints,
//       +/- zoom, mouse: L-rotate / R-pan, Esc quit
// Mouse wheel zoom is available on FREEGLUT

#if defined(WIN32)
#pragma warning(disable:4996)
#include <GL/freeglut.h>
#elif defined(__APPLE__) || defined(MACOSX)
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cfloat>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <ctime>
#include <cmath>

#define GLM_FORCE_RADIANS
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/euler_angles.hpp>

struct EdgeKey {
  int a,b;
  EdgeKey(int i,int j){ if(i<j){a=i;b=j;} else {a=j;b=i;} }
  bool operator==(const EdgeKey& o) const { return a==o.a && b==o.b; }
};
struct EdgeKeyHash {
  std::size_t operator()(const EdgeKey& k) const noexcept {
    return (std::hash<int>()(k.a)*73856093u) ^ (std::hash<int>()(k.b)*19349663u);
  }
};

struct Mesh {
  std::vector<glm::vec3> V;         // vertices
  std::vector<glm::ivec3> F;        // triangle indices (0-based)
  glm::vec3 bbmin{FLT_MAX}, bbmax{-FLT_MAX};
};

struct Graph { // adjacency list on 1-skeleton
  std::vector<std::vector<std::pair<int,float>>> adj; // (neighbor, length)
};

// Dijkstra header
#include "dijkstra.h"
// New: Selective Refinement stub header (delegates to Dijkstra for now)
#include "selective_refinement.h"

static Mesh gMesh;
static Graph gGraph;

// camera & interaction
static int winW=1280, winH=800;
static float camDist = 3.0f;
static glm::vec3 sceneCenter(0.0f);
static float yaw=0.0f, pitch=0.0f; // radians
static float panX=0.0f, panY=0.0f;
static int lastX=-1, lastY=-1;
static bool lbtn=false, rbtn=false;

// method selection & overlay
enum class Method { Dijkstra, SelectiveRefinement };
static Method gActive = Method::Dijkstra;
static bool   gOverlay = false;

// paths
static int srcIdx=-1, dstIdx=-1;
static std::vector<int> pathVerts;      // kept for backward compatibility (unused for drawing now)
static std::vector<int> pathDij;        // Dijkstra result
static std::vector<int> pathSR;         // Selective Refinement result (currently same as Dijkstra)
static std::mt19937 rng((unsigned)std::time(nullptr));

static void computeBounds(Mesh& M){
  M.bbmin = glm::vec3( FLT_MAX);
  M.bbmax = glm::vec3(-FLT_MAX);
  for(auto& p: M.V){
    M.bbmin = glm::min(M.bbmin, p);
    M.bbmax = glm::max(M.bbmax, p);
  }
}

static bool loadOBJ(const std::string& path, Mesh& M){
  std::ifstream ifs(path);
  if(!ifs) { std::cerr<<"[ERR] cannot open "<<path<<"\n"; return false; }
  std::string line;
  std::vector<glm::vec3> verts;
  std::vector<glm::ivec3> faces;
  while(std::getline(ifs,line)){
    if(line.empty()||line[0]=='#') continue;
    std::istringstream ss(line);
    std::string tag; ss>>tag;
    if(tag=="v"){
      float x,y,z; ss>>x>>y>>z;
      verts.emplace_back(x,y,z);
    } else if(tag=="f"){
      auto readIndex = [&](const std::string& tok)->int{
        std::stringstream s(tok);
        std::string a; std::getline(s,a,'/');
        if(a.empty()) return -1;
        int idx = std::stoi(a);
        if(idx<0) idx = (int)verts.size()+idx+1;
        return idx-1;
      };
      std::string a,b,c; ss>>a>>b>>c;
      if(a.empty()||b.empty()||c.empty()) continue;
      int ia=readIndex(a), ib=readIndex(b), ic=readIndex(c);
      if(ia>=0&&ib>=0&&ic>=0) faces.emplace_back(ia,ib,ic);
    }
  }
  if(verts.empty()||faces.empty()){
    std::cerr<<"[ERR] no verts or faces.\n"; return false;
  }
  M.V.swap(verts); M.F.swap(faces);
  computeBounds(M);
  return true;
}

static float edgeLen(const glm::vec3& a, const glm::vec3& b){
  return glm::length(a-b);
}

static void buildGraph(const Mesh& M, Graph& G){
  G.adj.assign(M.V.size(), {});
  std::unordered_set<EdgeKey,EdgeKeyHash> used;
  used.reserve(M.F.size()*3);
  for(auto &f: M.F){
    int id[3]={f.x,f.y,f.z}; // triangle vtx0, vtx1, vtx2
    for(int e=0;e<3;++e){
      int i=id[e], j=id[(e+1)%3];
      EdgeKey key(i,j);
      if(used.find(key)==used.end()){
        used.insert(key);
        float w=edgeLen(M.V[i], M.V[j]);
        G.adj[i].push_back({j,w});
        G.adj[j].push_back({i,w});
      }
    }
  }
}

// New: draw a path given a vertex index polyline and color/width
static void drawPathWith(const std::vector<int>& path, float r, float g, float b, float width){
  if(path.size()<2) return;
  glDisable(GL_LIGHTING);
  glLineWidth(width);
  glColor3f(r,g,b);
  glBegin(GL_LINE_STRIP);
  for(int vid : path){
    auto &p = gMesh.V[vid];
    glVertex3f(p.x,p.y,p.z);
  }
  glEnd();

  // endpoints as small points
  glPointSize(8.0f);
  glBegin(GL_POINTS);
  glVertex3f(gMesh.V[path.front()].x, gMesh.V[path.front()].y, gMesh.V[path.front()].z);
  glVertex3f(gMesh.V[path.back()].x,  gMesh.V[path.back()].y,  gMesh.V[path.back()].z);
  glEnd();
}

// Legacy helper (kept): draw global pathVerts in red
static void drawPath(){
  drawPathWith(pathVerts, 1.0f,0.1f,0.1f, 4.0f);
}

static void pickRandomEndpoints(){
  std::uniform_int_distribution<int> uni(0,(int)gMesh.V.size()-1);
  int a=uni(rng), b=uni(rng);
  while(b==a) b=uni(rng);
  srcIdx=a; dstIdx=b;
}

static void recomputePath(){
  if(srcIdx<0 || dstIdx<0) pickRandomEndpoints();

  // Baseline Dijkstra
  pathDij.clear();
  bool ok1 = shortestPath(gGraph, srcIdx, dstIdx, pathDij);

  // Selective Refinement (stub -> calls Dijkstra for now)
  pathSR.clear();
  bool ok2 = selectiveRefinementPath(gMesh, gGraph, srcIdx, dstIdx, pathSR, SRParams{});

  // Keep legacy pathVerts to not break older drawPath(); we mirror the active one
  pathVerts = (gActive == Method::Dijkstra ? pathDij : pathSR);

  // If either failed (disconnected graph), re-pick endpoints a few times
  if(!(ok1 && ok2)){
    for(int tries=0; tries<20 && !(ok1 && ok2); ++tries){
      pickRandomEndpoints();
      ok1 = shortestPath(gGraph, srcIdx, dstIdx, pathDij);
      ok2 = selectiveRefinementPath(gMesh, gGraph, srcIdx, dstIdx, pathSR, SRParams{});
      pathVerts = (gActive == Method::Dijkstra ? pathDij : pathSR);
    }
  }
}

static void initCamera(){
  sceneCenter = 0.5f*(gMesh.bbmin+gMesh.bbmax);
  float diag = glm::length(gMesh.bbmax - gMesh.bbmin);
  camDist = (diag>0.0f ? diag*1.8f : 3.0f);
  yaw = 0.4f; pitch = 0.2f;
  panX = panY = 0.0f;
}

static void setProjection(){
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  float aspect = (float)winW/(float)winH;
  gluPerspective(45.0, aspect, 0.001, 1e6);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glm::mat4 I(1.0f);
  glm::mat4 Tpan = glm::translate(I, glm::vec3(panX, panY, 0.0f));
  glm::mat4 R = glm::yawPitchRoll(yaw, pitch, 0.0f);
  glm::vec3 eye = sceneCenter + glm::vec3(R*glm::vec4(0,0,camDist,1));
  glm::mat4 V = glm::lookAt(eye, sceneCenter, glm::vec3(0,1,0));
  glm::mat4 MV = V * Tpan;

  glLoadMatrixf(&MV[0][0]);
}

static void drawMeshWire(){
  glDisable(GL_LIGHTING);
  glColor3f(0.8f,0.8f,0.8f);
  glLineWidth(1.0f);
  glBegin(GL_LINES);
  for(auto &f: gMesh.F){
    int id[3]={f.x,f.y,f.z};
    for(int e=0;e<3;++e){
      glm::vec3 a=gMesh.V[id[e]];
      glm::vec3 b=gMesh.V[id[(e+1)%3]];
      glVertex3f(a.x,a.y,a.z);
      glVertex3f(b.x,b.y,b.z);
    }
  }
  glEnd();
}

static void display(){
  glViewport(0,0,winW,winH);
  glClearColor(0.07f,0.08f,0.10f,1.0f);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  setProjection();
  drawMeshWire();

  // Draw either active method or both (overlay)
  if (gOverlay) {
    // Dijkstra = red, SR = yellow
    drawPathWith(pathDij, 1.0f,0.1f,0.1f, 4.0f);
    drawPathWith(pathSR,  1.0f,0.9f,0.1f, 3.0f);
  } else {
    const std::vector<int>& activePath = (gActive==Method::Dijkstra) ? pathDij : pathSR;
    drawPathWith(activePath, 1.0f,0.1f,0.1f, 4.0f);
  }

  glutSwapBuffers();
}

static void reshape(int w,int h){
  winW=w; winH=h>0?h:1;
  glutPostRedisplay();
}

// ---------- Keyboard ----------
static void keyboard(unsigned char key,int,int){
  switch(key){
    case 27: std::exit(0); break;
    case 'r': case 'R':
      pickRandomEndpoints(); recomputePath(); glutPostRedisplay(); break;
    case '+': camDist *= 0.9f; glutPostRedisplay(); break;
    case '-': camDist *= 1.1f; glutPostRedisplay(); break;

    // New: method switching and overlay toggle
    case '1':
      gActive = Method::Dijkstra;
      // keep legacy global in sync
      pathVerts = pathDij;
      glutPostRedisplay();
      break;
    case '2':
      gActive = Method::SelectiveRefinement;
      pathVerts = pathSR;
      glutPostRedisplay();
      break;
    case ' ':
      gOverlay = !gOverlay;
      glutPostRedisplay();
      break;

    default: break;
  }
}

static void mouse(int button,int state,int x,int y){
  if (button == GLUT_LEFT_BUTTON)  lbtn = (state == GLUT_DOWN);
  if (button == GLUT_RIGHT_BUTTON) rbtn = (state == GLUT_DOWN);
  lastX = x; lastY = y;
}

static void motion(int x,int y){
  int dx = x - lastX;
  int dy = y - lastY;
  lastX = x; lastY = y;

  if (lbtn) { // rotate(yaw/pitch)
    yaw   += dx * 0.005f;
    pitch += dy * 0.005f;
    const float lim = 1.55f;
    if (pitch >  lim) pitch =  lim;
    if (pitch < -lim) pitch = -lim;
  }
  if (rbtn) { // pan
    float s = 0.002f * camDist;
    panX += dx * s;
    panY -= dy * s;
  }
  glutPostRedisplay();
}

#if defined(FREEGLUT)
static void mouse_wheel(int wheel, int direction, int x, int y){
  camDist *= (direction > 0) ? 0.9f : 1.1f;
  if (camDist < 1e-3f) camDist = 1e-3f;
  glutPostRedisplay();
}
#endif

static void usage(){
  std::puts("=== Geodesic (Shortest Path on Mesh) ===");
  std::puts("\nMouse:");
  std::puts("  Left-drag   : rotate camera (yaw/pitch)");
  std::puts("  Right-drag  : pan");
#if defined(FREEGLUT)
  std::puts("  Wheel       : zoom");
#endif

  std::puts("\nKeyboard:");
  std::puts("  1           : select Dijkstra");
  std::puts("  2           : select Selective Refinement (stub)");
  std::puts("  Space       : overlay on/off (draw both)");
  std::puts("  R           : re-sample random endpoints");
  std::puts("  + / -       : zoom in / out");
  std::puts("  ESC         : quit");

  std::puts("\nNotes:");
  std::puts("  - The SR method currently delegates to Dijkstra (stub for future work).");
  std::puts("  - Later, SR will add Steiner points and on-face edges around the path.");
  std::puts("");
}

int main(int argc,char** argv){
  if(argc<2){
    std::cerr<<"Usage: "<<argv[0]<<" model.obj\n";
    return 1;
  }
  if(!loadOBJ(argv[1], gMesh)) return 1;
  buildGraph(gMesh, gGraph);
  initCamera();
  pickRandomEndpoints();
  recomputePath();

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(winW, winH);
  glutCreateWindow("Shortest Path on Mesh (Dijkstra & Selective Refinement)");

  usage();

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
#if defined(FREEGLUT)
  glutMouseWheelFunc(mouse_wheel);
#endif

  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  glutMainLoop();
  return 0;
}
