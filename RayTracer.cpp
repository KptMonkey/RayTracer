#include <iostream>
#include <fstream>
#include <string>
#include <glm/glm.hpp>
#include <random>
#include <chrono>
#include <limits>
struct HitRecord{
  float t;
  glm::vec3 p;
  glm::vec3 normal;
};


struct Ray{
  Ray( glm::vec3 const a, glm::vec3 const b) :
    Origin(a), Direction(b)
  {}
  glm::vec3
  getPointAtParameter(float t){
    return Origin + t*Direction;
  }
  glm::vec3 Origin, Direction;
};

class Hitable{
public:
  virtual bool
  hit( Ray & r, float t_min, float t_max, HitRecord & ref) = 0;
};

class Sphere : public Hitable {
public:
  Sphere(){}
  Sphere(glm::vec3 const & center, float radius)
  : m_Center(center), m_Radius(radius)
  {}
  virtual bool
  hit( Ray & r, float t_min, float t_max, HitRecord & ref);

  glm::vec3 m_Center;
  float m_Radius;
};

bool
Sphere::hit( Ray &r, float t_min, float t_max, HitRecord &ref)
{
  auto cc = r.Origin - m_Center;
  //Reminder ax² + 2bx +c = 0 Mitternachtsformel^^
  float a = glm::dot(r.Direction, r.Direction);
  float b = glm::dot(cc,r.Direction);
  float c = glm::dot(cc,cc) - m_Radius*m_Radius;
  float diskriminante = b*b - a* c;
  if (diskriminante > 0.0f){
    float temp = (-b + sqrt(diskriminante)/ a);
    if(temp > t_min && temp < t_max){
      ref.t = temp;
      ref.p = r.getPointAtParameter(ref.t);
      ref.normal = glm::normalize(ref.p - m_Center);
      return true;
    }
    temp = (-b - sqrt(diskriminante)/ a);
    if(temp > t_min && temp < t_max){
      ref.t = temp;
      ref.p = r.getPointAtParameter(temp);
      ref.normal = glm::normalize(ref.p - m_Center);
      return true;
    }
   }
   return false;
}

class HitableList : public Hitable{
public:
  HitableList(){}
  HitableList( Hitable **l, int n)
    : m_List(l), m_Size(n)
  {}
  virtual bool
  hit( Ray & r, float t_min, float t_max, HitRecord & ref);
  Hitable ** m_List;
  int m_Size;
};

bool
HitableList::hit(Ray &r, float t_min, float t_max, HitRecord &ref){

HitRecord tempRec;
bool hitAnything = false;
float closest = t_max;
for ( int i = 0; i<m_Size; ++i) {
  if(m_List[i]->hit( r, t_min, closest, tempRec)) {
    hitAnything = true;
    closest = tempRec.t;
    ref = tempRec;
  }
 }
 return hitAnything;
}

glm::vec3 randomInUnitSphere(){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    glm::vec3 p(2.0f);
    while( (glm::length(p)*glm::length(p)) >= 1.0f ){
        p = 2.0f * glm::vec3(dist(mt),dist(mt),dist(mt)) - glm::vec3(1.0f);
    }
    return p;
}

glm::vec3
colorFunc(Ray & ray, Hitable * world){
  HitRecord hitRec;
  if( world != nullptr && world->hit(ray, 0.001f,std::numeric_limits<float>::max(), hitRec )) {
     auto target = hitRec.p + hitRec.normal + randomInUnitSphere();
     auto bla = Ray(hitRec.p, target-hitRec.p);
     return 0.5f * colorFunc( bla,world);
  }
  else{
    glm::vec3 normRay = glm::normalize(ray.Direction);
    auto t = 0.5f * (normRay.y + 1.0f);
    return (1.0f-t)*glm::vec3(1.0f) + t * glm::vec3(0.5f, 0.7f, 1.0f);
  }
}

class Camera {
public:
  Camera(glm::vec3 blc, glm::vec3 hor, glm::vec3 ver, glm::vec3 orig)
    : m_BlCorner(blc),
      m_Horizontal(hor),
      m_Vertical(ver),
      m_Origin(orig)
  {}
  Ray getRay( float u, float v) {
    return Ray(m_Origin, m_BlCorner + u*m_Horizontal+v*m_Vertical);
  }
  glm::vec3 m_BlCorner;
  glm::vec3 m_Horizontal;
  glm::vec3 m_Vertical;
  glm::vec3 m_Origin;

};


int main() {
 std::string fileName = "rayImage.ppm";
 std::ofstream out(fileName, std::ios::out);
float x = 200.0f;
float y = 100.0f;
float s = 50.0f;
glm::vec3 blCorner(-2.0f, -1.0f, -1.0f);
glm::vec3 horizontal(4.0f, 0.0f, 0.0f);
glm::vec3 vertical(0.0f,2.0f,0.0f);
glm::vec3 origin(0.0f,0.0f,0.0f);
Camera cam(blCorner, horizontal, vertical, origin);
out << "P3\n" << x << " " << y << "\n255\n";
Hitable *list[2];
glm::vec3 center1(0.0f, 0.0f, -1.0f);
glm::vec3 center2(0.0f, -100.5f, -1.0f);
list[0] = new Sphere(center1,0.5f);
list[1] = new Sphere(center2,100.0f);
Hitable * world = new HitableList(list, 2);

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<float> dist(0.0f, 1.0f);
for ( int j = y-1; j >= 0; --j){
  std::cout << "Status: " << j << "\n";
  for (int i = 0; i < x; ++i){
    glm::vec3 col(0.0f, 0.0f, 0.0f);
    for ( int k = 0; k < s; k++) {
//      std::cout << "RandomNumber: " << dist(mt) << "\n";
      float u = static_cast<float>((i+dist(mt))/x);
      float v = static_cast<float>((j+dist(mt))/y);
      auto ray = cam.getRay(u,v);
      col += colorFunc(ray,world);
    }
    col = col/s;
    col = glm::vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
    int ir = static_cast<int>(255.99f * col[0]);
    int ig = static_cast<int>(255.99f * col[1]);
    int ib = static_cast<int>(255.99f * col[2]);

    out << ir << " " << ig << " " << ib << "\n";
  }
}
 out.close();

}
