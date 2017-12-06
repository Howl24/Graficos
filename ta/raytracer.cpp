// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
//
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// [/ignore]
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstring>

#if defined __linux__ || defined __APPLE__
    // "Compiled for Linux
#else
    // Windows doesn't define these values by default, Linux does
#define M_PI 3.141592653589793
#define INFINITY 1e8
#endif

float kEpsilon = 1e-8; 

template<typename T>
class Vec3
{
public:
    T x, y, z;
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
    Vec3& normalize()
    {
        T nor2 = length2();
        if (nor2 > 0) {
            T invNor = 1 / sqrt(nor2);
            x *= invNor, y *= invNor, z *= invNor;
        }
        return *this;
    }
    T dotProduct(const Vec3<T> &v) const
    { return x * v.x + y * v.y + z * v.z; }
    Vec3 crossProduct(const Vec3<T> &v) const
    { return Vec3<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); } 
    Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
    Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
    T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
    Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
    Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
    Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
    Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
    T length2() const { return x * x + y * y + z * z; }
    T length() const { return sqrt(length2()); }
    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
    friend Vec3 operator * (const T &r, const Vec3 &v)
    { return Vec3<T>(v.x * r, v.y * r, v.z * r); }
    friend Vec3 operator / (const T &r, const Vec3 &v)
    { return Vec3<T>(r / v.x, r / v.y, r / v.z); } 
};

typedef Vec3<float> Vec3f;

class Triangle
{
public:
    Vec3f center;
    Vec3f v0, v1, v2;
    Vec3f surfaceColor, emissionColor;
    float transparency, reflection;

    Triangle(const Vec3f &center,
             const Vec3f &v0, const Vec3f &v1, const Vec3f &v2,
             const Vec3f &sc,
             const float &refl = 0, const float &transp = 0,
             const Vec3f &ec = 0):

        center(center), v0(v0), v1(v1), v2(v2),
        surfaceColor(sc), emissionColor(ec), 
        transparency(transp), reflection(refl)

    { /* empty */ }

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t) const
    {
         //compute plane's normal
         Vec3f v0v1 = v1 - v0;
         Vec3f v0v2 = v2 - v0;
         // no need to normalize
         Vec3f N = v0v1.crossProduct(v0v2); // N
         float area2 = N.length();
         
         // Step 1: finding P
         
         // check if ray and plane are parallel ?
         float NdotRayDirection = N.dotProduct(dir);
         if (fabs(NdotRayDirection) < kEpsilon) // almost 0
         return false; // they are parallel so they don't intersect !
         
         // compute d parameter using equation 2
         float d = N.dotProduct(v0);
         
         // compute t (equation 3)
         t = (N.dotProduct(orig) + d) / NdotRayDirection;
         // check if the triangle is in behind the ray
         if (t < 0) return false; // the triangle is behind
         
         // compute the intersection point using equation 1
         Vec3f P = orig + t * dir;
         
         // Step 2: inside-outside test
         Vec3f C; // vector perpendicular to triangle's plane
         
         // edge 0
         Vec3f edge0 = v1 - v0;
         Vec3f vp0 = P - v0;
         C = edge0.crossProduct(vp0);
         if (N.dotProduct(C) < 0) return false; // P is on the right side
         
         // edge 1
         Vec3f edge1 = v2 - v1;
         Vec3f vp1 = P - v1;
         C = edge1.crossProduct(vp1);
         if (N.dotProduct(C) < 0) return false; // P is on the right side
         
         // edge 2
         Vec3f edge2 = v0 - v2;
         Vec3f vp2 = P - v2;
         C = edge2.crossProduct(vp2);
         if (N.dotProduct(C) < 0) return false; // P is on the right side;
         
         return true; // this ray hits the triangle 
    }

    void move(const Vec3f offset ){
        v0 += offset;
        v1 += offset;
        v2 += offset;
    }
};

class Sphere
{
public:
    Vec3f center;                           /// position of the triangle
    float radius, radius2;                  /// triangle radius and radius^2
    Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
    float transparency, reflection;         /// surface transparency and reflectivity
    Sphere(
        const Vec3f &c,
        const float &r,
        const Vec3f &sc,
        const float &refl = 0,
        const float &transp = 0,
        const Vec3f &ec = 0) :
        center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
        transparency(transp), reflection(refl)
    { /* empty */ }
    //[comment]
    // Compute a ray-triangle intersection using the geometric solution
    //[/comment]
    bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
    {
        Vec3f l = center - rayorig;
        float tca = l.dot(raydir);
        if (tca < 0) return false;
        float d2 = l.dot(l) - tca * tca;
        if (d2 > radius2) return false;
        float thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;

        return true;
    }
};

std::vector<Triangle> leerOFF(const char* filename){
    FILE* fid = fopen(filename, "rt");

    //Leer formato
    char buffer[1024];
    fscanf(fid, "%s", buffer);

    if(strcmp(buffer, "OFF")!=0) {
        printf("Error de formato\n");
        exit(EXIT_FAILURE);
    }

    int nverts, ntriang, nedges;
    fscanf(fid, "%d %d %d", &nverts, &ntriang, &nedges);
    printf("%d, %d, %d\n", nverts, ntriang, nedges);

    std::vector<Triangle> mesh;
    std::vector<Vec3f> vertex;

    int i;
    float x, y, z;
    for(i = 0; i < nverts; i++){
        fscanf(fid, "%f %f %f",
               &x, &y, &z);
        vertex.push_back(Vec3f(x,y,z));
    }

    for(i = 0; i < ntriang; i++){
        int nv;
        int i1, i2, i3; //indexes
        fscanf(fid, "%d %d %d %d", &nv, &i1, &i2, &i3);

        mesh.push_back(Triangle(Vec3f(5.0, 0, -20),
                                vertex[i1], vertex[i2], vertex[i3], 
                                Vec3f(1.00, 0.32, 0.36),
                                1, 0.0));
    }

    fclose(fid);

    for (i=0;i< mesh.size();i++){
        mesh[i].move(Vec3f(0.0, 0.0, -20.0));
    }



    return mesh;
}

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
#define MAX_RAY_DEPTH 5

float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

Vec3f trace(
    const Vec3f &rayorig,
    const Vec3f &raydir,
    const std::vector<Triangle> &triangles,
    const int &depth)
{
    float tnear = INFINITY;
    const Triangle* triangle = NULL;
    // find intersection of this ray with the triangle in the scene
    for (unsigned i = 0; i < triangles.size(); ++i){
        float t = INFINITY;

        if (triangles[i].intersect(rayorig, raydir, t)){
            if (t < tnear) {
                tnear = t;
                triangle = &triangles[i];
            }
        }
    }

    if (!triangle) return Vec3f(2);
    Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
    Vec3f phit = rayorig + raydir * tnear; // point of intersection
    Vec3f nhit = phit - triangle->center; // normal at the intersection point
    nhit.normalize(); // normalize normal direction
    // If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the triangle so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want
    // positive.
    float bias = 1e-4; // add some bias to the point from which we will be tracing
    bool inside = false;
    if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
    if ((triangle->transparency > 0 || triangle->reflection > 0) && depth < MAX_RAY_DEPTH) {
        float facingratio = -raydir.dot(nhit);
        // change the mix value to tweak the effect
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
        // compute reflection direction (not need to normalize because all vectors
        // are already normalized)
        Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
        refldir.normalize();
        Vec3f reflection = trace(phit + nhit * bias, refldir, triangles, depth + 1);
        Vec3f refraction = 0;
        // if the triangle is also transparent compute refraction ray (transmission)
        if (triangle->transparency) {
            float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
            float cosi = -nhit.dot(raydir);
            float k = 1 - eta * eta * (1 - cosi * cosi);
            Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
            refrdir.normalize();
            refraction = trace(phit - nhit * bias, refrdir, triangles, depth + 1);
        }
        // the result is a mix of reflection and refraction (if the triangle is transparent)
        surfaceColor = (
            reflection * fresneleffect +
            refraction * (1 - fresneleffect) * triangle->transparency) * triangle->surfaceColor;
    }
    else {
        // it's a diffuse object, no need to raytrace any further
        for (unsigned i = 0; i < triangles.size(); ++i) {
            if (triangles[i].emissionColor.x > 0) {
                // this is a light
                Vec3f transmission = 1;
                Vec3f lightDirection = triangles[i].center - phit;
                lightDirection.normalize();
                for (unsigned j = 0; j < triangles.size(); ++j) {
                    if (i != j) {
                        float t;
                        if (triangles[j].intersect(phit + nhit * bias, lightDirection, t)) {
                            transmission = 0;
                            break;
                        }
                    }
                }
                surfaceColor += triangle->surfaceColor * transmission *
                std::max(float(0), nhit.dot(lightDirection)) * triangles[i].emissionColor;
            }
        }
    }

    return surfaceColor + triangle->emissionColor;
}

void render(const std::vector<Triangle> &triangles)
{
    unsigned width = 640, height = 480;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);

    for (unsigned y = 0; y < height; ++y){
        for (unsigned x = 0 ; x< width; ++x, ++ pixel){
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -1);
            raydir.normalize();
            *pixel = trace(Vec3f(0), raydir, triangles, 0);
        }
    }

    std::ofstream ofs("./cubo.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width*height; ++ i){
        ofs << (unsigned char)(std::min(float(1), image[i].x)*255) <<
               (unsigned char)(std::min(float(1), image[i].y)*255) <<
               (unsigned char)(std::min(float(1), image[i].z)*255);
    }
    ofs.close();
    delete[] image;
}

//[comment]
// In the main function, we will create the scene which is composed of 5 triangles
// and 1 light (which is also a triangle). Then, once the scene description is complete
// we render that scene, by calling the render() function.
//[/comment]
int main(int argc, char **argv)
{
    printf("Hello\n");
    std::vector<Triangle> triangles;
    triangles = leerOFF("avion.off");
    printf("triangles loaded\n");
    render(triangles);

    return 0;
}
