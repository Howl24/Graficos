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
};




typedef Vec3<float> Vec3f;

class Sphere
{
public:
    Vec3f center;                           /// position of the sphere
    float radius, radius2;                  /// sphere radius and radius^2
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
    // Compute a ray-sphere intersection using the geometric solution
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


//Estructuras para objeto OFF
typedef struct Vertex{
    float x, y, z;
}Vertex;

typedef struct Triangle{
    unsigned int indices[3];
}Triangle;

typedef struct Mesh{
    //Informacion de estructura
    int numVertices;
    int numTriangles;
    Vertex* vertices;
    Triangle* triangles;

    //Información para transformación inicial
    Vertex center;
    float scale;

    bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const{
        return true;
    }

    //Matriz de transformación
    //glm::mat4 model_transform;

    //Buffers para graficado
    //GLfloat* object_vertices;
    //GLushort* object_indexes;

    //Id's para buffers
    //GLuint vbo_object;
    //GLuint ibo_object;
}Mesh;

Mesh leerOFF(const char* filename){
    FILE* fid = fopen(filename, "rt");

    //Leer formato
    char buffer[1024];
    fscanf(fid, "%s", buffer);

    if(strcmp(buffer, "OFF")!=0){
        printf("Error de formato\n");
        exit(EXIT_FAILURE);
    }

    int nverts, ntriang, nedges;
    fscanf(fid, "%d %d %d", &nverts, &ntriang, &nedges);
    printf("%d, %d, %d\n", nverts, ntriang, nedges);

    Mesh* mesh = new Mesh;
    mesh->numVertices = nverts;
    mesh->numTriangles = ntriang;

    mesh->vertices = new Vertex[nverts];
    mesh->triangles = new Triangle[ntriang];
    mesh->center.x = 0.0;
    mesh->center.y = 0.0;
    mesh->center.z = 0.0;

    int i;
    for(i = 0; i < nverts; i++){
        fscanf(fid, "%f %f %f", &mesh->vertices[i].x, &mesh->vertices[i].y, &mesh->vertices[i].z);
        mesh->center.x += mesh->vertices[i].x;
        mesh->center.y += mesh->vertices[i].y;
        mesh->center.z += mesh->vertices[i].z;
    }

    for(i = 0; i < ntriang; i++){
        int nv;
        fscanf(fid, "%d %d %d %d", &nv, &mesh->triangles[i].indices[0],
                                        &mesh->triangles[i].indices[1],
                                        &mesh->triangles[i].indices[2]);
    }

    fclose(fid);
    mesh->center.x /= nverts;
    mesh->center.y /= nverts;
    mesh->center.z /= nverts;

    float maxx = -1.0e-10, maxy= -1.0e-10, maxz= -1.0e-10;
    float minx = 1.0e10, miny= 1.0e10, minz= 1.0e10;

    for(int i = 0; i < mesh->numVertices; i++){
        if(mesh->vertices[i].x < minx)
            minx = mesh->vertices[i].x;
        if(mesh->vertices[i].x > maxx)
            maxx = mesh->vertices[i].x;
        if(mesh->vertices[i].y < miny)
            miny = mesh->vertices[i].y;
        if(mesh->vertices[i].y > maxy)
            maxy = mesh->vertices[i].y;
        if(mesh->vertices[i].z < minz)
            minz = mesh->vertices[i].z;
        if(mesh->vertices[i].z > maxz)
            maxz = mesh->vertices[i].z;
    }

    float diag = sqrt((maxx-minx)*(maxx-minx) + (maxy-miny)*(maxy-miny)+(maxz-minz)*(maxz-minz));
    mesh->scale = 2.0/diag;

    //mesh->model_transform = glm::mat4(1.0f);
    return *mesh;
}

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
#define MAX_RAY_DEPTH 5

float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}


//Vec3f trace(
//    const Vec3f &rayorig,
//    const Vec3f &raydir,
//    const Mesh* &mesh,
//    const int &depth)
//{
//    float tnear = INFINITY;
//    const Mesh* mesh = NULL;
//
//    for (unsigned i = 0; i< mesh->numTriangles; ++i){
//        float t0 = INFINITY, t1 = INFINITY;
//        if (mesh->triangles[i].intersect(rayorig, raydir, t0, t1)){
//            if (t0 < 0) t0 = t1;
//            if (t0 < tnear) {
//                tnear = t0;
//                triangle = mesh->triangles[i];
//            }
//        }
//    }
//
//
//    float t0 = INFINITY, t1 = INFINITY;
////    if (mesh.
//}
//

Vec3f trace(
    const Vec3f &rayorig,
    const Vec3f &raydir,
    const std::vector<Mesh> &scene,
    const int &depth)
{
    float tnear = INFINITY;
    const Mesh* mesh = NULL;
    // find intersection of this ray with the sphere in the scene
    for (unsigned i = 0; i < scene.size(); ++i){
        float t0 = INFINITY, t1 = INFINITY;
        if ( scene[i].intersect(rayorig, raydir, t0, t1)){
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                mesh = &scene[i];
            }
        }
    }
    // if there's no intersection return black or background color
    //if (#!sphere) return Vec3f(2);
    //else return Vec3f(3);
    //Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
    //Vec3f phit = rayorig + raydir * tnear; // point of intersection
    //Vec3f nhit = phit - sphere->center; // normal at the intersection point
    //nhit.normalize(); // normalize normal direction
    //// If the normal and the view direction are not opposite to each other
    //// reverse the normal direction. That also means we are inside the sphere so set
    //// the inside bool to true. Finally reverse the sign of IdotN which we want
    //// positive.
    //float bias = 1e-4; // add some bias to the point from which we will be tracing
    //bool inside = false;
    //if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
    //if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
    //    float facingratio = -raydir.dot(nhit);
    //    // change the mix value to tweak the effect
    //    float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
    //    // compute reflection direction (not need to normalize because all vectors
    //    // are already normalized)
    //    Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
    //    refldir.normalize();
    //    Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
    //    Vec3f refraction = 0;
    //    // if the sphere is also transparent compute refraction ray (transmission)
    //    if (sphere->transparency) {
    //        float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
    //        float cosi = -nhit.dot(raydir);
    //        float k = 1 - eta * eta * (1 - cosi * cosi);
    //        Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
    //        refrdir.normalize();
    //        refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
    //    }
    //    // the result is a mix of reflection and refraction (if the sphere is transparent)
    //    surfaceColor = (
    //        reflection * fresneleffect +
    //        refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor;
    //}
    //else {
    //    // it's a diffuse object, no need to raytrace any further
    //    for (unsigned i = 0; i < spheres.size(); ++i) {
    //        if (spheres[i].emissionColor.x > 0) {
    //            // this is a light
    //            Vec3f transmission = 1;
    //            Vec3f lightDirection = spheres[i].center - phit;
    //            lightDirection.normalize();
    //            for (unsigned j = 0; j < spheres.size(); ++j) {
    //                if (i != j) {
    //                    float t0, t1;
    //                    if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
    //                        transmission = 0;
    //                        break;
    //                    }
    //                }
    //            }
    //            surfaceColor += sphere->surfaceColor * transmission *
    //            std::max(float(0), nhit.dot(lightDirection)) * spheres[i].emissionColor;
    //        }
    //    }
    //}

    //return surfaceColor + sphere->emissionColor;
}




void render(const std::vector<Mesh> &scene)
{
    unsigned width = 640, height = 480;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);

    for (unsigned y = 0; y < height; ++y){
        for (unsigned x = 0 ; x< width; ++x, ++ pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -1);
            raydir.normalize();
            *pixel = trace(Vec3f(0), raydir, scene, 0); 
        }
    }

    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width*height; ++ i){
        ofs << (unsigned char)(std::min(float(1), image[i].x)*255) <<
               (unsigned char)(std::min(float(1), image[i].y)*255) <<
               (unsigned char)(std::min(float(1), image[i].z)*255);
    }
    ofs.close();
    delete[] image;
}


//Render function
//void render(Mesh* mesh){
//    unsigned width = 640, height=480;
//    Vec3f *image = new Vec3f[width * height], *pixel = image;
//    float invWidth = 1 / float(width), invHeight = 1 / float(height);
//    float fov = 30, aspectratio = width / float(height);
//    float angle = tan(M_PI * 0.5 * fov / 180.);
//
//    for (unsigned y = 0; y < height; ++y) {
//        for (unsigned x = 0; x < width; ++x, ++pixel){
//            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
//            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
//            Vec3f raydir(xx, yy, -1);
//            raydir.normalize();
//            *pixel = trace(Vec3f(0), raydir, mesh, 0);
//        }
//    }
//
//    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
//    ofs << "P6\n" << width << " " << height << "\n255\n";
//    for (unsigned i = 0; i < width * height; ++i) {
//        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
//               (unsigned char)(std::min(float(1), image[i].y) * 255) <<
//               (unsigned char)(std::min(float(1), image[i].z) * 255);
//    }
//    ofs.close();
//    delete [] image;
//}

//[comment]
// In the main function, we will create the scene which is composed of 5 spheres
// and 1 light (which is also a sphere). Then, once the scene description is complete
// we render that scene, by calling the render() function.
//[/comment]
int main(int argc, char **argv)
{
    printf("Hello\n");
    Mesh mesh = leerOFF("avion.off");
    printf("Mesh loaded\n");
    std::vector<Mesh> scene;
    scene.push_back(mesh);
    printf("Rendering scene\n");
    render(scene);
    return 0;
}

