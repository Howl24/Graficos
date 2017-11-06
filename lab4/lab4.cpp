#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glew.h>

#include <GL/freeglut.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "shader_utils.h"
#include "SOIL.h"

#include "res_texture.c"

#include <vector>
#include <iostream>

#include <typeinfo>

int screen_width=800, screen_height=600;
GLuint vbo_cube_vertices;
GLuint ibo_cube_elements;
GLuint program;

GLint attribute_coord3d;
GLint uniform_mvp;

GLuint vbo_cube_texcoords;
GLint attribute_texcoord;

GLuint texture_id;
GLint uniform_texture;

unsigned char* imag;


// Tipos de datos para representar el objeto
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

    //Matriz de transformación
    glm::mat4 model_transform;

    //Buffers para graficado
    GLfloat* object_vertices;
    GLushort* object_indexes;

    //Id's para buffers
    GLuint vbo_object;
    GLuint ibo_object;
}Mesh;

typedef struct Scene{
    int numMeshes;
    Mesh* meshes[0];
}Scene;

Scene scene;


//Funcion para lectura de archivo obj


Vertex leerVertice(FILE* fid){
  float x, y, z;
  fscanf(fid, "%f %f %f", &x, &y, &z);
 
  Vertex v;
  v.x = x; v.y = y; v.z = z;
  return v;
}

Mesh* leerOBJ(const char* filename){

  FILE* fid = fopen(filename, "rt");

  std::vector<Vertex> vertices;
  std::vector<Vertex> normales;
  std::vector<Vertex> texturas;

  std::vector<Triangle> triangulos;
  char line_type[10];
  while (fscanf(fid, "%s", &line_type)!=EOF){
    if (strcmp(line_type, "v")==0){
      vertices.push_back(leerVertice(fid));
    }

    if (strcmp(line_type, "vn") == 0){
      normales.push_back(leerVertice(fid));
    }

    if (strcmp(line_type, "vt") == 0){
      // Solo se utilizan las dos primera coordenadas (u, v);
      texturas.push_back(leerVertice(fid));
    }

    if (strcmp(line_type, "g") == 0 || strcmp(line_type, "s") == 0){
      char foo[100];
      fscanf(fid, "%s", &foo);
    }

    if (strcmp(line_type, "f") == 0){
      int iv, it, in;
      char last_char;
      std::vector< std::vector<int> > side_indexes;
      std::cout << "Foo\n";
      while(fscanf(fid, "%d %d %d%c", &iv, &it, &in, &last_char)){
        std::vector<int> indexes;
        indexes.push_back(iv);
        indexes.push_back(it);
        indexes.push_back(in);

        std::cout << "V: " << iv << it << in << std::endl;
        side_indexes.push_back(indexes);
        if (last_char == 13) break;
      }

      if (int(side_indexes.size()) == 3){
        Triangle t;
        t.indices[0] = side_indexes[0][0];
        t.indices[1] = side_indexes[1][0];
        t.indices[2] = side_indexes[2][0];

        triangulos.push_back(t);
      }

      if (int(side_indexes.size()) == 4){
        Triangle t1, t2;

        t1.indices[0] = side_indexes[0][0];
        t1.indices[1] = side_indexes[1][0];
        t1.indices[2] = side_indexes[2][0];

        t2.indices[0] = side_indexes[0][0];
        t2.indices[1] = side_indexes[2][0];
        t2.indices[2] = side_indexes[3][0];

        triangulos.push_back(t1);
        triangulos.push_back(t2);
      }
    }
  }

  Mesh* mesh = new Mesh;
  mesh->numVertices = vertices.size();
  mesh->numTriangles = triangulos.size();

  mesh->vertices = new Vertex[vertices.size()];
  mesh->triangles = new Triangle[triangulos.size()];
  mesh->center.x = 0.0;
  mesh->center.y = 0.0;
  mesh->center.z = 0.0;

  for (int i=0;i<vertices.size() ;i++){
    mesh->vertices[i] = vertices[i];
    mesh->center.x += mesh->vertices[i].x;
    mesh->center.y += mesh->vertices[i].y;
    mesh->center.z += mesh->vertices[i].z;
    //printf("%f\n", mesh->vertices[i].x);
  }




  return mesh;
}

int init_resources()
{
  scene.meshes[0] = leerOBJ("House_3_AO.obj");

  GLfloat cube_vertices[] = {
    // front
    -1.0, -1.0,  1.0,
     1.0, -1.0,  1.0,
     1.0,  1.0,  1.0,
    -1.0,  1.0,  1.0,
    // top
    -1.0,  1.0,  1.0,
     1.0,  1.0,  1.0,
     1.0,  1.0, -1.0,
    -1.0,  1.0, -1.0,
    // back
     1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0,
    -1.0,  1.0, -1.0,
     1.0,  1.0, -1.0,
    // bottom
    -1.0, -1.0, -1.0,
     1.0, -1.0, -1.0,
     1.0, -1.0,  1.0,
    -1.0, -1.0,  1.0,
    // left
    -1.0, -1.0, -1.0,
    -1.0, -1.0,  1.0,
    -1.0,  1.0,  1.0,
    -1.0,  1.0, -1.0,
    // right
     1.0, -1.0,  1.0,
     1.0, -1.0, -1.0,
     1.0,  1.0, -1.0,
     1.0,  1.0,  1.0,
  };

  glGenBuffers(1, &vbo_cube_vertices);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_cube_vertices);
  glBufferData(GL_ARRAY_BUFFER, sizeof(cube_vertices), cube_vertices, GL_STATIC_DRAW);

  GLfloat cube_texcoords[48];

  cube_texcoords[0] = 0.0;  cube_texcoords[1] = 0.0;
  cube_texcoords[2] = 1.0;  cube_texcoords[3] = 0.0;
  cube_texcoords[4] = 1.0;  cube_texcoords[5] = 1.0;
  cube_texcoords[6] = 0.0;  cube_texcoords[7] = 1.0;

  for(int i = 8; i < 48; i++){
    cube_texcoords[i] = cube_texcoords[i%8];
  }

  glGenBuffers(1, &vbo_cube_texcoords);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_cube_texcoords);
  glBufferData(GL_ARRAY_BUFFER, sizeof(cube_texcoords),
               cube_texcoords, GL_STATIC_DRAW);

    GLushort cube_elements[] = {
    // front
    0,  1,  2,
    2,  3,  0,
    // top
    4,  5,  6,
    6,  7,  4,
    // back
    8,  9, 10,
    10, 11,  8,
    // bottom
    12, 13, 14,
    14, 15, 12,
    // left
    16, 17, 18,
    18, 19, 16,
    // right
    20, 21, 22,
    22, 23, 20,
  };

  glGenBuffers(1, &ibo_cube_elements);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_cube_elements);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cube_elements), cube_elements, GL_STATIC_DRAW);

  glGenTextures(1, &texture_id);
  glBindTexture(GL_TEXTURE_2D, texture_id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                  GL_LINEAR);

  int width, height, channel;
  imag = SOIL_load_image("pattern2.jpg", &width, &height, &channel, SOIL_LOAD_AUTO);

  glTexImage2D(GL_TEXTURE_2D,
               0,
               GL_RGB,
               width,
               height,
               0,
               GL_RGB,
               GL_UNSIGNED_BYTE,
               imag);


  GLint link_ok = GL_FALSE;

  GLuint vs, fs;
  if ((vs = create_shader("cube.v.glsl", GL_VERTEX_SHADER))   == 0) return 0;
  if ((fs = create_shader("cube.f.glsl", GL_FRAGMENT_SHADER)) == 0) return 0;

  program = glCreateProgram();
  glAttachShader(program, vs);
  glAttachShader(program, fs);
  glLinkProgram(program);
  glGetProgramiv(program, GL_LINK_STATUS, &link_ok);
  if (!link_ok) {
    fprintf(stderr, "glLinkProgram:");
    print_log(program);
    return 0;
  }

  const char* attribute_name;
  attribute_name = "coord3d";
  attribute_coord3d = glGetAttribLocation(program, attribute_name);
  if (attribute_coord3d == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return 0;
  }

  attribute_texcoord = glGetAttribLocation(program,
                                           "texcoord");
  if(attribute_texcoord == -1){
    fprintf(stderr, "Could not bind attribute texcoord\n");
    return 0;
  }

  const char* uniform_name;
  uniform_name = "mvp";
  uniform_mvp = glGetUniformLocation(program, uniform_name);
  if (uniform_mvp == -1) {
    fprintf(stderr, "Could not bind uniform %s\n", uniform_name);
    return 0;
  }

  uniform_texture = glGetUniformLocation(program, "mytexture");
  if (uniform_texture == -1) {

    fprintf(stderr, "Could not bind uniform mytexture\n");
    return 0;
  }

  return 1;
}

void onIdle() {
  float angle = glutGet(GLUT_ELAPSED_TIME) / 1000.0 * glm::radians(15.0);  // base 15° per second
  glm::mat4 anim = \
    glm::rotate(glm::mat4(1.0f), angle*3.0f, glm::vec3(1, 0, 0)) *  // X axis
    glm::rotate(glm::mat4(1.0f), angle*2.0f, glm::vec3(0, 1, 0)) *  // Y axis
    glm::rotate(glm::mat4(1.0f), angle*4.0f, glm::vec3(0, 0, 1));   // Z axis

  glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, -4.0));
  glm::mat4 view = glm::lookAt(glm::vec3(0.0, 2.0, 0.0), glm::vec3(0.0, 0.0, -4.0), glm::vec3(0.0, 1.0, 0.0));
  glm::mat4 projection = glm::perspective(45.0f, 1.0f*screen_width/screen_height, 0.1f, 10.0f);

  glm::mat4 mvp = projection * view * model * anim;
  glUseProgram(program);
  glUniformMatrix4fv(uniform_mvp, 1, GL_FALSE, glm::value_ptr(mvp));
  glutPostRedisplay();
}

void onDisplay()
{
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  glUseProgram(program);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, texture_id);
  glUniform1i(uniform_texture, 0);

  glEnableVertexAttribArray(attribute_coord3d);
  // Describe our vertices array to OpenGL (it can't guess its format automatically)
  glBindBuffer(GL_ARRAY_BUFFER, vbo_cube_vertices);
  glVertexAttribPointer(
    attribute_coord3d, // attribute
    3,                 // number of elements per vertex, here (x,y,z)
    GL_FLOAT,          // the type of each element
    GL_FALSE,          // take our values as-is
    0,                 // no extra data between each position
    0                  // offset of first element
  );

  glEnableVertexAttribArray(attribute_texcoord);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_cube_texcoords);

  glVertexAttribPointer(
    attribute_texcoord,
    2,
    GL_FLOAT,
    GL_FALSE,
    0,0
    );

  /* Push each element in buffer_vertices to the vertex shader */
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_cube_elements);
  int size;  glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
  glDrawElements(GL_TRIANGLES, size/sizeof(GLushort), GL_UNSIGNED_SHORT, 0);

  glDisableVertexAttribArray(attribute_coord3d);
  glDisableVertexAttribArray(attribute_texcoord);
  glutSwapBuffers();
}

void onReshape(int width, int height) {
  screen_width = width;
  screen_height = height;
  glViewport(0, 0, screen_width, screen_height);
}

void free_resources()
{
  glDeleteProgram(program);
  glDeleteBuffers(1, &vbo_cube_vertices);
  glDeleteBuffers(1, &ibo_cube_elements);
  glDeleteBuffers(1, &vbo_cube_texcoords);
  glDeleteTextures(1, &texture_id);
  SOIL_free_image_data(imag);
}


int main(int argc, char* argv[]) {
  glutInit(&argc, argv);
  glutInitContextVersion(2,0);
  glutInitDisplayMode(GLUT_RGBA|GLUT_ALPHA|GLUT_DOUBLE|GLUT_DEPTH);
  glutInitWindowSize(screen_width, screen_height);
  glutCreateWindow("My Textured Cube");

  GLenum glew_status = glewInit();
  if (glew_status != GLEW_OK) {
    fprintf(stderr, "Error: %s\n", glewGetErrorString(glew_status));
    return 1;
  }

  if (!GLEW_VERSION_2_0) {
    fprintf(stderr, "Error: your graphic card does not support OpenGL 2.0\n");
    return 1;
  }

  if (init_resources()) {
    glutDisplayFunc(onDisplay);
    glutReshapeFunc(onReshape);
    glutIdleFunc(onIdle);
    glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glutMainLoop();
  }

  free_resources();
  return 0;
}
