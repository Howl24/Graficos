#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include <GL/glew.h>
#include <GL/freeglut.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "shader_utils.h"

GLuint program;
GLint attribute_coord;
GLint uniform_mvp;
GLint uniform_model;

//Variables para el movimiento del avion
GLfloat propeller_speed = 0.0f;
GLfloat airplane_speed = 0.0f;
GLfloat scale = 0.2f;
int num_airplanes = 5;

GLfloat t;

int screen_width = 800, screen_height = 800;

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
    Mesh* meshes[5];
}Scene;


Scene scene;
int numEdges;

Mesh* leerOFF(const char* filename){
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

    mesh->model_transform = glm::mat4(1.0f);
    return mesh;
}

void init_buffers(Mesh* mesh){
    mesh->object_vertices = new GLfloat[mesh->numVertices * 3];
    //mesh->object_color = new GLfloat[mesh->numVertices * 3];
    mesh->object_indexes = new GLushort[mesh->numTriangles * 3];

    int i;

    for(i = 0; i < mesh->numVertices; i++){
        mesh->object_vertices[3 * i] = mesh->vertices[i].x;
        mesh->object_vertices[3 * i + 1] = mesh->vertices[i].y;
        mesh->object_vertices[3 * i + 2] = mesh->vertices[i].z;

    }

    for(i = 0; i < mesh->numTriangles; i++){
        mesh->object_indexes[3 * i] = mesh->triangles[i].indices[0];
        mesh->object_indexes[3 * i + 1] = mesh->triangles[i].indices[1];
        mesh->object_indexes[3 * i + 2] = mesh->triangles[i].indices[2];
    }

    glGenBuffers(1, &mesh->vbo_object);
    glBindBuffer(GL_ARRAY_BUFFER, mesh->vbo_object);
    glBufferData(GL_ARRAY_BUFFER, mesh->numVertices * 3 * sizeof(GLfloat), mesh->object_vertices, GL_STATIC_DRAW);

    glGenBuffers(1, &mesh->ibo_object);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->ibo_object);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh->numTriangles * 3 * sizeof(GLushort), mesh->object_indexes, GL_STATIC_DRAW);
}


bool init_resources(){

    // Leer archivos off para el numero de aviones
    for (int i=0;i< num_airplanes;i++){
        scene.meshes[2*i] = leerOFF("avion.off");
        scene.meshes[2*i+1] = leerOFF("helice.off");
    }

    scene.numMeshes = 2 * num_airplanes;

    // Inicializar buffers segun la cantidad de aviones
    for (int i=0;i< num_airplanes;i++){
        init_buffers(scene.meshes[2*i]);
        init_buffers(scene.meshes[2*i + 1]);
    }

    GLint link_ok = GL_FALSE;
    GLuint vs, fs;
    if((vs = create_shader("basic3.v.glsl", GL_VERTEX_SHADER))==0) return false;
    if((fs = create_shader("basic3.f.glsl", GL_FRAGMENT_SHADER))==0) return false;

    program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);

    glGetProgramiv(program, GL_LINK_STATUS, &link_ok);
    if(!link_ok){
        std::cout << "Problemas con el Shader" << std::endl;
        return false;
    }

    attribute_coord = glGetAttribLocation(program, "coord3d");
    if(attribute_coord == -1){
        std::cout << "No se puede asociar el atributo coord" << std::endl;
        return false;
    }

    uniform_mvp = glGetUniformLocation(program, "mvp");
    if(uniform_mvp == -1){
        std::cout << "No se puede asociar el uniform mvp" << std::endl;
        return false;
    }

    uniform_model = glGetUniformLocation(program, "model");
    if(uniform_model == -1){
        std::cout << "No se puede asociar el uniform model" << std::endl;
        return false;
    }

    return true;
}

void graficarObjeto(Mesh* mesh){
    //Creamos matrices de modelo, vista y proyeccion
    glm::mat4 model =   mesh->model_transform;

    glm::mat4 view  = glm::lookAt(glm::vec3(10.0f, 10.0f, 10.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    glm::mat4 projection = glm::perspective(45.0f, 1.0f*screen_width/screen_height, 0.1f, 100.0f);
    glm::mat4 mvp = projection * view ;

    glUseProgram(program);

    //Enviamos la matriz que debe ser usada para cada vertice
    glUniformMatrix4fv(uniform_mvp, 1, GL_FALSE, glm::value_ptr(mvp));
    glUniformMatrix4fv(uniform_model,1, GL_FALSE, glm::value_ptr(model));

    glEnableVertexAttribArray(attribute_coord);
    glBindBuffer(GL_ARRAY_BUFFER, mesh->vbo_object);

    glVertexAttribPointer(
        attribute_coord,
        3,
        GL_FLOAT,
        GL_FALSE,
        0, 0
    );

    //Dibujar las primitivas
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->ibo_object);
    int size;   glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);

    //Dibujar los triánglos
    glDrawElements(GL_TRIANGLES, size/sizeof(GLushort), GL_UNSIGNED_SHORT, 0);

    glDisableVertexAttribArray(attribute_coord);
}

void onDisplay(){

    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    for(int i = 0; i < scene.numMeshes; i++)
        graficarObjeto(scene.meshes[i]);

    glutSwapBuffers();
}

void onReshape(int w, int h){
    screen_width = w;
    screen_height = h;

    glViewport(0,0,screen_width, screen_height);
}

// Avion del medio
void transform_airplane(Mesh* airplane, Mesh* propeller,
                        GLfloat airplane_speed, GLfloat propeller_speed, GLfloat scale){
    airplane -> model_transform = 
        glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, airplane_speed * -1.0f)) * 
        glm::scale(glm::mat4(1.0f), glm::vec3(scale, scale, scale));

    propeller -> model_transform =
        glm::rotate(glm::mat4(1.0f), propeller_speed, glm::vec3(0.0f, 0.0f, 1.0f)) *
        glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, airplane_speed * -1.0f)) *
        glm::scale(glm::mat4(1.0f), glm::vec3(scale, scale, scale));
}

// Aviones de los costados
void transform_side_airplane(Mesh* airplane, Mesh* propeller,
                             GLfloat airplane_speed, GLfloat propeller_speed, GLfloat scale,
                             int idx){
    float side_sep = 2.0f;
    float front_sep = 2.0f;
    float factor = idx;

    if (idx % 2){
        side_sep = 2.0f;
    }else{
        side_sep = -2.0f;
        factor--;
    }

    airplane -> model_transform = 
        glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, airplane_speed * -1.0f)) * 
        glm::translate(glm::mat4(1.0f), glm::vec3(side_sep * factor, 0.0f, front_sep * factor)) *
        glm::scale(glm::mat4(1.0f), glm::vec3(scale, scale, scale));


    propeller -> model_transform =
        glm::translate(glm::mat4(1.0f), glm::vec3(side_sep * factor, 0.0f, front_sep * factor)) *
        glm::rotate(glm::mat4(1.0f), propeller_speed, glm::vec3(0.0f, 0.0f, 1.0f)) *
        glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, airplane_speed * -1.0f)) *
        glm::scale(glm::mat4(1.0f), glm::vec3(scale, scale, scale));
}

void onIdle(){
    GLfloat t1 = 1.0f * glutGet(GLUT_ELAPSED_TIME);
    GLfloat delta = (t1-t);

    t = t1;
    float animation_factor = 40.0f;

    airplane_speed += (delta/animation_factor) * glm::radians(1.0f) * 5.0f;
    propeller_speed += (delta/animation_factor) * glm::radians(1.0f) * 50.0f;
    scale = 0.2f;

    // Transformacion del avion del centro
    transform_airplane(scene.meshes[0], scene.meshes[1], airplane_speed, propeller_speed, scale);
    for (int i=1; i<num_airplanes;i++){
        // Transformacion de aviones de los costados
        transform_side_airplane(scene.meshes[2*i],
                                scene.meshes[2*i+1], 
                                airplane_speed, 
                                propeller_speed, 
                                scale, i);
    }

    glutPostRedisplay();
}

void free_resources(){
    glDeleteProgram(program);

    for(int i = 0; i < scene.numMeshes; i++){
        glDeleteBuffers(1, &scene.meshes[i]->vbo_object);
        glDeleteBuffers(1, &scene.meshes[i]->ibo_object);
        delete[] scene.meshes[i]->object_vertices;
        delete[] scene.meshes[i]->object_indexes;
        delete[] scene.meshes[i]->vertices;
        delete[] scene.meshes[i]->triangles;
        delete scene.meshes[i];
    }
}

int main(int argc, char* argv[]){


    glutInit(&argc, argv);
    glutInitContextVersion(2,0);
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(screen_width, screen_height);
    glutCreateWindow("OpenGL");

    t = 1.0f*glutGet(GLUT_ELAPSED_TIME);

    GLenum glew_status = glewInit();
    if(glew_status != GLEW_OK){
        std::cout << "Error inicializando GLEW" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(!GLEW_VERSION_2_0){
        std::cout << "Tu tarjeta grafica no soporta OpenGL 2.0" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(init_resources()){
        glutDisplayFunc(onDisplay);
        glutReshapeFunc(onReshape);
        glEnable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glutIdleFunc(onIdle);
        glutMainLoop();
    }

    free_resources();
    exit(EXIT_SUCCESS);
}
