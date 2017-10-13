#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstdio>

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "shader_utils.h"

GLuint vbo_surface;
GLuint vbo_color;
GLint uniform_mvp;

GLuint program;

GLint attribute_coord3d;
GLint attribute_color;

GLfloat* surface_vertices;
GLfloat* surface_color;

int screen_width = 800, screen_height = 800;

int numPoints;

int nSteps = 200;
int factorial[4] = {1, 1, 2, 6};

float binomial(int n, int i){
    return 1.0*(factorial[n]/(factorial[i] * factorial[n-i]));
}

float bernstein(int n, int i, float t){
    return binomial(n,i) * powf(t,i) * powf(1-t,n-i);
}

std::vector<float> build_patch_points(GLfloat cp[4][4][3]){
    float u = 0.0;
    float v;

    std::vector<float> p;
    
    for(int i = 0; i <= nSteps ; i++){ //Recorrer el parámetro u
        v = 0.0;
        for(int j = 0; j <= nSteps; j++){ //Recorrer el parámetro v
            float x=0.0, y=0.0, z=0.0;
            
            for(int k = 0; k < 4; k++){ //Sumatorias en evaluación de Bezier
                float bi = bernstein(3, k, u);
                for(int l = 0; l < 4; l++){
                    float bj = bernstein(3, l, v);
                    x += bi * bj * cp[k][l][0];
                    y += bi * bj * cp[k][l][1];
                    z += bi * bj * cp[k][l][2];
                }
            }

            p.push_back(x);
            p.push_back(y);
            p.push_back(z);

            v += 1.0/nSteps;
        }

        u += 1.0/nSteps;
    }
    return p;
}


bool init_resources(){

    /*EXAMEN: Asignar memoria y valores a surface_vertices y surface_color*/

    FILE* fin = fopen("control.txt", "rt");
    
    int num_control_points;
    int num_control_patches;
    fscanf(fin, "%d %d", &num_control_points, &num_control_patches);

    // Save control points in an auxiliar variable
    GLfloat control_points[num_control_points][3];
    for (int i=0;i<num_control_points;i++){
        fscanf(fin, "%f %f %f", &control_points[i][0], &control_points[i][1], &control_points[i][2]);
    }

    // Control points grouped by patch
    GLfloat cp[num_control_patches][4][4][3];
    for (int i=0;i<num_control_patches;i++){
        for (int j=0;j<4;j++){
            for (int k=0;k<4;k++){
                int idx;
                fscanf(fin, "%d", &idx);
                for (int l=0;l<3;l++){
                    cp[i][j][k][l] = control_points[idx][l];
                }
            }
        }
    }

    // p vector contains the points to draw
    std::vector<float> p;
    for (int i=0;i<num_control_patches;i++){
        std::vector<float> patch_points = build_patch_points(cp[i]);
        // Concatenate patch points to the end of p vector
        p.insert(p.end(), patch_points.begin(), patch_points.end());
    }

    surface_vertices =  new GLfloat[p.size()];
    surface_color = new GLfloat[p.size()];

    numPoints = p.size()/3;

    // Assign point to surface vertices array;
    for(int i = 0; i < p.size(); i++){
        surface_vertices[i] = p[i];
    }

    // Assign color to surface color array
    for(int i = 0; i < p.size()/3; i++){
        surface_color[3 * i] = 0.0;
        surface_color[3 * i + 1] = 1.0;
        surface_color[3 * i + 2] = 0.0;
    }


    glGenBuffers(1, &vbo_surface);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_surface);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*p.size(), surface_vertices, GL_STATIC_DRAW);

    glGenBuffers(1, &vbo_color);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*p.size(), surface_color, GL_STATIC_DRAW);

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

    attribute_coord3d = glGetAttribLocation(program, "coord3d");
    if(attribute_coord3d == -1){
        std::cout << "No se puede asociar el atributo coord3d" << std::endl;
        return false;
    }

    attribute_color = glGetAttribLocation(program, "color");
    if(attribute_color == -1){
        std::cout << "No se puede asociar el atributo color" << std::endl;
        return false;
    }

    uniform_mvp = glGetUniformLocation(program, "mvp");
    if(uniform_mvp == -1){
        std::cout << "No se puede asociar el uniform mvp" << std::endl;
        return false;
    }

    return true;
}

void onDisplay(){
    //Creamos matrices de modelo, vista y proyeccion
    glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, -1.5));
    glm::mat4 view  = glm::lookAt(glm::vec3(4.0f, 4.0f, 4.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    glm::mat4 projection = glm::perspective(45.0f, 1.0f*screen_width/screen_height, 0.1f, 10.0f);
    glm::mat4 mvp = projection * view * model;

    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    glUseProgram(program);
    //Enviamos la matriz que debe ser usada para cada vertice
    glUniformMatrix4fv(uniform_mvp, 1, GL_FALSE, glm::value_ptr(mvp));

    glEnableVertexAttribArray(attribute_coord3d);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_surface);

    glVertexAttribPointer(
        attribute_coord3d,
        3,
        GL_FLOAT,
        GL_FALSE,
        0, 0
    );

    glEnableVertexAttribArray(attribute_color);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_color);

    glVertexAttribPointer(
        attribute_color,
        3,
        GL_FLOAT,
        GL_FALSE,
        0, 0
    );

    glBindBuffer(GL_ARRAY_BUFFER, vbo_surface);

    glDrawArrays(GL_POINTS, 0, numPoints);
    //glDrawArrays(GL_POINTS, 0, /*EXAMEN: Numero de puntos*/);

    glDisableVertexAttribArray(attribute_coord3d);
    glDisableVertexAttribArray(attribute_color);
    glutSwapBuffers();
}

void onReshape(int w, int h){
    screen_width = w;
    screen_height = h;

    glViewport(0,0,screen_width, screen_height);
}

void free_resources(){
    glDeleteProgram(program);
    glDeleteBuffers(1, &vbo_surface);
    glDeleteBuffers(1, &vbo_color);


    /*EXAMEN - OPCIONAL: elimine la memoria que se utilizó en los buffers. Use delete o free en los arreglos
                surface_vertices y surface_color
    */

    delete[] surface_vertices;
    delete[] surface_color;
}

int main(int argc, char* argv[]){
    glutInit(&argc, argv);
    glutInitContextVersion(2,0);
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(screen_width, screen_height);
    glutCreateWindow("OpenGL");

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
        glutMainLoop();
    }

    free_resources();
    exit(EXIT_SUCCESS);
}
