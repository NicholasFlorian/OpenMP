
//==================================================================================================
// Written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is distributed
// without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication along
// with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==================================================================================================

#include "camera.h"
#include "hittable_list.h"
#include "material.h"
#include "random.h"
#include "sphere.h"

#include <float.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <omp.h>



vec3 color(const ray& r, hittable *world, int depth) {
    hit_record rec;
    if (world->hit(r, 0.001, MAXFLOAT, rec)) {
        ray scattered;
        vec3 attenuation;
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
             return attenuation*color(scattered, world, depth+1);
        }
        else {
            return vec3(0,0,0);
        }
    }
    else {
        vec3 unit_direction = unit_vector(r.direction());
        float t = 0.5*(unit_direction.y() + 1.0);
        return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
    }
}


hittable *random_scene() {
    int n = 500;
    hittable **list = new hittable*[n+1];
    list[0] =  new sphere(vec3(0,-1000,0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
    int i = 1;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            float choose_mat = random_double();
            vec3 center(a+0.9*random_double(),0.2,b+0.9*random_double());
            if ((center-vec3(4,0.2,0)).length() > 0.9) {
                if (choose_mat < 0.8) {  // diffuse
                    list[i++] = new sphere(
                        center, 0.2,
                        new lambertian(vec3(random_double()*random_double(),
                                            random_double()*random_double(),
                                            random_double()*random_double()))
                    );
                }
                else if (choose_mat < 0.95) { // metal
                    list[i++] = new sphere(
                        center, 0.2,
                        new metal(vec3(0.5*(1 + random_double()),
                                       0.5*(1 + random_double()),
                                       0.5*(1 + random_double())),
                                  0.5*random_double())
                    );
                }
                else {  // glass
                    list[i++] = new sphere(center, 0.2, new dielectric(1.5));
                }
            }
        }
    }

    list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
    list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
    list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

    return new hittable_list(list,i);
}

int* initSplits(int threadTotal, int xSize){

    // variables
    double splitSize;       // the ideal split size for each section
    int currentSize;        // tracks the size of the current split
    int currentSplit;       // current split position
    int processSize;        // the current size of the process

    int* splits;

    // allocate memory for splits
    processSize = threadTotal;
    splits = (int*) malloc(sizeof(int) * (processSize + 1));


    // assign the first and last memebrs of the split
    splits[0] = 0;
    splits[processSize] = xSize;
    
    // assign the rest of the split positions
    splitSize =     xSize / processSize;    // calculate the split size
    currentSize =   0;                      // assign the current size to 0 to start
    currentSplit =  1;                      // assign currentSplit to 1 before starting the loop
    
    
    // calculate the splits by summing approxamate splits
    for(int i = 0; i < xSize; i++) {

        // increase the current size
        currentSize++;

        // currentSize 
        if(currentSize >= splitSize){
        
            // update the split index the position
            splits[currentSplit++] = i + 1;

            // reset the current split
            currentSize = 0;
        }
    }


    // assign the last split to 0
    splits[processSize] = xSize;

    for(int i = 0; i < threadTotal + 1; i++)
         std::cout << splits[i] << "\n";
    



    return splits;
}

int main(int argumentSize, char* argumentArray[]) {
    
    // thread variables
    int threadTotal = 0;
    int* splits;

    // file variables
    int doOutput = 0;
    std::ofstream file;

    // standard variables
    int nx = 640;
    int ny = 480;

    int ns = 10;

    // the world itself
    hittable *world = random_scene();

    // graphics variables
    vec3 lookfrom(13,2,3);
    vec3 lookat(0,0,0);
    float dist_to_focus = 10.0;
    float aperture = 0.1;

    // camera object
    camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, dist_to_focus);


    // handle arguments
    for(int i = 1; i < argumentSize; i++){

        if(strcmp(argumentArray[i], "-size") == 0){

            nx = atoi(argumentArray[++i]);
            ny = atoi(argumentArray[++i]);
        }
        else if(strcmp(argumentArray[i], "-output") == 0){

            doOutput = 1;
        }
        else if(strcmp(argumentArray[i], "-threads") == 0){
            
            threadTotal = atoi(argumentArray[++i]);
        }
        else{

            std::cout << "Invalid Arguments\n";
            exit(1);
        }

    }

    // print out arguments
    std::cout << "Size {" << nx << ", " << ny << "}, Output " << doOutput << ", Threads " << threadTotal << "\n";

    // create splits
    splits = initSplits(threadTotal, nx);

    // create file
    if(doOutput){

        // create the file and print the first line
        file.open("image.ppm");
        file << "P3\n" << nx << " " << ny << "\n255\n";
    }


    // TEST CODE ?????////
    int check = 0;

    // start running in parallel
    for (int j = ny-1; j >= 0; j--) {
        
        // run pragma on inner loop with block
        #pragma omp parallel for num_threads(threadTotal)
        for (int i = 0; i < nx; i++) {
            
            int currentThread;
            currentThread = omp_get_thread_num();

            vec3 col(0, 0, 0);
            for (int s=0; s < ns; s++) {
                float u = float(i + random_double()) / float(nx);
                float v = float(j + random_double()) / float(ny);
                ray r = cam.get_ray(u, v);
                col += color(r, world,0);
            }
            col /= float(ns);
            col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
            int ir = int(255.99*col[0]);
            int ig = int(255.99*col[1]);
            int ib = int(255.99*col[2]);
            
           // #pragma omp barrier

            while(currentThread != check);

            #pragma omp critical
            {

                std::cout << currentThread << "\n";
                if(doOutput)
                    file << ir << " " << ig << " " << ib << "\n";
            
            
                if(currentThread == threadTotal - 1)
                    check = 0;
                else
                    check++;   
            } 
        }

    }
    
    // close the file
    if(doOutput)
        file.close();
}