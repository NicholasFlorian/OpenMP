
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
#include <time.h>


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

hittable *random_scene(int threadTotal) {
    int n = 500;
    hittable **list = new hittable*[n+1];
    list[0] =  new sphere(vec3(0,-1000,0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
    int i = 1;
    
    // run in parallel
    #pragma omp parallel for num_threads(threadTotal)
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            float choose_mat = random_double();
            vec3 center(a+0.9*random_double(),0.2,b+0.9*random_double());
            if ((center-vec3(4,0.2,0)).length() > 0.9) {
                
                // race condition is preformed in a critical section
                // only one thread should have access to list at a time
                #pragma omp critical
                {
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
    }

    list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
    list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
    list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

    return new hittable_list(list,i);
}

int main(int argumentSize, char* argumentArray[]) {
    
    // thread variables
    int threadTotal = 1;    // total threads default 1
    int check = 0;          // global locking variable

    // file variables
    int doOutput = 0;       // control output
    std::ofstream file;     // output file object

    // standard variables
    int nx = 640;           // width of the image
    int ny = 480;           // height of the image

    int ns = 10;

    // timing variables
    time_t start;           // start time
    time_t end;             // end time
    double totalTime;       // total time


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

            // kill the program if arguments are invalid
            std::cout << "Invalid Arguments\n";
            exit(1);
        }
    }

    // print out arguments
    std::cout << "Size {" << nx << ", " << ny << "}, Output " << doOutput << ", Threads " << threadTotal << "\n";


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    // run main simulation
   
    // start timing
    time(&start);

    // simulation variables
    hittable *world = random_scene(threadTotal);
    vec3 lookfrom(13,2,3);
    vec3 lookat(0,0,0);
    float dist_to_focus = 10.0;
    float aperture = 0.1;

    // camera object
    camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, dist_to_focus);


    // create file
    if(doOutput){

        // create the file and print the first line
        file.open("image.ppm");
        file << "P3\n" << nx << " " << ny << "\n255\n";
    }

    // start running in parallel
    for (int j = ny-1; j >= 0; j--) {
        
        // run pragma on inner loop with block
        #pragma omp parallel for schedule(dynamic, 1) num_threads(threadTotal) 
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
            
            // manual scheduled barrier 
            while(currentThread != check);

            // race condition is preformed in a critical section
            #pragma omp critical
            {
                // prefrom output
                if(doOutput)
                    file << ir << " " << ig << " " << ib << "\n";
            
                // unlock next thread
                if(currentThread == threadTotal - 1)
                    check = 0;
                else
                    check++;   
            } 
        }

    }
    

    // end timing
    time(&end);
    totalTime = double(end - start); 
    std::cout << "Time elapsed: " << totalTime << "\n";

    // close the file
    if(doOutput)
        file.close();
}