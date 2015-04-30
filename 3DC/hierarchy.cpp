#include <iostream>
#include <iomanip>
#include <memory>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <fstream>
#include <unordered_map>
#include <string>
#include <GL/gl.h>





#include "FirstLayer.h"

using namespace cv;
using namespace std;


int main( int argc, char** argv )
{


    if( argc < 2)
    {
     cout <<" Usage: display_image ImageToLoadAndDisplay" << endl;
     return -1;
    }

    String fileName = "";
    int kernelSize = 15;



    if( argv[1] == std::string("cone")){
        fileName = "cone.ply";

    }

    else if( argv[1] == std::string("cube")){
        cout << argv[1];
        fileName = "../simple_objects/my_cube.png";
    }

    else{
        cout <<  "No such image" << std::endl ;
        return -1;
    }
    if (argc > 2){
        if ((argv[2] == std::string("-k")) || (argv[2] == std::string("--kernel"))) {
            cout << kernelSize << endl;
            kernelSize = atof(argv[3]); // Increment 'i' so we don't get the argument as the next argv[i].
            cout << kernelSize << endl;
        }
    }

    FirstLayer first;  //Base class of our Hierarchy

    first.readDiscretization("crude_400.txt");

    std::vector<std::vector<Eigen::Vector3f> > pointsSeq;
    std::vector<std::vector<Eigen::Vector3f> > normalsSeq;
    std::vector<std::vector<Eigen::Vector4i> > facesSeq;
    std::vector<std::vector<Eigen::Vector3f> > quantizedNormalsSeq;
    for(int i = 0; i < 1; i++){
        std::vector<Eigen::Vector3f> points;
        std::vector<Eigen::Vector3f> normals;
        std::vector<Eigen::Vector4i> faces;
        std::vector<Eigen::Vector3f> quantizedNormals;
        first.readPlyFile(fileName,points,normals,faces);
        first.quantizeNormals(normals,quantizedNormals);
        pointsSeq.push_back(points);
        normalsSeq.push_back(normals);
        facesSeq.push_back(faces);
        quantizedNormalsSeq.push_back(quantizedNormals);
        first.writePlyToFile("coneQuantized.ply",points,quantizedNormals,faces);
    }



    waitKey(0);                                          // Wait for a keystroke in the window

    return 0;
}

