#include "VolumeViewer.h"
#include <iostream>
#include <QtGui>

int main(int argc, char * argv[])
{
    // initialize OSPRay
    ospInit(&argc, (const char **)argv);

    // initialize Qt
    QApplication * app = new QApplication(argc, argv);

    // parse commandline arguments
    // we expect: <filename> <dim1> <dim2> <dim3> <type> <format>
    if(argc < 7)
    {
        std::cerr << "usage: " << argv[0] << " <filename> <dim1> <dim2> <dim3> <format> <layout> [options]" << std::endl;
        std::cerr << "options:" << std::endl;
        std::cerr << "  -dt <dt>                             : use ray cast sample step size 'dt'" << std::endl;
        std::cerr << "  -rotate <rate>                       : automatically rotate view according to 'rate'" << std::endl;
        std::cerr << "  -benchmark <warm-up frames> <frames> : run benchmark and report overall frame rate" << std::endl;
        std::cerr << "  -viewsize <width>x<height>           : force OSPRay view size to 'width'x'height'" << std::endl;

        return 1;
    }

    std::string filename = argv[1];

    osp::vec3i dimensions;
    dimensions.x = atoi(argv[2]);
    dimensions.y = atoi(argv[3]);
    dimensions.z = atoi(argv[4]);

    std::string format = argv[5];
    std::string layout = argv[6];

    std::cout << "got filename = " << filename << ", dimensions = (" << dimensions.x << ", " << dimensions.y << ", " << dimensions.z << "), format = " << format << ", layout = " << layout << std::endl;

    // parse optional command line arguments
    float dt = 0.f;                 // viewer will auto-set dt for dt == 0
    float rotationRate = 0.f;
    int benchmarkWarmUpFrames = 0;
    int benchmarkFrames = 0;
    int viewSizeWidth = 0;
    int viewSizeHeight = 0;

    for(int i=7; i<argc; i++)
    {
        std::string arg = argv[i];

        if(arg == "-dt")
        {
            if(i+1 >= argc)
            {
                throw std::runtime_error("missing <dt> argument");
            }

            dt = atof(argv[++i]);

            std::cout << "got dt = " << dt << std::endl;
        }
        else if(arg == "-rotate")
        {
            if(i+1 >= argc)
            {
                throw std::runtime_error("missing <rate> argument");
            }

            rotationRate = atof(argv[++i]);

            std::cout << "got rotationRate = " << rotationRate << std::endl;
        }
        else if(arg == "-benchmark")
        {
            if(i+2 >= argc)
            {
                throw std::runtime_error("missing <warm-up frames> <frames> arguments");
            }

            benchmarkWarmUpFrames = atoi(argv[++i]);
            benchmarkFrames = atoi(argv[++i]);

            std::cout << "got benchmarkWarmUpFrames = " << benchmarkWarmUpFrames << ", benchmarkFrames = " << benchmarkFrames << std::endl;
        }
        else if(arg == "-viewsize")
        {
            if(i+1 >= argc)
            {
                throw std::runtime_error("missing <width>x<height> argument");
            }

            std::string arg2(argv[++i]);

            size_t pos = arg2.find("x");

            if(pos != std::string::npos)
            {
                arg2.replace(pos, 1, " ");

                std::stringstream ss(arg2);
                ss >> viewSizeWidth >> viewSizeHeight;

                std::cout << "got viewSizeWidth = " << viewSizeWidth << ", viewSizeHeight = " << viewSizeHeight << std::endl;
            }
            else
            {
                throw std::runtime_error("improperly formatted <width>x<height> argument");
            }
        }
        else
        {
            throw std::runtime_error("unknown parameter " + arg);
        }
    }

    // create volume viewer window
    VolumeViewer * volumeViewer = new VolumeViewer();

    // load the volume from the commandline arguments
    volumeViewer->loadVolume(filename, dimensions, format, layout, dt);

    // set rotation rate
    volumeViewer->getQOSPRayWindow()->setRotationRate(rotationRate);

    // set benchmarking parameters
    volumeViewer->getQOSPRayWindow()->setBenchmarkParameters(benchmarkWarmUpFrames, benchmarkFrames);

    // set OSPRay view size if specified
    if(viewSizeWidth != 0 && viewSizeHeight != 0)
    {
        volumeViewer->getQOSPRayWindow()->setFixedSize(viewSizeWidth, viewSizeHeight);
    }

    // enter Qt event loop
    app->exec();

    // cleanup
    delete volumeViewer;

    return 0;
}
