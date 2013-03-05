//============================================================================
// Name        : Information.hpp
// Author      : Jonathan Webster & Max Revitt
// Description : 2nd Assessed Coursework for ACA 2013
//============================================================================

#ifndef INFORMATION_HPP_
#define INFORMATION_HPP_

#include <iostream>
#include <iomanip>

/* Information about current system & GPUs */
void platformInfo() {

    unsigned i, j;
    char* value;
    size_t valueSize;
    cl_uint platformCount;
    cl_platform_id* platforms;
    cl_uint deviceCount;
    cl_device_id* devices;
    cl_uint maxComputeUnits;

    // get all platforms
    clGetPlatformIDs(0, NULL, &platformCount);
    platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
    clGetPlatformIDs(platformCount, platforms, NULL);

    for (i = 0; i < platformCount; i++) {

        // get all devices
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
        devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);

        // for each device print critical attributes
        for (j = 0; j < deviceCount; j++) {

            // print device name
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
            printf("%d. Device: %s\n", j+1, value);
            free(value);

            // print hardware device version
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
            printf(" %d.%d Hardware version: %s\n", j+1, 1, value);
            free(value);

            // print software driver version
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
            printf(" %d.%d Software version: %s\n", j+1, 2, value);
            free(value);

            // print c version supported by compiler for device
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
            printf(" %d.%d OpenCL C version: %s\n", j+1, 3, value);
            free(value);

            // print parallel compute units
            clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
                    sizeof(maxComputeUnits), &maxComputeUnits, NULL);
            printf(" %d.%d Parallel compute units: %d\n", j+1, 4, maxComputeUnits);

        }

        free(devices);

    }

    free(platforms);
}

void reportSmoothHeaders(Quality initialQ) {

	
	std::cout	<< std::setw(10) << "Initial";
	std::cout	<< " #"; // separator

	std::cout	<< std::setw(5) << "   ";
	std::cout	<< " #"; // separator

	std::cout	<< std::setw(10) << initialQ.mean;
	std::cout	<< " #"; // separator

	std::cout	<< std::setw(10) << initialQ.min;
	std::cout	<< " #\n"; // separator

	std::cout << " Version   #      #    Qmin   #   Qmean   #  TEST #   SCORE    #\n";

}


void reportSmooth(Mesh* mesh, Timer* time, std::string name) {

		std::cout	<< std::setw(10) << name;
		std::cout	<< " #"; // separator
		
		std::cout	<< std::setw(5) << "   ";
		std::cout	<< " #"; // separator

	Quality q = mesh->get_mesh_quality();

		std::cout	<< std::setw(10) << q.mean;
		std::cout	<< " #"; // separator

		std::cout	<< std::setw(10) << q.min;
		std::cout	<< " #"; // separator

  if((q.mean>0.90)&&(q.min>0.55)) {
		std::cout	<< std::setw(6) << "PASS";
		std::cout	<< " #"; // separator
  }else{
		std::cout	<< std::setw(6) << "FAIL";
		std::cout	<< " #"; // separator
  }

		std::cout	<< std::setw(10) << time->getTime();
		std::cout	<< "s #\n"; // separator


}


#endif