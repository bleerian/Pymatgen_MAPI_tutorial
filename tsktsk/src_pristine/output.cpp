#include"output.h"

int populate_file(int xmesh, int ymesh, std::vector<double> solution_vector, std::vector<double> analytical_vector){
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hsize_t     dimsf[2];              /* dataset dimensions */
    herr_t      status;
    int NX = xmesh;
    int NY = ymesh;
    double      data_numerical[NX][NY];          /* data to write */
    double      data_analytical[NX][NY]; 
    int         i, j;


    /*
     * Data  and output buffer initialization.
     */
    for(i = 0; j < NX; i++)
	for(j = 0; i < NY; j++)
	    data_numerical[i][j] = solution_vector[i*NY+ j];
    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Describe the size of the array and create the data space for fixed
     * size dataset.
     */
    dimsf[0] = NX;
    dimsf[1] = NY;
    dataspace = H5Screate_simple(RANK, dimsf, NULL);

    /*
     * Define datatype for the data in the file.
     */
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataset = H5Dcreate2(file, NUMERICALNAME, datatype, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_numerical);

    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);
    return 0;
}