/*-------------------------------------------------------------------
 *
 ***   write_data.c   ***
 *
 * This script is used to write data to the hdf5 file format.
 *
 * This file contains two functions:
 *     1.)  create_file
 *     2.)  write_file
 *
 * 1.) create_file
 *     INPUTS:
 *            o fname:     filename
 *            o rows:      Number of rows in data array
 *            o cols:      Number of columns in data array
 *            o num_scen:  Number of groups in hdf5 file
 *
 *     OUTPUTS:
 *             o status: Error message from hdf
 *
 * 2.) write_file
 *     INPUTS:
 *            o time:      Time data
 *            o sols:      Solution data
 *            o fname:     Filename
 *            o rows:      Number of rows in data array
 *            o cols:      Number of columns in data array
 *            o scen:      Current group (i.e. scenario)
 *            o phi:       Attribute data (equivalence ratio)
 *            o scen_par:  Attribute data (heating rate for now)
 *
 *     OUTPUTS:
 *             o status: Error message from hdf
 *
 * Description of data format
 * ---------------------------
 * The data is organized by scenarios.  Each hdf group corresponds 
 * to one scenario.  A scenario is a tuple of parameters
 * (phi, heating rate) or (phi, initial temperature).  Each group
 * contains two datasets and one attribute set.  The first dataset
 * is the time points at which that data was taken (a vector of time).
 * The second dataset is the actual data stored in the format
 * ntimes X n_species + 1.  The last column is temperature.  The
 * attribute data contains the equivalence ratio and the heating rate.
 * This will eventually be expanded to also include the initial
 * temperature.
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/
#include "write_data.h"
#include <hdf5.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>

int create_file(const char* fname, int rows, int cols)
{

    hid_t     file_id;
    hid_t     dataset_id_fx, dataset_id_fT;
    hid_t     dataspace_id_fx, dataspace_id_fT;
    hid_t     group_id;
    hsize_t   dims[2], dimE[1];
    herr_t    status = 0;

    char group_name[5];      // Need 5 for the null terminator
    char dataset_name_fx[7];
    char dataset_name_fT[7];

    // Create new file using default properties
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Set up dimensions for data space */

    // For species rhs
    dims[0] = rows;
    dims[1] = cols;

    // For T rhs
    dimE[0] = rows;

    // Create groups
    strcpy(group_name, "/rhs");
    group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create the data space for the species rhs
    dataspace_id_fx = H5Screate_simple(2, dims, NULL);

    // Create the data space for energy rhs
    dataspace_id_fT = H5Screate_simple(1, dimE, NULL);

    // Create the dataset for the species rhs
    strcpy(dataset_name_fx, group_name);
    strcat(dataset_name_fx, "/fx");

    dataset_id_fx = H5Dcreate(file_id, dataset_name_fx, H5T_NATIVE_DOUBLE, dataspace_id_fx,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create the dataset for the energy rhs
    strcpy(dataset_name_fT, group_name);
    strcat(dataset_name_fT, "/fT");

    dataset_id_fT = H5Dcreate(file_id, dataset_name_fT, H5T_NATIVE_DOUBLE, dataspace_id_fT,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Terminate access to the data space
    status += H5Sclose(dataspace_id_fx);
    status += H5Sclose(dataspace_id_fT);

    // End access to the dataset and release resources used by it
    status += H5Dclose(dataset_id_fx);
    status += H5Dclose(dataset_id_fT);

    // Close the groups
    status += H5Gclose(group_id);

    // Close the file
    status += H5Fclose(file_id);

    return status;
}


int write_file(std::vector<double>& fx, std::vector<double>& fT, const char* fname, int rows, int cols)
{

    hid_t    file_id, dataset_id_fx, dataset_id_fT;
    herr_t   status = 0;

    // Define some file characters
    char group_name[5];      // Need 5 for the null terminator
    char dataset_name_fx[7];
    char dataset_name_fT[7];

    // Initialize solution array
    double species_rhs[rows][cols] = {0.0};
    double energy_rhs[rows] = {0.0};

    // Set up solution array
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            species_rhs[i][j] = fx[cols * i + j];
        }
        energy_rhs[i] = fT[i]; // Copy time from vector to array 
                               // (HDF5 C API requires this?)
    }

    // Open the file
    file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

    /* Open the dataset */
    strcpy(group_name, "/rhs");

    // Get species rhs group name
    strcpy(dataset_name_fx, group_name);
    strcat(dataset_name_fx, "/fx");

    // Get energy rhs group name
    strcpy(dataset_name_fT, group_name);
    strcat(dataset_name_fT, "/fT");

    // Open datasets and attributes
    dataset_id_fx = H5Dopen(file_id, dataset_name_fx, H5P_DEFAULT);
    dataset_id_fT = H5Dopen(file_id, dataset_name_fT, H5P_DEFAULT);

    // Write the dataset
    status += H5Dwrite(dataset_id_fx, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      species_rhs);
    status += H5Dwrite(dataset_id_fT, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      energy_rhs);

    // Close the dataset
    status += H5Dclose(dataset_id_fx);
    status += H5Dclose(dataset_id_fT);

    // Close the file
    status += H5Fclose(file_id);

    return status;

}
