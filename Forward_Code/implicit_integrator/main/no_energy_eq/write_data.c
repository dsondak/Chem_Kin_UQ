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

int create_file(const char* fname, int rows, int cols, int num_scen)
{

    hid_t     file_id;
    hid_t     dataset_id_time, dataset_id_truth, dataset_id_ignition;
    hid_t     dataspace_id_time, dataspace_id_truth;
    hid_t     dataspace_id_ignition, dataspace_id_att;
    hid_t     group_id, attribute_id;
    hsize_t   dims[2], dimsi[1], dimt[1], dimsA[1];
    herr_t    status = 0;

    char group_name[11]; // Need 11 for the null terminator
    char dataset_name_time[15];
    char dataset_name_truth[16];
    char dataset_name_ignition[19];

    char scen_str[2];

    // Create new file using default properties
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Set up dimensions for data space */

    // For truth data
    dims[0] = rows;
    dims[1] = cols;

    // For ignition data
    dimsi[0] = 2;

    // For time data
    dimt[0] = rows;

    // For attributes
    dimsA[0] = 2;

    for (int ns = 0; ns < num_scen; ns++)
    {

        // Create groups
        strcpy(group_name, "/Scenario");
        sprintf(scen_str, "%d", ns+1);
        strcat(group_name, scen_str);
        group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Create the data space for the truth data 
        dataspace_id_truth = H5Screate_simple(2, dims, NULL);

        // Create the data space for the ignition data 
        dataspace_id_ignition = H5Screate_simple(1, dimsi, NULL);

        // Create the data space for time
        dataspace_id_time = H5Screate_simple(1, dimt, NULL);

        // Create the dataset for the time data
        strcpy(dataset_name_time, group_name);
        strcat(dataset_name_time, "/time");

        dataset_id_time = H5Dcreate(file_id, dataset_name_time, H5T_NATIVE_DOUBLE, dataspace_id_time,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Create the dataset for the truth data
        strcpy(dataset_name_truth, group_name);
        strcat(dataset_name_truth, "/truth");

        dataset_id_truth = H5Dcreate(file_id, dataset_name_truth, H5T_NATIVE_DOUBLE, dataspace_id_truth,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Create the dataset for the ignition data
        strcpy(dataset_name_ignition, group_name);
        strcat(dataset_name_ignition, "/ignition");

        dataset_id_ignition = H5Dcreate(file_id, dataset_name_ignition, H5T_NATIVE_DOUBLE, dataspace_id_ignition,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Create dataspace for attributes
        dataspace_id_att = H5Screate_simple(1, dimsA, NULL);

        // Create attribute dataset
        attribute_id = H5Acreate(dataset_id_truth, "Scenario Parameters", H5T_NATIVE_DOUBLE, dataspace_id_att,
                                 H5P_DEFAULT, H5P_DEFAULT);

        // Close the attribute
        status = H5Aclose(attribute_id);

        // Terminate access to the data space
        status += H5Sclose(dataspace_id_time);
        status += H5Sclose(dataspace_id_truth);
        status += H5Sclose(dataspace_id_ignition);
        status += H5Sclose(dataspace_id_att);

        // End access to the dataset and release resources used by it
        status += H5Dclose(dataset_id_time);
        status += H5Dclose(dataset_id_truth);
        status += H5Dclose(dataset_id_ignition);

        // Close the groups
        status += H5Gclose(group_id);

    }

    // Close the file
    status += H5Fclose(file_id);

    return status;
}


int write_file(std::vector<double>& time, std::vector<double>& sols, const char* fname, int rows, int cols, int scen,
                double phi, double scen_par)
{

    hid_t    file_id, dataset_id_truth, dataset_id_time, dataset_id_ignition, attribute_id;
    herr_t   status = 0;

    // Define some file characters
    char group_name[11]; // Need 11 for the null terminator
    char dataset_name_time[16]; // was 15
    char dataset_name_truth[17]; // was 16
    char dataset_name_ignition[20]; // was 19

    char scen_str[2];

    // Initialize solution array
    double truth_data[rows][cols] = {0.0};
    double ignition_data[2] = {0.0};
    double times[rows] = {0.0};

    // Set up solution array
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            truth_data[i][j] = sols[cols * i + j];
        }
    }

    ignition_data[0] = sols[rows * cols];
    ignition_data[1] = sols[rows * cols + 1];

    // Copy time from vector to array (HDF5 C API requires this?)
    for (int i = 0; i < rows; i++)
    {
        times[i] = time[i];
    }

    // Get attribute data
    double attr_data[2];
    attr_data[0] = phi;
    attr_data[1] = scen_par;

    // Open the file
    file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

    /* Open the dataset */
    strcpy(group_name, "/Scenario");

    // Get group name
    sprintf(scen_str, "%d", scen+1);
    strcat(group_name, scen_str);

    // Get time group name
    strcpy(dataset_name_time, group_name);
    strcat(dataset_name_time, "/time");

    // Get truth group name
    strcpy(dataset_name_truth, group_name);
    strcat(dataset_name_truth, "/truth");

    // Get ignition group name
    strcpy(dataset_name_ignition, group_name);
    strcat(dataset_name_ignition, "/ignition");

    // Open datasets and attributes
    dataset_id_time     = H5Dopen(file_id, dataset_name_time, H5P_DEFAULT);
    dataset_id_truth    = H5Dopen(file_id, dataset_name_truth, H5P_DEFAULT);
    dataset_id_ignition = H5Dopen(file_id, dataset_name_ignition, H5P_DEFAULT);
    attribute_id        = H5Aopen(dataset_id_truth, "Scenario Parameters", H5P_DEFAULT);

    // Write the attribute
    status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, attr_data);

    // Write the dataset
    status += H5Dwrite(dataset_id_time, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      times);
    status += H5Dwrite(dataset_id_truth, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      truth_data);
    status += H5Dwrite(dataset_id_ignition, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      ignition_data);

    // Close the dataset
    status += H5Dclose(dataset_id_time);
    status += H5Dclose(dataset_id_truth);
    status += H5Dclose(dataset_id_ignition);
    status += H5Aclose(attribute_id);

    // Close the file
    status += H5Fclose(file_id);

    return status;

}
