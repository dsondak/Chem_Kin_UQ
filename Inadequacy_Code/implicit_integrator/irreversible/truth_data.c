/*-------------------------------------------------------------------
 *
 ***   truth_data.c   ***
 *
 * Truth data class.
 *
 * INPUTS:
 *        o fname:       filename
 *        o n_phis:      Number of equivalence ratios
 *        o n_scen:      Number of heating rates (or
 *                       rarely initial temperatures)
 *        o n_times:     Times at which samples were
 *                       taken.
 *        o num_fields:  Number of species + 1 for
 *                       temperature
 *
 * MEMBERS:
 *         o observation_data:  The truth data for species and temperature
 *         o sample_points:     Times at which truth data was sampled
 *         o scenario_params:   Tuple of each scenario in the form
 *                              (phi, heating rate).  Will eventually also
 *                              include initial temperature.
 *
 * Contact:   David Sondak, dsondak@ices.utexas.edu
 *
 *-------------------------------------------------------------------*/

#include "truth_data.h"
#include <hdf5.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>

truth_data::truth_data(const char* fname, int n_phis, int n_scen, int n_times, int num_fields)
:
    observation_data((double *) malloc(n_phis * n_scen * n_times * num_fields * sizeof(double))), // Observation data
    ignition_data((double *) malloc(2 * n_phis * n_scen * sizeof(double))),                       // Ignition time and temperature
    sample_points((double *) malloc(n_scen * n_phis * n_times * sizeof(double))),                 // Sample points
    scenario_params((double *) malloc(2* n_phis * n_scen * sizeof(double)))                       // Scenario parameters (1 for phi and 1 for Q)
{


    // Define hdf5 variables
    hid_t    file_id, dataset_id_truth, dataset_id_time, dataset_id_ignition, attribute_id;
    herr_t   status;

    double * rdata;
    rdata = (double *) malloc(n_times * num_fields * sizeof(double));

    double * idata;
    idata = (double *) malloc(2 * sizeof(double));

    double * tdata;
    tdata = (double *) malloc(n_times * sizeof(double));

    double * adata;
    adata = (double *) malloc(2.0 * sizeof(double));

    char group_name[11];
    char dataset_name_truth[16]; // was 17
    char dataset_name_ignition[19]; // was 20
    char dataset_name_time[15]; // Was 16
    char attribute_name[20];
    char scen_str[2];

    int scen;

    // Get the attribute dataset name
    strcpy(attribute_name, "Scenario Parameters");

    // Open file
    file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

    for (int ii = 0; ii < n_phis; ii++)
    {
        for (int jj = 0; jj < n_scen; jj++)
        {

            scen = ii * n_scen + jj;

            // Get the group name
            strcpy(group_name, "/Scenario");
            sprintf(scen_str, "%d", scen+1);
            strcat(group_name, scen_str);

            // Get the time dataset name
            strcpy(dataset_name_time, group_name);
            strcat(dataset_name_time, "/time");

            // Open time dataset
            dataset_id_time = H5Dopen(file_id, dataset_name_time, H5P_DEFAULT);

            // Read time dataset
            status = H5Dread(dataset_id_time, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                             tdata);

            // Now put tdata into sample_points
            for (int i = 0; i < n_times; i++)
            {
                sample_points[scen * n_times + i] = tdata[i];
            }
            // Get the truth dataset name
            strcpy(dataset_name_truth, group_name);
            strcat(dataset_name_truth, "/truth");

            // Now open the truth dataset
            dataset_id_truth = H5Dopen(file_id, dataset_name_truth, H5P_DEFAULT);

            // Read dataset
            status = H5Dread(dataset_id_truth, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                             rdata);

            // Now put rdata into observation_data
            for (int k = 0; k < n_times; k++)
            {
                for (int l = 0; l < num_fields; l++)
                {
                    observation_data[num_fields*n_times*n_scen*ii + num_fields*n_times*jj + num_fields*k + l] = rdata[num_fields*k + l];
                }
            }

            // Get the ignition dataset name
            strcpy(dataset_name_ignition, group_name);
            strcat(dataset_name_ignition, "/ignition");

            // Now open the ignition dataset
            dataset_id_ignition = H5Dopen(file_id, dataset_name_ignition, H5P_DEFAULT);

            // Read ignition dataset
            status = H5Dread(dataset_id_ignition, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                             idata);

            // Now put idata into ignition_data
            ignition_data[2*scen]     = idata[0];
            ignition_data[2*scen + 1] = idata[1];

            // Open attribute dataset
            attribute_id = H5Aopen(dataset_id_truth, attribute_name, H5P_DEFAULT);

            // Now read in the attributes
            status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, adata);

            // Now put adata into scenario_params
            scenario_params[2*n_scen*ii + 2*jj    ] = adata[0];
            scenario_params[2*n_scen*ii + 2*jj + 1] = adata[1];

            // Close time dataset
            status += H5Dclose(dataset_id_time);

            // Close attribute dataset
            status += H5Aclose(attribute_id);

            // Close truth dataset
            status += H5Dclose(dataset_id_truth);

            // Close ignition dataset
            status += H5Dclose(dataset_id_ignition);

        }
    }

    // Close file
    status += H5Fclose(file_id);

}
