#include <stdio.h>
#include <iostream>
#include <ctime>
#include <vector>
#include <array>
#include <typeinfo>
#include <fstream>

#include "openmc/bank.h"
#include "openmc/constants.h"
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/settings.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
#include "openmc/error.h"
#include "openmc/constants.h"
#include "openmc/locations.h"
#include "../vendor/cr2s/cR2Smodules.h"
#include "openmc/cR2Sheader.h"
#include "openmc/particle.h"
#include "openmc/cR2Sarray.h"
#include "openmc/message_passing.h"
#include "openmc/position.h"

#define FUNCTION_NAME "MCR2SSrc:"
namespace openmc{
//***************************************************************************//

double * get_random_numbers(int, long);

static void* src_pointer;
static int spatial_sampling = 0;
static int energy_sampling = 0;
static int remax;
static int nmeshes;
static int verbose_int = 0;
static int ignore_cell_int = 0;
static double total_src;


// Read source file //
void ReadCDGS()
{

  char c_error[MCR2S_MAX_ERROR_LENGTH];


  src_pointer = get_c_ptr();

  std::cout << "cR2S Source _reading data" << std::endl;

  // Read the common file //
  c_read_common_file("commonfile", src_pointer, &nmeshes, &total_src, verbose_int, c_error);
  if(c_error[0] != '\0')
  {
      fatal_error(c_error);
  }

  // Find the correct transform number //
  // Transforms not currently implemented//

  std::cout << "cR2S Source _ finished reading data" << std::endl;

}

void MCR2SSrc(long src, double *x, double *y, double *z, double *u,  double *v, 
	          double *w, double *E, double *wgt, double *t, long id)

{  
  long ptr, np;
  double x0, y0, z0, xmin, xmax, ymin, ymax, zmin, zmax, u0, v0, w0, E0, rnd, pitch;
  const double *params;

  //*************************************************************************//


  //*************************************************************************//

  // Subroutine variables //
  int selected_voxel;
  int selected_mesh_element;
  int i, j, k;
  int mesh_idx;
  int resample_in_voxel;
  int resamc;
  int selected_cell;
  long cell;
  long mat;
  double * random_numbers;
  long ID;
  char dummyString[4*MAX_STR];
  int print_source_particles;



  // Set values from user parameters //

  /*   spatial sampling mode
              0 = uniform sampling of active voxels.
                  Weight correction for voxel selection and cell under voxel selection <default>
              1 = analogue sampling of active voxels
                  Weight correction for cell under voxel selection.
                  cell under mesh selection. */
  spatial_sampling = int(settings::mcr2s_spatial_sampling);

  /* energy sampling mode
               0 = analogue sampling of energy bins.
               1 = uniform sampling of active energy bins */
  energy_sampling = int(settings::mcr2s_energy_sampling);

  /* If > 0, ignore cell under mesh cell numbers (e.g. dummy cell data where cell=0). */
  ignore_cell_int = int(IGNORE_CELL_UNDER_MESH);

  /* number of spatial retries before abandoning voxel [optional, default = 1000] */
  remax = int(settings::mcr2s_remax);

  verbose_int = 0;
  print_source_particles = 1;

  prn_set_stream(STREAM_SOURCE);                // Left untouched //
  rnd = prn();                                  // Initialised random number //



  *wgt = 1.0;
  cell = 0;
  mat = 0;
  mesh_idx = 0;
  random_numbers = NULL;

  c_mesh_selection(src_pointer, spatial_sampling, rnd, total_src, &mesh_idx, wgt, verbose_int);

 // Restart point, if re-sampling can't find acceptable cell within limit of attemps //

  for (;;)
  {

    prn_set_stream(STREAM_SOURCE);
    rnd = prn();

    c_voxel_selection(src_pointer, mesh_idx, spatial_sampling, rnd, &selected_voxel, &selected_mesh_element, &i, &j, &k, wgt, verbose_int);

    if (verbose_int == 1) printf("Active voxel index : %d\n",selected_mesh_element);
    if (verbose_int == 1) printf("Corresponding element : %d\n",selected_voxel );

    for (resamc=0;resamc<remax;resamc++)
    {   
        random_numbers = get_random_numbers(3, ID);
        if(c_get_meshtype(src_pointer,mesh_idx) == 1)
        {
            c_rec_pos_sample(src_pointer, mesh_idx, i,  j,  k, random_numbers, x, y, z);
        }
        else if (c_get_meshtype(src_pointer,mesh_idx) == 2)
        {
            c_cyl_pos_sample(src_pointer, mesh_idx, i,  j,  k, random_numbers, x, y, z);
        }    
        else if (c_get_meshstruc(src_pointer,mesh_idx) == 2)
        {
            fatal_error(std::string("Unstructure usersrc not supported yet."));
        }        
        free(random_numbers);


        if (verbose_int == 1) printf("Selected point : %f, %f, %f\n", *x, *y, *z);

        bool found = false;
        int n_reject = 0;
        static int n_accept = 0;

        double xyz[] {*x,*y,*z};

        int32_t cell_index, instance;
        openmc_find_cell(xyz, &cell_index, &instance);
        cell = cell_index;
        const auto& c = model::cells[cell_index];
        auto mat_index = c->material_.size() == 1
        ? c->material_[0] : c->material_[instance];
        if (mat_index == MATERIAL_VOID) {
            mat = 0; // If material not found, assume to be void = 0 //
          }
        else {
            mat = mat_index;
        }
        c_check_position(src_pointer, mesh_idx, selected_mesh_element, cell, (double)1.0, (int)mat, ignore_cell_int, &resample_in_voxel, &selected_cell, verbose_int);

        if (resample_in_voxel == 0) 
        {
            break;
        }

    }

    //*** If failed to find a material and cell, restart, i.e. reached max count ***//
    if (resamc == remax)
    {
      continue;
    }  

    if (c_get_cellid(src_pointer, mesh_idx, selected_mesh_element, selected_cell) > 0 &&
        c_get_cellnum(src_pointer, mesh_idx, selected_mesh_element) > 1 &&
        ignore_cell_int == 0)
    {
        c_cell_adjust_weight(src_pointer, mesh_idx, selected_mesh_element, selected_cell, wgt);
    }

    if (verbose_int == 1) printf("Sampling energy\n");


    random_numbers = get_random_numbers(2,ID);

    c_energy_selection(src_pointer, mesh_idx, selected_mesh_element, &selected_cell, 
                       energy_sampling, random_numbers, ignore_cell_int, E, wgt);
    free(random_numbers);


    if (print_source_particles == 1)
    {
	    printf("%d,%f,%f,%f,%ld,%f,%f\n", selected_voxel,*x,*y,*z,cell,*E,*wgt);
	    std::ofstream ofile;
	    ofile.open("SourceParticles.txt", std::ios::app);
	    ofile << selected_voxel << " " << *x << " " << *y << " " << *z << " " << cell << " " << *E << " " << wgt << std::endl;
	    ofile.close();
    }
    break;
  }
}

double * get_random_numbers(int quant, long ID)
{

    int i;
    double * random_numbers;
    prn_set_stream(STREAM_SOURCE);
    random_numbers = (double*)calloc(quant,sizeof(double));

    for (i=0;i<quant;i++)
    {
        random_numbers[i] = prn();
    }

    return random_numbers;
}
}
