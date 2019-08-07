#include <stdio.h>
#include <iostream>
#include <ctime>
#include <vector>
#include <array>

#include "openmc/bank.h"
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
#include "/home/mmavis/openmc/vendor/cr2s/cR2Smodules.h"
#include "openmc/cR2Sheader.h"
#include "openmc/particle.h"
#include "openmc/cR2Sarray.h"
#include "openmc/message_passing.h"

#define FUNCTION_NAME "MCR2SSrc:"
namespace openmc{
//***************************************************************************//

double * get_random_numbers(int, long);

static void* src_pointer;
static int spatial_sampling = 0;
static int energy_sampling = 0;
static int remax;
static int nmeshes;
static int verbose_int = 1;
static int ignore_cell_int = 0;
static double total_src;





//*** Params explanation ***//
/*
  params[0]  = spatial sampling mode
              0 = uniform sampling of active voxels.
                  Weight correction for voxel selection and cell under voxel selection <default>
              1 = analogue sampling of active voxels
                  Weight correction for cell under voxel selection.
                  cell under mesh selection.
  params[1]  = energy sampling mode
               0 = analogue sampling of energy bins.
               1 = uniform sampling of active energy bins
  params[2] = If > 0, ignore cell under mesh cell numbers (e.g. dummy cell data where cell=0).
  params[3]  = number of spatial retries before abandoning voxel [optional, default = 1000]
*/  


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
  std::cout << src_pointer << std::endl;
  std::cout << &nmeshes << std::endl;
  std::cout << &total_src << std::endl;
  std::cout << verbose_int << std::endl;
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
  
  
  std::cout << "Set up user parameters" << std::endl;
  // Set values from user parameters //
  //spatial_sampling = (int)params[0];
  spatial_sampling = 0;
  //energy_sampling = (int)params[1];
  energy_sampling = 0;
  //ignore_cell_int = (int)params[2];
  ignore_cell_int = 0;
  remax = 1000;

   /* if (params[3] > 0)
  {
      remax = (long)params[3];
  }
  else
  {
      remax = 1000;
  } */
 
  //verbose_int = (int)params[4]  ;
  verbose_int = 1  ;
  
  std::cout << "Setting random number" << std::endl;
  prn_set_stream(STREAM_VOLUME);                // Left untouched //
  rnd = prn();                                  // Initialised random number //
  


  std::cout << "Setting initial variables" << std::endl;
  *wgt = 1.0;
  cell = 0;
  mat = 0;
  mesh_idx = 0;
  random_numbers = NULL;
  
  std::cout << "Starting mesh_selection" << std::endl;
  c_mesh_selection(src_pointer, spatial_sampling, rnd, total_src, &mesh_idx, wgt, verbose_int);

 // Restart point, if re-sampling can't find acceptable cell within limit of attemps //
 
  for (;;)
  {
      
    prn_set_stream(STREAM_SOURCE);
    rnd = prn();
    std::cout << rnd << std::endl;
    c_voxel_selection(src_pointer, mesh_idx, spatial_sampling, rnd, &selected_voxel, &selected_mesh_element, &i, &j, &k, wgt, verbose_int);
  
    if (verbose_int == 1) printf("Active voxel index : %d\n",selected_mesh_element);
    if (verbose_int == 1) printf("Corresponding element : %d\n",selected_voxel );

    for (resamc=0;resamc<remax;resamc++)
    {   
        std::cout << "Get random Number" << std::endl;
        random_numbers = get_random_numbers(3, ID);
        std::cout << "Mesh type" << std::endl;
        if(c_get_meshtype(src_pointer,mesh_idx) == 1)
        {
            std::cout << "Rec" << std::endl;
            c_rec_pos_sample(src_pointer, mesh_idx, i,  j,  k, random_numbers, x, y, z);
            std::cout << "Sampled Rec" << std::endl;
        }
        else if (c_get_meshtype(src_pointer,mesh_idx) == 2)
        {
            std::cout << "Cyl" << std::endl;
            c_cyl_pos_sample(src_pointer, mesh_idx, i,  j,  k, random_numbers, x, y, z);
            std::cout << "Sampled Cyl" << std::endl;
        }    
        else if (c_get_meshstruc(src_pointer,mesh_idx) == 2)
        {
            fatal_error(std::string("Unstructure usersrc not supported yet."));
        }        
        free(random_numbers);
        
        
        if (verbose_int == 1) printf("Selected point : %f, %f, %f\n", *x, *y, *z);
        Particle p;
        Particle::Bank site;
        site.particle = Particle::Type::photon;
        settings::photon_transport = true;
        // Set weight to one by default
        site.wgt = 1.0;
        // Repeat sampling source location until a good site has been found
        bool found = false;
        int n_reject = 0;
        static int n_accept = 0;
       
        // Set particle type
        /* particle_ = Particle::Type::photon;
        settings::photon_transport = true; */
        double xpos,ypos,zpos;
        xpos = *x;
        ypos = *y;
        zpos = *z;
        site.r.x = xpos;
        site.r.y = ypos;
        site.r.z = zpos;
        
        std::cout << typeid(site.r).name() << std::endl;
        std::cout << typeid(site.r.x).name() << std::endl;
        std::cout << typeid(site.r.y).name() << std::endl;
        std::cout << typeid(site.r.z).name() << std::endl;
        // Sample spatial distribution
            
       /*  std::cout << *x << std::endl;
        std::cout << *y << std::endl;
        std::cout << *z << std::endl;
        double xpos = *x;
        double ypos = *y;
        double zpos = *z;
        std::cout << xpos << std::endl;
        std::cout << ypos << std::endl;
        std::cout << zpos << std::endl;
        const double xyz[] {xpos,ypos,zpos}; */
        std::cout << site.r.x << std::endl;
        double xyz[] {site.r.x, site.r.y, site.r.z};
    // Now search to see if location exists in geometry
        //Position pos {*x,*y,*z};
        //find_cell(pos,false);
        std::cout << xyz << std::endl;
        int32_t cell_index, instance;
        std::cout << "Find Cell" << std::endl;
        int err = openmc_find_cell(xyz, &cell_index, &instance); //seg faulting
        std::cout << "Cell Found" << std::endl;
        cell = cell_index;
        const auto& c = model::cells[cell_index];
        auto mat_index = c->material_.size() == 1
        ? c->material_[0] : c->material_[instance];
        if (mat_index = MATERIAL_VOID) {
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
   

    if (verbose_int == 1) printf("%d,%f,%f,%f,%ld,%f,%f\n", selected_voxel,*x,*y,*z,cell,*E,*wgt);  

    break;
  }
}

double * get_random_numbers(int quant, long ID)
{

    int i;
    double * random_numbers;
    prn_set_stream(STREAM_VOLUME);
    random_numbers = (double*)calloc(quant,sizeof(double));
    
    for (i=0;i<quant;i++)
    {
        random_numbers[i] = prn();
    }
    
    return random_numbers;
}
}

