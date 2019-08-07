module cR2S_modules_api_m

use iso_c_binding
use cR2S_modules
use utils_api_m

implicit none
private

public :: c_read_common_file
public :: c_mesh_selection
public :: c_voxel_selection
public :: c_rec_pos_sample
public :: c_cyl_pos_sample
public :: c_energy_selection
public :: c_cell_adjust_weight
public :: c_check_position
public :: c_transform
public :: c_get_meshtype
public :: c_get_meshstruc
public :: c_get_meshtransform
public :: c_set_meshtransform
public :: c_get_cellid
public :: c_get_cellnum
public :: get_c_ptr


! Type to wrap the src
type :: source_wrapper_t
    type(source_attributes), allocatable :: src(:)
end type source_wrapper_t

type(source_wrapper_t), target :: STATIC_SOURCETYPE

integer(c_int), public, bind(C, name="MCR2S_MAX_ERROR_LENGTH") :: MCR2S_MAX_ERROR_LENGTH = MAX_ERROR_LENGTH+1

contains


    ! Function to get the c location of the src data type
    function get_c_ptr() result(ptr) bind(c, name='get_c_ptr')
        type(c_ptr) :: ptr
        
        ptr = c_loc(STATIC_SOURCETYPE)
        
    end function get_c_ptr
    
    function get_logical(logical_int) result(logical_return)
        logical :: logical_return
        integer(c_int) :: logical_int
        
        ! Create logical
        if(logical_int /=0) then
            logical_return = .True.
        else
            logical_return = .False.
        endif
        
    end function get_logical
        
   
    
    ! Subroutine to read the commonfile
    subroutine c_read_common_file(c_common_filename, src_pointer, nmeshes, total_src, verbose_test, c_error)&
                                  bind(c,name='c_read_common_file')
    
        character(kind=c_char), dimension(*), intent(in) :: c_common_filename
        type(c_ptr), value, intent(in) :: src_pointer
        integer(c_int), intent(out) :: nmeshes
        real(c_double), intent(out) :: total_src
        integer(c_int), intent(in), value :: verbose_test

        character(kind=c_char), intent(inout) :: c_error(MCR2S_MAX_ERROR_LENGTH)
        
        logical :: verbose
        character (len=MCR2S_MAX_ERROR_LENGTH) :: error
        type(source_wrapper_t), pointer :: src_wrapper_ptr

        error=""
                   
        ! Turn c pointer into fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)
        
        ! Create logical
        verbose = get_logical(verbose_test)

        call read_common_file(c_f_string(c_common_filename, MCR2S_MAX_ERROR_LENGTH), src_wrapper_ptr%src, &
                              nmeshes, total_src, verbose, error)


        c_error = f_c_string(error)

    end subroutine c_read_common_file

    ! Subroutine to select a mesh
    subroutine c_mesh_selection(src_pointer, sample_type, random_num, total_src, mesh_idx, wgt, verbose_int)&
                                bind(c,name='c_mesh_selection')

        type(c_ptr), intent(in), value      :: src_pointer 
        integer(c_int), intent(in), value   :: sample_type   ! Type of sampling mode (see note below)
        integer(c_int), intent(in), value   :: verbose_int   ! Write debug info to screen
        real(c_double), intent(in), value   :: random_num    ! Random number between 0-1 used to pick mesh
        real(c_double), intent(in), value   :: total_src     ! Total source strength
    
        integer(c_int), intent(out)         :: mesh_idx      ! Picked mesh IDX value
        real(c_double), intent(inout)       :: wgt           ! Weight of the starting particle for the picked mesh         
    
        type(source_wrapper_t), pointer :: src_wrapper_ptr
        logical :: verbose
        
        
        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)
        

        
        verbose = get_logical(verbose_int)
        call mesh_selection(src_wrapper_ptr%src, sample_type, random_num, total_src, mesh_idx, wgt, verbose)
        
   
    end subroutine c_mesh_selection
    
    ! Subroutine to select voxel
    subroutine c_voxel_selection(src_pointer, mesh_idx, sample_type, random_num, selected_voxel,&
                                selected_mesh_element, i, j, k, wgt, verbose_int)               &
                                bind(c,name='c_voxel_selection')
                                
        type(c_ptr), value, intent(in)      :: src_pointer           ! Chosen source mesh
        integer(c_int), value, intent(in)          :: sample_type           ! Type of sampling mode (sample_type = 1 analogue sampling, sample_type /= 1 uniform sampling)
        real(c_double), value, intent(in)          :: random_num            ! Random number between 0-1 used to pick voxel
        integer(c_int), value, intent(in)          :: verbose_int           ! Write debug info to screen    
        integer(c_int), value, intent(in)          :: mesh_idx
        
        integer(c_int), intent(out)         :: selected_voxel        ! Selected voxel ID
        integer(c_int), intent(out)         :: selected_mesh_element ! Reference to the selected element in the mesh voxel array 
        integer(c_int), intent(out)         :: i, j, k               ! x, y and z index of the voxel

        real(c_double), intent(inout)       :: wgt                   ! Weight of the starting particle    
    
        type(source_wrapper_t), pointer :: src_wrapper_ptr
        logical :: verbose
        
        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)
        
        verbose = get_logical(verbose_int)
        call voxel_selection(src_wrapper_ptr%src(mesh_idx), sample_type, random_num, selected_voxel, &
                             selected_mesh_element, i, j, k, wgt, verbose)
    
    end subroutine c_voxel_selection

    
    ! Pick a position in a rectangular voxel
    subroutine c_rec_pos_sample(src_pointer, mesh_idx, i, j, k, random_nums, x, y, z) bind(c,name='c_rec_pos_sample')
    
        type(c_ptr), value, intent(in)      :: src_pointer           ! Chosen source mesh
        integer(c_int), value, intent(in)   :: mesh_idx
        integer(c_int), value, intent(in)   :: i, j, k              ! x, y and z index of the voxel
        real(c_double), intent(in)          :: random_nums(1:3)     ! Random number between 0-1 used to pick voxel
        real(c_double), intent(out)         :: x, y, z              ! X, Y, Z position in voxel

        type(source_wrapper_t), pointer :: src_wrapper_ptr        

        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)
        
        call rec_pos_sample(src_wrapper_ptr%src(mesh_idx),i, j, k, random_nums, x, y, z)

    end subroutine c_rec_pos_sample

    
    ! Pick a position in a rectangular voxel
    subroutine c_cyl_pos_sample(src_pointer, mesh_idx, i, j, k, random_nums, x, y, z) bind(c,name='c_cyl_pos_sample')
    
        type(c_ptr), value, intent(in)      :: src_pointer          ! Chosen source mesh
        integer(c_int), value, intent(in)   :: mesh_idx
        integer(c_int), value, intent(in)   :: i, j, k              ! x, y and z index of the voxel
        real(c_double), intent(in)          :: random_nums(1:3)     ! Random number between 0-1 used to pick voxel
        real(c_double), intent(out)         :: x, y, z              ! X, Y, Z position in voxel

        type(source_wrapper_t), pointer :: src_wrapper_ptr        
        
        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)
        
        call cyl_pos_sample(src_wrapper_ptr%src(mesh_idx),i, j, k, random_nums, x, y, z)

    end subroutine c_cyl_pos_sample
    
    ! Sample the particle energy 
    subroutine c_energy_selection(src_pointer, mesh_idx, selected_mesh_element, selected_cell, energy_sample_type, &
                                  random_nums, ignore_cell_info_int, erg, wgt)                                                &  
                                  bind(c,name='c_energy_selection')
                                  
        type(c_ptr), value, intent(in)      :: src_pointer           ! Chosen source mesh
        integer(c_int), value, intent(in)   :: mesh_idx        
        integer(c_int), value, intent(in)   :: selected_mesh_element ! Selected voxel ID
        integer(c_int), intent(inout)       :: selected_cell         ! Selected voxel Cell
        integer(c_int), value, intent(in)   :: energy_sample_type    ! Type of sampling mode (energy_sample_type = 0 analogue sampling, energy_sample_type /= 1 uniform sampling)
        real(c_double), intent(in)          :: random_nums(1:2)      ! Random number between 0-1 used to pick energy
        integer(c_int), value, intent(in)   :: ignore_cell_info_int  ! Flag to say to ignore the cell info (True: ignore_cell_info_int = 1)
        real(c_double), intent(out)         :: erg                   ! Energy of the starting particle
        real(c_double), intent(inout)       :: wgt                   ! Weight of the starting particle 
        
        type(source_wrapper_t), pointer :: src_wrapper_ptr
        logical :: ignore_cell_info

        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)

        ignore_cell_info = get_logical(ignore_cell_info_int)
        
        call energy_selection(src_wrapper_ptr%src(mesh_idx)%voxel(selected_mesh_element), selected_cell, &
                              src_wrapper_ptr%src(mesh_idx)%mesh%ebin, energy_sample_type, random_nums,  &
                              ignore_cell_info, erg, wgt)
        
    end subroutine c_energy_selection
    
    ! Subroutine to adjust the weight if cell under voxel rejection is used
    subroutine c_cell_adjust_weight(src_pointer, mesh_idx, selected_mesh_element, selected_cell, wgt) &
                                    bind(c,name='c_cell_adjust_weight')
                                  
        type(c_ptr), value, intent(in)      :: src_pointer           ! Chosen source mesh                                  
        integer(c_int), value, intent(in)   :: mesh_idx        
        integer(c_int), value, intent(in)   :: selected_mesh_element ! Selected voxel ID
        integer(c_int), value, intent(in)   :: selected_cell         ! Selected voxel Cell        
        real(c_double), intent(inout)       :: wgt                   ! Weight of the starting particle 
        
        type(source_wrapper_t), pointer :: src_wrapper_ptr

        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)

        call cell_adjust_weight(src_wrapper_ptr%src(mesh_idx)%voxel(selected_mesh_element), selected_cell, wgt)       
        
    end subroutine c_cell_adjust_weight
    
    ! Subroutine to check the selected spatial point 
    subroutine c_check_position(src_pointer, mesh_idx, selected_mesh_element, picked_cell, cell_importance, material, &
                                ignore_cell_info_int, resample_in_voxel_int, selected_cell, verbose_int)              &
                                bind(c,name='c_check_position') 

        type(c_ptr), value, intent(in)      :: src_pointer           ! Chosen source mesh                                  
        integer(c_int), value, intent(in)   :: mesh_idx        
        integer(c_int), value, intent(in)   :: selected_mesh_element ! Selected voxel ID  
        integer(c_int), value, intent(in)   :: picked_cell           ! Selected Cell ID
        integer(c_int), intent(out)         :: selected_cell     
        real(c_double), value, intent(in)   :: cell_importance
        integer(c_int), value, intent(in)   :: material         
        integer(c_int), value, intent(in)   :: ignore_cell_info_int  ! Flag to say to ignore the cell info (True: ignore_cell_info_int = 1)
        integer(c_int), value, intent(in)   :: verbose_int
        integer(c_int), intent(out)         :: resample_in_voxel_int
         
        type(source_wrapper_t), pointer :: src_wrapper_ptr         
        logical :: verbose
        logical :: resample_in_voxel
        logical :: ignore_cell_info

        resample_in_voxel_int = 0 
        
        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)

        verbose = get_logical(verbose_int)
        resample_in_voxel = get_logical(resample_in_voxel_int)
        
        call check_position(src_wrapper_ptr%src, mesh_idx, selected_mesh_element, picked_cell, cell_importance, material, &
                            ignore_cell_info, resample_in_voxel, selected_cell, verbose)
        
        if(resample_in_voxel) then
            resample_in_voxel_int = 1
        else
            resample_in_voxel_int = 0
        end if
        
    end subroutine c_check_position     

    ! Subroutine to apply a source rotation and translation if applied
    subroutine c_transform(x,y,z,u,v,w,cosines,translations,verbose_int) bind(c,name='c_transform') 
    
        real(c_double), intent(in) :: cosines(1:9)
        real(c_double), intent(in) :: translations(1:3)
        real(c_double), intent(inout) :: x, y, z
        real(c_double), intent(inout) :: u, v, w
        integer(c_int), value, intent(in)    :: verbose_int
        
        logical :: verbose
        
        verbose = get_logical(verbose_int)
        
        call transform(x,y,z,u,v,w,cosines,translations,verbose)
        
        
    end subroutine c_transform
    
    ! Get the type of mesh
    function c_get_meshtype(src_pointer, mesh_idx) result(mesh_type) bind(c,name='c_get_meshtype') 
    
        type(c_ptr), value, intent(in)   :: src_pointer
        integer, value, intent(in)       :: mesh_idx
        integer                          :: mesh_type     
        
        type(source_wrapper_t), pointer :: src_wrapper_ptr      
        
        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)

        mesh_type = src_wrapper_ptr%src(mesh_idx)%mesh%geom_id
    
    end function c_get_meshtype
    
    ! Get the structure of mesh
    function c_get_meshstruc(src_pointer, mesh_idx) result(struc_type) bind(c,name='c_get_meshstruc') 
    
        type(c_ptr), value, intent(in)   :: src_pointer
        integer, value, intent(in)       :: mesh_idx
        integer                          :: struc_type     
        
        type(source_wrapper_t), pointer :: src_wrapper_ptr      
        
        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)

        struc_type = src_wrapper_ptr%src(mesh_idx)%mesh%struc_ID
    
    end function c_get_meshstruc
    
    ! Get the transform of the mesh
    function c_get_meshtransform(src_pointer, mesh_idx) result(trans_num) bind(c,name='c_get_meshtransform') 
    
        type(c_ptr), value, intent(in)   :: src_pointer
        integer, value, intent(in)       :: mesh_idx
        integer                          :: trans_num     
        
        type(source_wrapper_t), pointer :: src_wrapper_ptr      
        
        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)

        trans_num = src_wrapper_ptr%src(mesh_idx)%mesh%trans_num
    
    end function c_get_meshtransform
    
    ! Set the transform of the mesh
    subroutine c_set_meshtransform(src_pointer, mesh_idx, trans_num) bind(c,name='c_set_meshtransform') 
    
        type(c_ptr), value, intent(in)   :: src_pointer
        integer, value, intent(in)       :: mesh_idx
        integer, value, intent(in)       :: trans_num     
        
        type(source_wrapper_t), pointer :: src_wrapper_ptr      
        
        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)

        src_wrapper_ptr%src(mesh_idx)%mesh%trans_num = trans_num
    
    end subroutine c_set_meshtransform
    
    ! Get the cell id
    function c_get_cellid(src_pointer, mesh_idx, selected_mesh_element, selected_cell) result(cellid) bind(c,name='c_get_cellid') 
    
        type(c_ptr), value, intent(in)   :: src_pointer
        integer, value, intent(in)       :: mesh_idx
        integer, value, intent(in)       :: selected_mesh_element
        integer, value, intent(in)       :: selected_cell
        integer                          :: cellid     
        
        type(source_wrapper_t), pointer :: src_wrapper_ptr      
        
        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)

        cellid = src_wrapper_ptr%src(mesh_idx)%voxel(selected_mesh_element)%cel(selected_cell)%id
    
    end function c_get_cellid
    
    ! Get the number of cells in voxel
    function c_get_cellnum(src_pointer, mesh_idx, selected_mesh_element) result(cellnum) bind(c,name='c_get_cellnum') 
    
        type(c_ptr), value, intent(in)   :: src_pointer
        integer, value, intent(in)       :: mesh_idx
        integer, value, intent(in)       :: selected_mesh_element
        integer                          :: cellnum     
        
        type(source_wrapper_t), pointer :: src_wrapper_ptr      
        
        ! Turn the c pointer into a fortran native pointer
        call c_f_pointer(src_pointer, src_wrapper_ptr)

        cellnum = src_wrapper_ptr%src(mesh_idx)%voxel(selected_mesh_element)%ncel
    
    end function c_get_cellnum
    
    
end module cR2S_modules_api_m
