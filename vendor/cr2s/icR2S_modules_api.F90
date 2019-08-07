module icR2S_modules_api_m
use utils_api_m
use iso_c_binding
implicit none
public

    integer(c_int), parameter :: MCR2S_MAX_ERROR_LENGTH = MAX_ERROR_LENGTH + 1     

    interface

        function get_c_ptr() result(ptr) bind(C, name="get_c_ptr")
            use iso_c_binding
            implicit none
            type(c_ptr) :: ptr
        end function get_c_ptr

        subroutine c_read_common_file(c_common_filename, src_pointer, nmeshes, total_src, verbose_test, c_error)&
                                      bind(c,name='c_read_common_file')
        
            use iso_c_binding
            implicit none
            
            type(c_ptr), intent(in), value :: src_pointer
            integer(c_int), intent(out) :: nmeshes
            real(c_double), intent(out) :: total_src
            integer(c_int), intent(in), value :: verbose_test
            character(kind=c_char), dimension(*), intent(in) :: c_common_filename
            character(kind=c_char), dimension(*), intent(inout) :: c_error       
        
        end subroutine c_read_common_file
        
        
        subroutine c_mesh_selection(src_pointer, sample_type, random_num, total_src, mesh_idx, wgt, verbose_int)&
                                    bind(c,name='c_mesh_selection')

            use iso_c_binding
            implicit none
            
            type(c_ptr), intent(in), value      :: src_pointer 
            integer(c_int), intent(in), value   :: sample_type   ! Type of sampling mode (see note below)
            integer(c_int), intent(in), value   :: verbose_int   ! Write debug info to screen
            real(c_double), intent(in), value   :: random_num    ! Random number between 0-1 used to pick mesh
            real(c_double), intent(in), value   :: total_src     ! Total source strength
        
            integer(c_int), intent(out)         :: mesh_idx      ! Picked mesh IDX value
            real(c_double), intent(inout)       :: wgt           ! Weight of the starting particle for the picked mesh            
    
        end subroutine c_mesh_selection
     
     
        ! Subroutine to select voxel
        subroutine c_voxel_selection(src_pointer, mesh_idx, sample_type, random_num, selected_voxel, &
                                    selected_mesh_element, i, j, k, wgt, verbose_int) bind(c,name='c_voxel_selection')
        
            use iso_c_binding
            implicit none
            
            type(c_ptr), value, intent(in)      :: src_pointer           ! Chosen source mesh
            integer(c_int), value, intent(in)   :: sample_type           ! Type of sampling mode (sample_type = 1 analogue sampling, sample_type /= 1 uniform sampling)
            real(c_double), value, intent(in)   :: random_num            ! Random number between 0-1 used to pick voxel
            integer(c_int), value, intent(in)   :: verbose_int           ! Write debug info to screen    
            integer(c_int), value, intent(in)   :: mesh_idx
            
            integer(c_int), intent(out)         :: selected_voxel        ! Selected voxel ID
            integer(c_int), intent(out)         :: selected_mesh_element ! Reference to the selected element in the mesh voxel array 
            integer(c_int), intent(out)         :: i, j, k               ! x, y and z index of the voxel
    
            real(c_double), intent(inout)       :: wgt                   ! Weight of the starting particle    
        
        end subroutine c_voxel_selection
        
        
        ! Pick a position in a rectangular voxel
        subroutine c_rec_pos_sample(src_pointer, mesh_idx, i, j, k, random_nums, x, y, z) bind(c,name='c_rec_pos_sample')

            use iso_c_binding
            implicit none
       
            type(c_ptr), value, intent(in)      :: src_pointer          ! Chosen source mesh
            integer(c_int), value, intent(in)   :: mesh_idx
            integer(c_int), value, intent(in)   :: i, j, k              ! x, y and z index of the voxel
            real(c_double), intent(in)          :: random_nums(1:3)     ! Random number between 0-1 used to pick voxel
            real(c_double), intent(out)         :: x, y, z              ! X, Y, Z position in voxel
    
    
        end subroutine c_rec_pos_sample
        
        ! Pick a position in a rectangular voxel
        subroutine c_cyl_pos_sample(src_pointer, mesh_idx, i, j, k, random_nums, x, y, z) bind(c,name='c_cyl_pos_sample')

            use iso_c_binding
            implicit none
        
            type(c_ptr), value, intent(in)      :: src_pointer           ! Chosen source mesh
            integer(c_int), value, intent(in)   :: mesh_idx
            integer(c_int), value, intent(in)   :: i, j, k               ! x, y and z index of the voxel
            real(c_double), intent(in)          :: random_nums(1:3)        ! Random number between 0-1 used to pick voxel
            real(c_double), intent(out)         :: x, y, z              ! X, Y, Z position in voxel
    
    
        end subroutine c_cyl_pos_sample
        
        ! Sample the particle energy 
        subroutine c_energy_selection(src_pointer, mesh_idx, selected_mesh_element, selected_cell, energy_sample_type, &
                                    random_nums, ignore_cell_info_int, erg, wgt)                                                &  
                                    bind(c,name='c_energy_selection')

            use iso_c_binding
            implicit none
                                    
            type(c_ptr), value, intent(in)      :: src_pointer           ! Chosen source mesh
            integer(c_int), value, intent(in)   :: mesh_idx        
            integer(c_int), value, intent(in)   :: selected_mesh_element ! Selected voxel ID
            integer(c_int), intent(inout)       :: selected_cell         ! Selected voxel Cell
            integer(c_int), value, intent(in)   :: energy_sample_type    ! Type of sampling mode (energy_sample_type = 0 analogue sampling, energy_sample_type /= 1 uniform sampling)
            real(c_double), intent(in)          :: random_nums(1:2)        ! Two random number between 0-1 used to pick energy
            integer(c_int), value, intent(in)   :: ignore_cell_info_int  ! Flag to say to ignore the cell info (True: ignore_cell_info_int = 1)
            real(c_double), intent(out)         :: erg                   ! Energy of the starting particle
            real(c_double), intent(inout)       :: wgt                   ! Weight of the starting particle 
            
        end subroutine c_energy_selection
        
        ! Subroutine to adjust the weight if cell under voxel rejection is used
        subroutine c_cell_adjust_weight(src_pointer, mesh_idx, selected_mesh_element, selected_cell, wgt) &
                                        bind(c,name='c_cell_adjust_weight')

            use iso_c_binding
            implicit none
                                        
            type(c_ptr), value, intent(in)      :: src_pointer           ! Chosen source mesh                                  
            integer(c_int), value, intent(in)   :: mesh_idx        
            integer(c_int), value, intent(in)   :: selected_mesh_element ! Selected voxel ID
            integer(c_int), value, intent(in)   :: selected_cell         ! Selected voxel Cell        
            real(c_double), intent(inout)       :: wgt                   ! Weight of the starting particle 
            
        end subroutine c_cell_adjust_weight
        
        ! Subroutine to check the selected spatial point 
        subroutine c_check_position(src_pointer, mesh_idx, selected_mesh_element, picked_cell, cell_importance, material, &
                                    ignore_cell_info_int, resample_in_voxel_int, selected_cell, verbose_int)              &
                                    bind(c,name='c_check_position') 
    
            use iso_c_binding
            implicit none
                                
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
            
        end subroutine c_check_position 
        
        ! Subroutine to apply a source rotation and translation if applied
        subroutine c_transform(x,y,z,u,v,w,cosines,translations,verbose_int) bind(c,name='c_transform') 

            use iso_c_binding
            implicit none
        
            real(c_double), intent(in) :: cosines(1:9)
            real(c_double), intent(in) :: translations(1:3)
            real(c_double), intent(inout) :: x, y, z
            real(c_double), intent(inout) :: u, v, w
            integer(c_int), value, intent(in)    :: verbose_int
            
        end subroutine c_transform
        
        ! Get the type of mesh
        function c_get_meshtype(src_pointer, mesh_idx) result(mesh_type) bind(c,name='c_get_meshtype') 

            use iso_c_binding
            implicit none
        
            type(c_ptr), value, intent(in)   :: src_pointer
            integer, value, intent(in)       :: mesh_idx
            integer                          :: mesh_type
     
        end function c_get_meshtype
        
        ! Get the structure of mesh
        function c_get_meshstruc(src_pointer, mesh_idx) result(struc_type) bind(c,name='c_get_meshstruc') 

            use iso_c_binding
            implicit none
        
            type(c_ptr), value, intent(in)   :: src_pointer
            integer, value, intent(in)       :: mesh_idx
            integer                          :: struc_type     
        
        end function c_get_meshstruc
     
        ! Get the transform of the mesh
        function c_get_meshtransform(src_pointer, mesh_idx) result(trans_num) bind(c,name='c_get_meshtransform') 
        
            use iso_c_binding
            implicit none        
        
            type(c_ptr), value, intent(in)   :: src_pointer
            integer, value, intent(in)       :: mesh_idx
            integer                          :: trans_num     
        
        end function c_get_meshtransform  

        ! Set the transform of the mesh
        subroutine c_set_meshtransform(src_pointer, mesh_idx, trans_num) bind(c,name='c_set_meshtransform') 
        
            use iso_c_binding
            implicit none  
        
            type(c_ptr), value, intent(in)   :: src_pointer
            integer, value, intent(in)       :: mesh_idx
            integer, value, intent(in)       :: trans_num     
        
        end subroutine c_set_meshtransform        
        
        ! Get the cell id
        function c_get_cellid(src_pointer, mesh_idx, selected_mesh_element, selected_cell) result(cellid)&
                              bind(c,name='c_get_cellid') 
        
            use iso_c_binding
            implicit none
        
            type(c_ptr), value, intent(in)   :: src_pointer
            integer, value, intent(in)       :: mesh_idx
            integer, value, intent(in)       :: selected_mesh_element
            integer, value, intent(in)       :: selected_cell
            integer                          :: cellid     
        
        end function c_get_cellid
        
        ! Get the number of cells in voxel
        function c_get_cellnum(src_pointer, mesh_idx, selected_mesh_element) result(cellnum) bind(c,name='c_get_cellnum') 
        
            use iso_c_binding
            implicit none
        
            type(c_ptr), value, intent(in)   :: src_pointer
            integer, value, intent(in)       :: mesh_idx
            integer, value, intent(in)       :: selected_mesh_element
            integer                          :: cellnum     
        
        end function c_get_cellnum
        
     
    end interface

end module icR2S_modules_api_m
