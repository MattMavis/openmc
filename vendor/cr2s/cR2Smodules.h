#ifndef cR2Smodules
#define cR2Smodules

#ifdef __cplusplus
extern "C" {
#endif

    extern const int MCR2S_MAX_ERROR_LENGTH;

    void* get_c_ptr();

    void c_read_common_file(const char* c_common_filename, void* src_pointer, int* nmeshes, double* total_src, int verbose_int, char* c_error);

    void c_mesh_selection(void* src_pointer, int sample_type, double random_num, double total_src, int* mesh_idx, double* wgt, int verbose_int);

    void c_voxel_selection(void* src_pointer, int mesh_idx, int sample_type, double random_num, int* selected_voxel, int* selected_mesh_element, int* i, int* j, int* k, double* wgt, int verbose_int);

    void c_rec_pos_sample(void* src_pointer, int mesh_idx, int i, int j, int k, double random_nums[], double* x, double* y, double* z);

    void c_cyl_pos_sample(void* src_pointer, int mesh_idx, int i, int j, int k, double random_nums[], double* x, double* y, double* z);
    
    void c_energy_selection(void* src_pointer, int mesh_idx, int selected_mesh_element, int* selected_cell, int energy_sample_type, double random_nums[], int ignore_cell_info_int, double* erg, double* wgt);
    
    void c_cell_adjust_weight(void* src_pointer, int mesh_idx, int selected_mesh_element, int selected_cell, double* wgt);
    
    void c_check_position(void* src_pointer, int mesh_idx, int selected_mesh_element, int picked_cell, double cell_importance, int material, int ignore_cell_info_int, int* resample_in_voxel_int, int* selected_cell, int verbose_int);
    
    void c_transform(double* x, double* y, double* z, double* u, double* v, double* w, double cosines[], double translations[], int verbose_int);
    
    int c_get_meshtype(void* src_pointer, int mesh_idx);
    
    int c_get_meshstruc(void* src_pointer, int mesh_idx);
    
    int c_get_meshtransform(void* src_pointer, int mesh_idx);
    
    void c_set_meshtransform(void* src_pointer, int mesh_idx, int trans_num);
    
    int c_get_cellid(void* src_pointer, int mesh_idx, int selected_mesh_element, int selected_cell);
    
    int c_get_cellnum(void* src_pointer,int mesh_idx, int selected_mesh_element);
    
#ifdef __cplusplus
}
#endif
    

#endif