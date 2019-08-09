! COMMON DECAY GAMMA SOURCE MODULES
! Zamir Ghani, Andrew Turner, Tim Eade. UKAEA (CCFE). 2017


module cR2S_modules

  implicit NONE

  integer, parameter, private    :: dknd   = selected_real_kind(15,307)   != 8-byte real kind
  integer, parameter, private    :: i4knd = selected_int_kind( 9)        != 4-byte integer kind
  integer, parameter, private    :: i8knd = selected_int_kind(18)        != 8-byte integer kind
  real(dknd), parameter, private :: cutoff = 1E-30  
  real, parameter, private  :: pi=3.1415926535898_dknd
  

  ! module to read the gamma decay source defined in the common format
  ! not all variables defined below are read form the file. Some of them
  ! are defined to simplify sampling process.


  ! ========== derived types ==============

  ! attributes of energy spectrum for cell
  type,private :: erg_attributes
    integer             ::  id           ! bin identification
    real(dknd)          ::  str          ! intensity of bin
    real(dknd)          ::  err          ! statistical error of bin
    real(dknd)          ::  str_CDF      ! intensity of bin CDF

  endtype erg_attributes

  ! attributes of cell inside voxel
  type,private :: cell_attributes
    integer             ::  id           ! cell identification inside the voxel
    real(dknd)          ::  frac         ! volume fraction of the cell inside voxel 
    real(dknd)          ::  str          ! intensity emitted by the fraction of the cell inside the voxel
    integer             ::  ebins_active  !since we store zero bins, we need to know the number of non-zero bins
    type(erg_attributes), allocatable :: erg(:) 

  endtype cell_attributes


  ! type defining voxels features
  type,private :: voxel_attributes
    integer             ::  id           ! voxel id 
    integer             ::  ncel         ! number of listed cells inside voxel
    real(dknd)          ::  vol          ! volume of the voxel 
    real(dknd)          ::  str          ! total strength of the voxel
    real(dknd)          ::  str_CDF      ! cumulative strength of the voxel (CDF)
    real(dknd)          ::  pos(3)       ! position coordinates of voxel center
    real(dknd)          ::  size(3)      ! voxel size in each dimension
    real(dknd)          ::  tot_cell_frac = 0.0  ! total of the volume fraction of the listed cells
    type(cell_attributes), allocatable :: cel(:) 
  endtype voxel_attributes

  ! type defining the mesh structure of the source
  type, private :: mesh_structure
    real(dknd), allocatable  :: xrbin(:)    ! x/r bins 
    real(dknd), allocatable  :: ytbin(:)    ! y/t bins
    real(dknd), allocatable  :: zbin(:)     ! z bins
    real(dknd), allocatable  :: ebin(:)     ! energy bins
    integer                  :: nxr         ! number of mesh points in x/r
    integer                  :: nyt         ! number of mesh points in y/theta
    integer                  :: nz          ! number of mesh points in z
    integer                  :: ne          ! number of energy points
    integer                  :: struc_ID    ! Identify mesh structure (1=structured, 2=unstructured)
    integer                  :: geom_ID     ! Identify the geometry of structured mesh(1=rec, 2=cyl)
    integer                  :: nvoxel      ! number of activated voxels
    integer                  :: spec_ID     ! Identify the gamma spectra structure (1=group, 2=line)
    integer                  :: trans_num=0 ! Transform number associated with the mesh (0=no transform,>0=TR ID in MCNP model) 

    real(dknd)               :: vol         ! activated volume of the mesh
    real(dknd)               :: str         ! total strength of the mesh (gamma/s)
    real(dknd)               :: str_CDF     ! cumulative strength of the voxel (CDF)
    real(dknd)               :: ori(3)      ! Origin of a cylindrical mesh 
    real(dknd)               :: axs(3)      ! Axis of a cylindrical mesh
    real(dknd)               :: vec(3)      ! Vector to define the starting location of theta on a cylindrical mesh

    character (LEN=128)      :: unstruc_mesh ! name of the unstructured mesh file
  endtype mesh_structure

 
  ! type encapsulating all source features
  type,public :: source_attributes
    type(mesh_structure)               :: mesh 
    type(voxel_attributes),allocatable :: voxel(:) 
  endtype



  


  ! Subroutines and functions contained in the module
  contains


  !--------------------------------------------------------------------------
  ! read data from common decay file module
  !--------------------------------------------------------------------------
  subroutine read_common_file(common_filename, src, nmeshes, total_src, verbose, error)

    implicit none

    !************************************************************************************************
    ! Subroutine dummy variables
    type(source_attributes), intent(out) , allocatable :: src(:)            ! Source meshes
    integer (i4knd), intent(out)                       :: nmeshes           ! Number of meshes contained in the commonfile
    real(dknd), intent(out)                            :: total_src         ! total source rate, all meshes in file
    character (len=*), intent(in)                      :: common_filename   ! Common file name
    logical, intent(in)                                :: verbose           ! Verbose output
    character (len=*), intent(inout)                   :: error             ! Error handler (returns errors string if occurs)
    !************************************************************************************************
    
    !*************************************************************************************************
    ! Local variables
    real(dknd) :: str                                                ! Temp source strength
    real(dknd) :: vol                                                ! Temp volume
    real(dknd) :: frac                                               ! Temp cell volume fraction
    real(dknd) :: dummy_real                                         ! to do read    
    real(dknd), allocatable  ::  spec_temp(:,:)                      ! spectra, temporary variable
    real(dknd) :: decay_time                                         ! Cooling time after shutdown (s)
    
    integer(i4knd) :: io                                             ! Input/output error
    integer(i4knd) :: ivxl                                           ! Voxel number (looping variable)
    integer(i4knd) :: id                                             ! Temp voxel/cell ID
    integer(i4knd) :: ncel                                           ! Temp number of cells under voxel
    integer(i4knd) :: icl                                            ! Cell looping variable
    integer(i4knd) :: elements                                       ! Number of unstructured mesh elements
    integer(i4knd) :: ebins_active                                   ! Temp var for number of active ebins
    integer(i4knd) :: dummy_int                                      ! Dummy integer
    integer(i4knd) :: num_unstr                                      ! Number of unstructured meshes (only 1 supported)
    integer(i4knd) :: jt                                             ! Temp source transform variables
    integer(i4knd) :: iunit                                          ! Common file unit number
    integer(i4knd) :: loop, loop2                                    ! Looping variable
    integer(i4knd) :: mesh_idx                                       ! Index of the current mesh
    
    logical :: store                                                 ! flag whether this voxel data is to be read and stored, or just read and ignored
    logical :: got_total_src                                         ! Logical to record whether total_src present
    logical :: got_spec_ID                                           ! Logical to record whether spec_ID present
    logical :: got_ne                                                ! Logical to record whether ne present
    logical :: got_mesh_type                                         ! Logical to record whether mesh_type present
    logical :: got_mesh_boundaries                                   ! Logical to record whether mesh_boundaries present
    
    character(len=1000) :: dummy_string                              ! Character string to read line
    character(len=100)  :: input_label                               ! Input label/key
    character(len=100)  :: mesh_type                                 ! Temp mesh type
    character(len=100)  :: erg_type                                  ! Temp energy type
    character(len=80)   :: title                                     ! Common file title     
    !*************************************************************************************************

    integer :: zgtr, ntr  !for transform card processing 

    !Try to open commonfile
    open(iunit,file=common_filename, status='old', action='read', iostat=io)
    
    !If not present then exit
    if(io /=0 )then
      write(*,*)'Error'
      error = commonfile_error('No source file present in directory')
      write(*,*)trim(error)
      return
    end if


    ! Set number of unstrucutred meshes to 0
    num_unstr = 0
    
    ! ---------   read file header information ------------

    if (verbose) write(*,*) 'Reading commonfile data'
    if (verbose) write(*,*) 'Reading number of meshes'

    !read first line: num_meshes
    read(iunit,*,iostat=io) input_label, nmeshes
    if (io .ne. 0 .OR. input_label /= 'num_meshes') then
       error = commonfile_error('Error on line 1, num_meshes') 
       return 
    endif

    allocate(src(nmeshes))  !for a multi mesh source


    multimesh: do mesh_idx=1,nmeshes,1
    
       ! Set the logical variables to false for each mesh
       got_total_src=.FALSE.
       got_spec_ID = .FALSE.
       got_ne=.FALSE. 
       got_mesh_type=.FALSE.
       got_mesh_boundaries=.FALSE.    
    
       ! ---------   read per mesh header information ------------
       if (verbose) write(*,*) 'Reading mesh ', mesh_idx

       !must start with mesh_id, title, cooling time.

       read(iunit,*,iostat=io) input_label, dummy_int
       if (io .ne. 0 .OR. input_label /= 'mesh_id' .OR. dummy_int .ne. mesh_idx) then
          write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', mesh_id'
          error = commonfile_error(dummy_string) 
          return 
       endif

       !this is title, cooling time No optional parameters.
       read(iunit,'(A)') title

       !read  line: cooling_time
       read(iunit,*,iostat=io) input_label, decay_time
       if (io .ne. 0 .OR. input_label /= 'cooling_time') then  !EoF no longer valid end condition
          dummy_string = 'Error - expecting cooling_time'
          error = commonfile_error(dummy_string) 
          return 
       endif

       ! the other parameters can be out of order, since some are optional

       readlines: do 

          read(iunit,'(A)',iostat=io) dummy_string
          if (verbose) write(*,*) dummy_string

          if (io .ne. 0) then  !EoF no longer valid end condition
             dummy_string = 'Commonfile mesh header read error'
             error = commonfile_error(dummy_string) 
             return 
          endif
          dummy_string = trim(adjustl(dummy_string))
          if (len(dummy_string) .eq. 0) cycle readlines  ! handles blank lines

          read(dummy_string,*) input_label

          if (input_label == 'total_source') then
             read(dummy_string,*,iostat=io) input_label, src(mesh_idx)%mesh%str
             got_total_src = .TRUE.
             cycle readlines
          endif

          if (input_label == 'energy_type') then
             read(dummy_string,*,iostat=io) input_label, erg_type
             if (io .ne. 0) then
                write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', energy_format incorrect'
                error = commonfile_error(dummy_string) 
                return 
             endif
  
             got_spec_ID = .TRUE.

             erg_type = trim(adjustl(erg_type))
             if (erg_type == 'bins') then
                src(mesh_idx)%mesh%spec_ID = 1
             elseif (erg_type == 'lines') then
                !line sampling not yet supported
                write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,': Line energies not currently supported'
                error = commonfile_error(dummy_string) 
                return 
             else
                ! not a supported energy type keyword
                write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', erg_type incorrect'
                error = commonfile_error(dummy_string) 
                return 
             endif

             cycle readlines
          endif

          if (input_label == 'energy_boundaries') then
             if (.NOT. got_spec_ID) then  !this needs to have been read already
                write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,' - boundaries but energy_format not defined'
                error = commonfile_error(dummy_string) 
             endif

             read(dummy_string,*,iostat=io) input_label, src(mesh_idx)%mesh%ne
             if (io .ne. 0) then
                write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', number of energy groups incorrect'
                error = commonfile_error(dummy_string) 
                return 
             endif
             allocate(src(mesh_idx)%mesh%ebin(src(mesh_idx)%mesh%ne))

             !read boundaries
             read(iunit,*,iostat=io) src(mesh_idx)%mesh%ebin
             if (io .ne. 0) then
                write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', energy boundaries incorrect'
                error = commonfile_error(dummy_string) 
                return 
             endif
             ! Energy ne and bins read OK
             got_ne = .TRUE. 
             cycle readlines
          endif

          if (input_label == 'mesh_type') then
             read(dummy_string,*,iostat=io) input_label, mesh_type
             mesh_type = trim(adjustl(mesh_type))
             if (mesh_type == 'rec' .OR. mesh_type == 'rectangular') then
                src(mesh_idx)%mesh%struc_ID = 1
                src(mesh_idx)%mesh%geom_ID = 1
             else if (mesh_type == 'cyl' .OR. mesh_type == 'cylindrical') then
                src(mesh_idx)%mesh%struc_ID = 1
                src(mesh_idx)%mesh%geom_ID = 2
             else if (mesh_type == 'unstr' .OR. mesh_type == 'unstructured') then
                src(mesh_idx)%mesh%struc_ID = 2
                num_unstr = num_unstr+1
                if (num_unstr .gt. 1) then
                   ! only 1 unstructured mesh supported
                   write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', more than one unstructured mesh'
                   error = commonfile_error(dummy_string) 
                   return 
                endif
             else
                ! not a supported mesh type keyword
                write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', mesh_type incorrect'
                error = commonfile_error(dummy_string) 
                return 
             endif
             got_mesh_type = .TRUE.
             cycle readlines
          endif

          if (input_label == 'mesh_boundaries') then
             if (.NOT. got_mesh_type) then  ! meshtype must exist before mesh boundaries are read 
                write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', mesh_type not read'
                error = commonfile_error(dummy_string) 
                return 
             endif

             if (src(mesh_idx)%mesh%struc_ID == 1 ) then

               !read structured mesh dimension (XYZ/RTZ)
                read(dummy_string,*,iostat=io) input_label, src(mesh_idx)%mesh%nxr,src(mesh_idx)%mesh%nyt,src(mesh_idx)%mesh%nz
                if (io .ne. 0) then
                   write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', mesh boundaries incorrect'
                   error = commonfile_error(dummy_string) 
                   return 
                endif
 
                if(src(mesh_idx)%mesh%geom_ID==2) then !read 
                   read(iunit,*,iostat=io)src(mesh_idx)%mesh%ori(1:3),src(mesh_idx)%mesh%axs(1:3),src(mesh_idx)%mesh%vec(1:3)

                   if (io .ne. 0) then
                      write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', cylinder ori/axs/vec incorrect'
                      error = commonfile_error(dummy_string) 
                      return 
                   endif

                end if

                ! read first coordinate spacing
                allocate(src(mesh_idx)%mesh%xrbin(src(mesh_idx)%mesh%nxr))
                read(iunit,*,iostat=io) (src(mesh_idx)%mesh%xrbin(loop),loop=1,src(mesh_idx)%mesh%nxr)
                if (io .ne. 0) then
                   write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', xr mesh grid incorrect'
                   error = commonfile_error(dummy_string) 
                   return 
                endif

                ! read second coordinate spacing
                allocate(src(mesh_idx)%mesh%ytbin(src(mesh_idx)%mesh%nyt))
                read(iunit,*,iostat=io) (src(mesh_idx)%mesh%ytbin(loop),loop=1,src(mesh_idx)%mesh%nyt)
                if (io .ne. 0) then
                   write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', yt mesh grid incorrect'
                   error = commonfile_error(dummy_string) 
                   return 
                endif

                ! read third coordinate spacing
                allocate(src(mesh_idx)%mesh%zbin(src(mesh_idx)%mesh%nz))
                read(iunit,*,iostat=io) (src(mesh_idx)%mesh%zbin(loop),loop=1,src(mesh_idx)%mesh%nz)
                if (io .ne. 0) then
                   write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', z mesh grid incorrect'
                   error = commonfile_error(dummy_string) 
                   return 
                endif

                allocate(src(mesh_idx)%voxel(src(mesh_idx)%mesh%nxr*src(mesh_idx)%mesh%nyt*src(mesh_idx)%mesh%nz))

             else !else if unstructured mesh read in the number of elements and allocate arrays

                read(dummy_string,*,iostat=io) input_label, src(mesh_idx)%mesh%unstruc_mesh, elements !read path, and num elements
                if (io .ne. 0) then
                   write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,', unstructured mesh parameters incorrect'
                   error = commonfile_error(dummy_string) 
                   return 
                endif

                ! As src(mesh_idx)%mesh%nxr is actually the number of boundaries need to add one to the number of elements
                src(mesh_idx)%mesh%nxr=elements+1

                !set the y and z to 2 to allow the structured mesh routines to work on the unstructured mesh
                src(mesh_idx)%mesh%nyt=2
                src(mesh_idx)%mesh%nz=2

                !allocate mesh boundary bins and voxels
                allocate(src(mesh_idx)%mesh%xrbin(src(mesh_idx)%mesh%nxr))
                allocate(src(mesh_idx)%mesh%ytbin(src(mesh_idx)%mesh%nyt))
                allocate(src(mesh_idx)%mesh%zbin(src(mesh_idx)%mesh%nz))
                allocate(src(mesh_idx)%voxel(src(mesh_idx)%mesh%nxr-1))
             endif
             got_mesh_boundaries = .TRUE.
             cycle readlines
          endif

          if (input_label == 'source_data') then
             !end of the mesh header
             exit readlines
          endif

          ! --- optional mesh header variables ---

          if (input_label == 'mesh_transform') then
             read(dummy_string,*,iostat=io) input_label, src(mesh_idx)%mesh%trans_num
             cycle readlines
          endif

          ! --- Unprocessed header item --- list this last...
          ! only gets to here if a mesh entry is not recognised
          write(dummy_string,'(A,I2,A,A)') 'mesh ', mesh_idx,' - unrecognised: ',trim(adjustl(input_label))
          error = commonfile_error(dummy_string) 
          return 

       end do readlines

       !check program has successfully read all compulsory mesh variables
       if (.NOT. got_total_src) then
          write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,': total_source missing'
          error = commonfile_error(dummy_string) 
          return 
       endif

       if (.NOT. got_mesh_boundaries) then
          write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,': mesh_boundaries missing'
          error = commonfile_error(dummy_string) 
          return 
       endif

      if (.NOT. got_mesh_type) then
          write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,': mesh_type missing'
          error = commonfile_error(dummy_string) 
          return 
       endif

       if (.NOT. got_spec_ID) then
          write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,': energy_format missing'
          error = commonfile_error(dummy_string) 
          return 
       endif

       if (.NOT. got_ne) then
          write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,': energy_boundaries missing'
          error = commonfile_error(dummy_string) 
          return 
       endif


       ! ----------  read the mesh data for this mesh ---------------------------------


       if (verbose) write(*,*) 'mesh ', mesh_idx, ' header info read successfully...reading voxels'

       !-------------------------
       ! read voxel information !
       !-------------------------

       !  loop over voxels
       !-------------------
       ! the CDGS file should contain only non-zero voxels, however, filter these out at this stage
       ivxl=0
       voxel_loop: do

         read(iunit,'(A)',iostat=io) dummy_string  !for checking 'end_source_data'

         
        ! check read / EOF
         if (io .lt. 0) then       ! end of file
            write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,': End of file before end_source_data'
            error = commonfile_error(dummy_string) 
            return 
         endif
         if (io .gt. 0 ) then
            write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,': source data read error'
            error = commonfile_error(dummy_string) 
            return 
         endif
         if(verbose) write(*,*)dummy_string
         
         if (dummy_string .eq. 'end_source_data') then
           if(verbose) write(*,*)'Exiting voxel loop'
           exit voxel_loop      ! end of mesh
         end if

         ! this must be the start of a new voxel
         read(dummy_string,*,iostat=io) id,str,vol,ncel
         if (io .ne. 0) then      
            write(dummy_string,'(A,I2,A)') 'mesh ', mesh_idx,': read error id/str/vol/ncel'
            error = commonfile_error(dummy_string) 
            return 
         endif



         ! Read OK and not EOF...
         store=.TRUE.
         if (str/src(mesh_idx)%mesh%str .lt. cutoff) store=.FALSE.   !if intensity/total lower than cut-off then read over and ignore

         if (store) then
            ivxl=ivxl+1 
            src(mesh_idx)%voxel(ivxl)%id=id
            src(mesh_idx)%voxel(ivxl)%str=str
            src(mesh_idx)%voxel(ivxl)%vol=vol
            src(mesh_idx)%voxel(ivxl)%ncel=ncel
         endif

         ! read cell under voxel data and store if needed
         !  loop over cell under voxel
         !----------------------------
         ! it is in the definition of the source that only the 'active' / non-zero cells are given.
         ! however, in cases where cells with zero intensity are also given, this is supported as well: stored and rejected during sampling step.
         ! this uses random spatial sampling - as such the total volume fraction needs to exclude any later rejections

         if (store) then
            allocate(src(mesh_idx)%voxel(ivxl)%cel(ncel))

            cell_in_voxel_loop: do icl=1,ncel
               read(iunit,*,iostat=io) id,frac,str
               if (io .ne. 0) then      
                  write(dummy_string,'(A,I2,A,I9,A)') 'mesh ', mesh_idx,' voxel: ',src(mesh_idx)%voxel(ivxl)%id,&
                                                      ': read error id/frac/str'
                  error = commonfile_error(dummy_string) 
                  return 
               endif

               src(mesh_idx)%voxel(ivxl)%cel(icl)%id=id
               src(mesh_idx)%voxel(ivxl)%cel(icl)%frac=frac
               src(mesh_idx)%voxel(ivxl)%cel(icl)%str=str

               if (str/src(mesh_idx)%mesh%str .ge. cutoff)  src(mesh_idx)%voxel(ivxl)%tot_cell_frac = &
                      src(mesh_idx)%voxel(ivxl)%tot_cell_frac + frac  !for counting total cell fraction, which is intiialised to zero

               allocate(spec_temp(src(mesh_idx)%mesh%ne - 1,2))     
     

               ! read photon spectra and associate error
               do loop = 1,2
                  read(iunit,*,iostat=io) (spec_temp(loop2,loop),loop2=1,src(mesh_idx)%mesh%ne - 1) 
               enddo 

               if (io .ne. 0) then      
                  write(dummy_string,'(A,I2,A,I9,A)') 'mesh ', mesh_idx,' voxel: ',src(mesh_idx)%voxel(ivxl)%id,&
                                                      ': read error in spectra'
                  error = commonfile_error(dummy_string) 
                  return 
               endif

               ! create 'active' ebins list for this cell, to avoid rejection step for uniform sampling

               ebins_active = 0  ! count them
               do loop=1,src(mesh_idx)%mesh%ne - 1 
                  if (spec_temp(loop,1)/src(mesh_idx)%mesh%str .ge. cutoff) ebins_active=ebins_active+1
               end do
               src(mesh_idx)%voxel(ivxl)%cel(icl)%ebins_active = ebins_active
               allocate(src(mesh_idx)%voxel(ivxl)%cel(icl)%erg(ebins_active))     

               ebins_active=0 !(reusing ebins_active as a counter here)
               do loop=1,src(mesh_idx)%mesh%ne - 1              
                  if (spec_temp(loop,1)/src(mesh_idx)%mesh%str .ge. cutoff) then
                     ebins_active=ebins_active+1
                     src(mesh_idx)%voxel(ivxl)%cel(icl)%erg(ebins_active)%id = loop
                     src(mesh_idx)%voxel(ivxl)%cel(icl)%erg(ebins_active)%str = spec_temp(loop,1)
                     src(mesh_idx)%voxel(ivxl)%cel(icl)%erg(ebins_active)%err = spec_temp(loop,2)
                  endif
               end do
               deallocate(spec_temp) !for re-use

               !renormalise to CDF for this cells emission spectra
               call construct_CDF(src(mesh_idx)%voxel(ivxl)%cel(icl)%erg(:)%str,src(mesh_idx)%voxel(ivxl)%cel(icl)%erg(:)%str_CDF)
            enddo cell_in_voxel_loop

         else !read over
            if (verbose) write(*,*) 'mesh ', mesh_idx, ' skipping voxel ',id

            cell_in_voxel_loop_read_only: do icl=1,ncel
               read(iunit,*,iostat=io) id,frac,str
               if (io .ne. 0) then      
                  write(dummy_string,'(A,I2,A,I9,A)') 'mesh ', mesh_idx,' voxel: ',src(mesh_idx)%voxel(ivxl)%id,&
                                                      ': read error id/frac/str'
                  error = commonfile_error(dummy_string) 
                  return 
               endif
     
               ! read photon spectra and associate error
               do loop = 1,2
                  read(iunit,*,iostat=io) (dummy_real,loop2 = 1,src(mesh_idx)%mesh%ne - 1)
               enddo 
               if (io .ne. 0) then      
                  write(dummy_string,'(A,I2,A,I9,A)') 'mesh ', mesh_idx,' voxel: ',src(mesh_idx)%voxel(ivxl)%id,&
                                                      ': read error in spectra'
                  error = commonfile_error(dummy_string) 
                  return 
               endif
            enddo cell_in_voxel_loop_read_only
             
         endif
     
       enddo voxel_loop

       src(mesh_idx)%mesh%nvoxel=ivxl   !store number of active voxels

       total_src = total_src + src(mesh_idx)%mesh%str  !accumulate a total source

       ! --- populate mesh voxel strength cumulative array ---
       call construct_CDF(src(mesh_idx)%voxel(:)%str,src(mesh_idx)%voxel(:)%str_CDF)     ! convert voxel strength data to a CDF

    end do multimesh

    ! --- populate mesh total strength cumulative array ---
    call construct_CDF(src(:)%mesh%str,src(:)%mesh%str_CDF)     ! convert voxel strength data to a CDF
    return
  end subroutine read_common_file
  !-----------------------------------------------------------------------------------

  
  !-----------------------------------------------------------------------------------  
  subroutine mesh_selection(src, sample_type, random_num, total_src, mesh_idx, wgt, verbose)
  
     implicit none

     
    !**********************************************************************************************************
    ! Subroutine dummy variables   
    type(source_attributes), intent(in) :: src(:)        ! Source meshes

    integer(i4knd), intent(in)          :: sample_type   ! Type of sampling mode (see note below)
 
    logical, intent(in)                 :: verbose       ! Write debug info to screen

    real(dknd), intent(in)              :: random_num    ! Random number between 0-1 used to pick mesh
    real(dknd), intent(in)              :: total_src     ! Total source strength
    
    integer(i4knd), intent(out)         :: mesh_idx      ! Picked mesh IDX value
    real(dknd), intent(inout)           :: wgt           ! Weight of the starting particle for the picked mesh 
    !**********************************************************************************************************
     
     
     !two modes, uniform/weight correction, and analogue/proportional to source intensity
     !chosen by the code to be the same as the setting for the voxel in mesh selection.
     !sets mesh_idx (sample_type = 1 analogue sampling, sample_type /= 1 uniform sampling)

     ! Subroutine returns the picked mesh idx and the weight of the starting particle for the mesh picked
     
     if (sample_type .eq. 1) then   ! analogue sampling proportional to intensity
     
        mesh_idx = sample_CDF_bin(random_num,src(:)%mesh%str_cdf)
        if(verbose) write(*,*) 'Mesh ', mesh_idx, ' selected'
        
     else      ! uniform sampling with weight correction
        mesh_idx = ceiling(random_num*size(src))   !Select a mesh between 1 and the highest numbered mesh

        if(verbose) write(*,*) 'Mesh ', mesh_idx, ' selected'

        !weighting == number of meshes * MESH INTENSITY / TOTAL INTENSITY           
        wgt=wgt*size(src)*src(mesh_idx)%mesh%str/total_src
        if(verbose) write(*,*) 'Mesh weight ', wgt

     endif

     return

  end subroutine mesh_selection
  !-----------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------
  ! Create error string
  function commonfile_error(msg) result(error_msg)
     implicit none
     character(len=*) :: msg
     character(1000) :: error_msg

     error_msg = 'cR2S read error: ' // trim(adjustl(msg))

  end function commonfile_error
  !-----------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------------
  ! subroutine to select a voxel and set the particles weight  
  !-------------------------------------------------------------------------------------------
  subroutine voxel_selection(src, sample_type, random_num, selected_voxel, selected_mesh_element, i, j, k, wgt, verbose)

    implicit none

    !**********************************************************************************************************
    ! Subroutine dummy variables      
    type(source_attributes), intent(in) :: src                   ! Chosen source mesh
    integer(i4knd), intent(in)          :: sample_type           ! Type of sampling mode (sample_type = 1 analogue sampling, sample_type /= 1 uniform sampling)
    real(dknd), intent(in)              :: random_num            ! Random number between 0-1 used to pick voxel
    logical, intent(in)                 :: verbose               ! Write debug info to screen    
    
    integer(i4knd), intent(out)         :: selected_voxel        ! Selected voxel ID
    integer(i4knd), intent(out)         :: selected_mesh_element ! Reference to the selected element in the mesh voxel array 
    integer(i4knd), intent(out)         :: i, j, k               ! x, y and z index of the voxel
    real(dknd), intent(inout)           :: wgt                   ! Weight of the starting particle 
    !**********************************************************************************************************

    !**********************************************************************************************************
    ! Local variables
    integer(i4knd)    :: size_x, size_y, size_z                  ! size of mesh in x,y,z dimensions
    !**********************************************************************************************************
    
    ! Select a voxel from the list (which are the NON-ZERO INTENSITY VOXELS)
    ! the alternative would be for the uniform sampling to add rejection step here.

    if (sample_type .eq. 1) then  !analogue spatial sampling
    
       selected_mesh_element = sample_CDF_bin(random_num,src%voxel(:)%str_cdf)
       
    else  !uniform spatial sampling
    
       selected_mesh_element = ceiling(random_num*src%mesh%nvoxel)   !Select an element between 1 and the highest active element
       ! weighting == no. of voxels*VOXEL INTENSITY / TOTAL INTENSITY           
       wgt=wgt*(src%voxel(selected_mesh_element)%str/src%mesh%str)*src%mesh%nvoxel

    endif

    ! Get the voxel ID from the selected element
    selected_voxel=src%voxel(selected_mesh_element)%id      

    ! translate nth voxel to i,j,k array element
    size_x= src%mesh%nxr-1
    size_y= src%mesh%nyt-1
    size_z= src%mesh%nz-1

    !Work out the i, j and k values for the voxel
    i=((selected_voxel-1)/(size_y*size_z))+1
    j=(((selected_voxel-1)-(i-1)*(size_y*size_z))/size_z)+1
    k= selected_voxel-((i-1)*(size_y*size_z))-((j-1)*size_z)
 
    return

  end subroutine voxel_selection
  !---------------------------------------------------------------------------------




  !-----------------------------------------------------------------------------------
  ! Pick a position in a rectangular voxel
  !---------------------------------------------------------------------------------
  subroutine rec_pos_sample(src, i, j, k, random_nums, x, y, z)
  
    !**********************************************************************************************************
    ! Subroutine dummy variables   
    type(source_attributes), intent(in) :: src              ! Source meshes
    integer(i4knd), intent(in)          :: i, j, k          ! Voxel index
    real(dknd), intent(in)              :: random_nums(1:3)  ! Three random numbers between 0-1 used to pick position
    
    real(dknd), intent(out)             :: x, y, z          ! X, Y, Z position in voxel
    !**********************************************************************************************************

    !**********************************************************************************************************
    ! Local variables    
    real(dknd) :: xmin,xmax,ymin,ymax,zmin,zmax  
    !**********************************************************************************************************
        
    !Set voxel bin boundary based on i,j,k values calculated in voxel select routine.
    !print *, "i =", i , "j=", j, "k=" , k
    xmin=src%mesh%xrbin(i)
    xmax=src%mesh%xrbin(i+1)
    ymin=src%mesh%ytbin(j)
    ymax=src%mesh%ytbin(j+1)
    zmin=src%mesh%zbin(k)
    zmax=src%mesh%zbin(k+1)
    !print *, "xmin=",xmin,"ymin=",ymin,"zmin=",zmin
    !print *, "xmax=",xmax,"ymax=",ymax,"zmax=",zmax

    !Choose point at random in the voxel
    x=xmin+random_nums(1)*(xmax-xmin)
    y=ymin+random_nums(2)*(ymax-ymin)
    z=zmin+random_nums(3)*(zmax-zmin)

  end subroutine rec_pos_sample
  !---------------------------------------------------------------------------------




  !---------------------------------------------------------------------------------
  ! Pick a position in a cylindrical voxel
  !---------------------------------------------------------------------------------
  subroutine cyl_pos_sample(src, i, j, k, random_nums, x, y, z)
  
    !**********************************************************************************************************
    ! Subroutine dummy variables   
    type(source_attributes), intent(in) :: src              ! Source meshes
    integer(i4knd), intent(in)          :: i, j, k          ! Voxel index
    real(dknd), intent(in)              :: random_nums(1:3) ! Three random numbers between 0-1 used to pick position
    
    real(dknd), intent(out)             :: x, y, z          ! X, Y, Z position in voxel
    !**********************************************************************************************************

    !**********************************************************************************************************
    ! Local variables    
    real(dknd) :: rmin,rmax,tmin,tmax,zmin,zmax
    real(dknd) :: XYZ(3), RTZ(3), TempVec(3), TempAxs(3), InitialXYZ(3), XYZ2(3)
    real(dknd) :: VecAngle
    !**********************************************************************************************************

    !set up the voxel bin boundaries based on picked i,j,k values from voxel select routine 
    rmin=src%mesh%xrbin(i)
    rmax=src%mesh%xrbin(i+1)
    tmin=src%mesh%ytbin(j)
    tmax=src%mesh%ytbin(j+1)
    zmin=src%mesh%zbin(k)
    zmax=src%mesh%zbin(k+1)


    ! Pick an R value (biased to pick uniformly in area)
    RTZ(1) = (random_nums(1)*(rmax**2-rmin**2)+rmin**2)**0.5
    ! Pick a z value uniformly
    RTZ(2) = random_nums(2) * (tmax - tmin) + tmin
    !Pick a theta value uniformly
    RTZ(3) = random_nums(3) * (zmax - zmin) + zmin

    !Transform point into cartesian coordinate for cylinder on z axis
    XYZ(1) = RTZ(1)*cos(RTZ(2)*2*pi)
    XYZ(2) = RTZ(1)*sin(RTZ(2)*2*pi)
    XYZ(3) = RTZ(3)
    !WRITE(*,*)rmin,rmax,tmin,tmax,zmin,zmax,RTZ,xyz
    
    !Save initial point
    ! InitialXYZ=XYZ

    !Calculate the angle between 1,0,0 and 'vec'
    TempVec(1) = 1
    TempVec(2) = 0
    TempVec(3) = 0

    !If not on axis then rotate    
    if((src%mesh%axs(1)/=0 .or. src%mesh%axs(2)/=0 .or. src%mesh%axs(3)/=1) .or.&
       (src%mesh%vec(1)/=1 .or. src%mesh%vec(2)/=1 .or. src%mesh%vec(3)/=1))then

      ! If cylinder is not aligned with z axis:
      if(src%mesh%axs(1)/=0 .or. src%mesh%axs(2)/=0)then
        ! Rotate tempvec to be in the same plane as vec 
        TempVec=Rotate(src%mesh%axs, TempVec)

        !If temp vector and cylinder vector are the same, angle between is 0 otherwise work it out.
        if(TempVec(1)==src%mesh%vec(1) .and. TempVec(2)==src%mesh%vec(2) .and. &
           TempVec(3)==src%mesh%vec(3)) then
          VecAngle=0
        else
          VecAngle=acos(dot_product(TempVec,src%mesh%vec)/&
          (sqrt(TempVec(1)**2+TempVec(2)**2+TempVec(3)**2)*sqrt(src%mesh%vec(1)**2+&
          src%mesh%vec(2)**2+src%mesh%vec(3)**2)))
        end if

      ! If cylinder is aligned with z axis:
      elseif(src%mesh%axs(1)==0 .and. src%mesh%axs(2)==0)then

        VecAngle=acos(dot_product(TempVec,src%mesh%vec)/&
        (sqrt(TempVec(1)**2+TempVec(2)**2+TempVec(3)**2)*sqrt(src%mesh%vec(1)**2+&
        src%mesh%vec(2)**2+src%mesh%vec(3)**2)))

      end if

    
      ! Rotate x and y values by VecAngle
      XYZ2(1) = cos(VecAngle)*XYZ(1)-sin(VecAngle)*XYZ(2)
      XYZ2(2) = sin(VecAngle)*XYZ(1)+cos(VecAngle)*XYZ(2)

      XYZ(1) = XYZ2(1)
      XYZ(2) = XYZ2(2)
      !IF(myid==1) WRITE(111,'(5f10.4)')InitialXYZ,XYZ,VecAngle


      !If cylinder is off z-axis rotate point to be in line with axis of the cylinder
      if(src%mesh%axs(1)/=0 .or. src%mesh%axs(2)/=0)then

        XYZ=Rotate(src%mesh%axs,XYZ)
        ! IF(myid==1) WRITE(105,'(A,9f10.4)')'Radius sampling',lowX,highX,RTZ(1),InitialXYZ,XYZ

      !If cylinder is in the -z direction flip z coordinate of point      
      elseif(src%mesh%axs(1)==0 .and. src%mesh%axs(2)==0 .and. src%mesh%axs(3)<0)then

        XYZ(3) = -XYZ(3)     

      end if
      
    end if

    !Shift the origin of the cylinder
    XYZ = XYZ + src%mesh%ori

    !Set the x,y,z mcnp variables
    x = XYZ(1)
    y = XYZ(2)
    z = XYZ(3)

    !WRITE(110,'(f10.5,i7,6f10.5)')pbl%r%wgt,selected_voxel,RTZ,pbl%r%x,pbl%r%y,pbl%r%z
  end subroutine cyl_pos_sample
  !---------------------------------------------------------------------------------





  !---------------------------------------------------------------------------------
  ! Function to rotate a point 
  !---------------------------------------------------------------------------------
  function rotate(CylAxis,XYZ) result(XYZout)
  
    !**********************************************************************************************************
    ! Subroutine dummy variables   
    real(dknd), intent(in) :: CylAxis(:),XYZ(:)
    real(dknd) :: XYZout(3)
    !**********************************************************************************************************
      
    !**********************************************************************************************************
    ! Local variables        
    real(dknd) :: Transform(1:3,1:3), Transform2(1:3,1:3), alpha,beta,gamma
    
    Transform(:,:) = 0
    Transform2(:,:) = 0
    
    ! Rotation in XZ plane
    Transform2(1,1)=CylAxis(3)/((CylAxis(1)**2+CylAxis(2)**2+CylAxis(3)**2)**0.5)
    Transform2(1,2)=0
    Transform2(1,3)=((CylAxis(1)**2+CylAxis(2)**2)**0.5)/((CylAxis(1)**2+CylAxis(2)**2+CylAxis(3)**2)**0.5)
    Transform2(2,1)=0
    Transform2(2,2)=1
    Transform2(2,3)=0
    Transform2(3,1)=-((CylAxis(1)**2+CylAxis(2)**2)**0.5)/((CylAxis(1)**2+CylAxis(2)**2+CylAxis(3)**2)**0.5)
    Transform2(3,2)=0
    Transform2(3,3)=CylAxis(3)/((CylAxis(1)**2+CylAxis(2)**2+CylAxis(3)**2)**0.5)

    ! Rotation away from XZ plane
    Transform(1,1)=CylAxis(1)/((CylAxis(1)**2+CylAxis(2)**2)**0.5)
    Transform(1,2)=-CylAxis(2)/((CylAxis(1)**2+CylAxis(2)**2)**0.5)
    Transform(1,3)=0
    Transform(2,1)=CylAxis(2)/((CylAxis(1)**2+CylAxis(2)**2)**0.5)
    Transform(2,2)=CylAxis(1)/((CylAxis(1)**2+CylAxis(2)**2)**0.5)
    Transform(2,3)=0
    Transform(3,1)=0
    Transform(3,2)=0
    Transform(3,3)=1

    !Perform rotations
    XYZout = matmul(Transform2,XYZ)
    XYZout = matmul(Transform,XYZout)

  end function rotate
  !---------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------
  ! Sample the particle energy  
  !---------------------------------------------------------------------------------
  subroutine energy_selection(voxel, selected_cell, ebins, energy_sample_type, random_num, ignore_cell_info, erg, wgt)

    implicit none

    !**********************************************************************************************************
    ! Subroutine dummy variables 
    type(voxel_attributes), intent(in)  :: voxel              ! Voxel to be sampled from
    real(dknd), intent(in)              :: ebins(:)           ! Gamma energy bin boundaries
    real(dknd), intent(in)              :: random_num(1:2)    ! Two random numbers between 0-1 used to pick position
    integer(i4knd), intent(in)          :: energy_sample_type ! Type of sampling mode (energy_sample_type = 0 analogue sampling, energy_sample_type /= 1 uniform sampling)
    integer(i4knd), intent(inout)       :: selected_cell      ! The selected cell in the voxel (passed in if cell under voxel used otherwise set to zero)
    
    logical, intent(in)                 :: ignore_cell_info   ! Logical to ignore the cell info and just use first cell under voxel
    
    
    real(dknd), intent(out)             :: erg                ! Energy of the starting particle
    real(dknd), intent(inout)           :: wgt                ! Weight of the starting particle 
    !**********************************************************************************************************
    
    
    !**********************************************************************************************************
    ! Local variables  
    real(dknd):: ergt,ergrn
    integer   :: selected_ebin
    !**********************************************************************************************************

    
    ! new format src(mesh_idx)%voxel(ivxl)%cel(icl)%erg(:)%str
    ! >>   the list on entering here is the non-zero groups only!  <<

    ! ergt for sampling energy group
    ergt=random_num(1)
    ! ergrn for sampling energy randomly within the group
    ergrn=random_num(2)

    if(voxel%cel(1)%id .eq. 0 .OR. ignore_cell_info) selected_cell=1
    !if we have valid cell data (1st cell > 0), selected_cell previously determined.  Else set it to 1.
    
    ! assign the energy of the particle
    ! find the energy bin and then sample from it.

    if (energy_sample_type .eq. 1) then   !uniform sampling of energy bins with weight.
       selected_ebin = ceiling(ergt*voxel%cel(selected_cell)%ebins_active)   !Select with equal probability from the list of active bins
      ! weighting == no. of --choosable-- ebins* bin intensity / total intensity of all bins (i.e. of cell)  	 
       wgt=wgt*voxel%cel(selected_cell)%ebins_active*&
                 voxel%cel(selected_cell)%erg(selected_ebin)%str/&
                 voxel%cel(selected_cell)%str 
    else   ! analogue sampling proportional to intensity (default)
       selected_ebin = sample_CDF_bin(ergt, voxel%cel(selected_cell)%erg(:)%str_CDF)
    endif

    ! randomly pick within the group bounds
    ! energy bins are displaced by 1 (1 = lower bound, 2 = upper bound bin 1 and so forth)

    ! convert to energy bin ID
    selected_ebin = voxel%cel(selected_cell)%erg(selected_ebin)%id
    erg=ebins(selected_ebin) + ergrn*(ebins(selected_ebin+1) - ebins(selected_ebin))

  end subroutine energy_selection
  !--------------------------------------------------------------------------


  !---------------------------------------------------------------------------------
  ! Subroutine to adjust the weight if cell under voxel rejection is used
  !---------------------------------------------------------------------------------
  subroutine cell_adjust_weight(voxel, selected_cell, wgt)

    !**********************************************************************************************************
    ! Subroutine dummy variables 
    type(voxel_attributes), intent(in)  :: voxel            ! Voxel to be sampled from
    integer(i4knd), intent(in   )       :: selected_cell    ! The selected cell in the voxel (passed in if cell under voxel used otherwise set to zero)    
    real(dknd), intent(inout)           :: wgt              ! Weight of the starting particle 
    !**********************************************************************************************************
      
      ! adjust weight for fractional intensity of cell sampled within voxel
      wgt=wgt*(voxel%cel(selected_cell)%str/&
               (voxel%str))
      ! adjust for volume fraction vs permitted sampling volume fraction (non-rejected)
      wgt=(wgt*voxel%tot_cell_frac)/&
                voxel%cel(selected_cell)%frac
    

  end subroutine cell_adjust_weight
  !---------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------
  ! Subroutine to check the selected spatial point 
  !---------------------------------------------------------------------------------
  subroutine check_position(src, mesh_idx, selected_mesh_element, picked_cell, cell_importance, material, &
                            ignore_cell_info, resample_in_voxel, selected_cell, verbose)
    
    !**********************************************************************************************************
    ! Subroutine dummy variables 
    type(source_attributes), intent(in) :: src(:)
    integer(i4knd), intent(in) :: mesh_idx
    integer(i4knd), intent(in) :: selected_mesh_element
    integer(i4knd), intent(in) :: picked_cell
    integer(i4knd), intent(in) :: material
    integer(i4knd), intent(out) :: selected_cell
    
    real(dknd), intent(in) :: cell_importance

    logical, intent(in) :: ignore_cell_info
    logical, intent(in) :: verbose
    logical, intent(out) :: resample_in_voxel
    !**********************************************************************************************************

    
    !**********************************************************************************************************
    ! Local variables   
    logical        :: cell_selection
    integer(i4knd) :: cell_search
    !**********************************************************************************************************
    
    cell_selection=.FALSE.
    resample_in_voxel=.FALSE.
    selected_cell = 0

    ! Check if the cell is on the list of cells in the source file (which should not be zero intensity, but potentially could be...)

    !if the first cell is non-zero cell number, and we have not turned off cell data
    if((src(mesh_idx)%voxel(selected_mesh_element)%cel(1)%id .gt. 0) .AND. (ignore_cell_info .eqv. .False.) ) then
      cell_selection=.FALSE.
    ! convert mcnp cell number to input deck cell number
      do cell_search=1, src(mesh_idx)%voxel(selected_mesh_element)%ncel
        if (picked_cell==src(mesh_idx)%voxel(selected_mesh_element)%cel(cell_search)%id) then
          selected_cell=cell_search
          cell_selection=.TRUE.
        endif

        !Exit loop if this is the correct cell
        if (cell_selection.EQV..TRUE.) exit
      enddo

      ! check intensity of emission of cell, to apply rejection if zero intensity cells listed by accident in the source file
      if (src(mesh_idx)%voxel(selected_mesh_element)%cel(selected_cell)%str/src(mesh_idx)%mesh%str .lt. cutoff) then
         cell_selection = .FALSE.
      endif


      if (cell_selection .eqv. .FALSE.) then
        resample_in_voxel=.TRUE.
        if (verbose) write(*,*) 'Sampled position not acceptable - cell not in active list or zero intensity'
      endif
      
    endif

    !  check for the importance of the cell
    if(cell_importance.eq.0) then
      resample_in_voxel=.TRUE.
      if (verbose) write(*,*) 'Sampled position not acceptable - importance zero location'
    endif

    ! Check the material of the cell
    if(material.eq.0) then
      resample_in_voxel=.TRUE.
      if (verbose) write(*,*) 'Sampled position not acceptable - void'
    endif

    if ((resample_in_voxel .eqv. .FALSE.) .AND. (verbose)   ) write(*,*) 'Sampled position acceptable'


  end subroutine check_position
  !--------------------------------------------------------------------------

  !---------------------------------------------------------------------------------
  ! Subroutine to apply a source rotation and translation if applied
  !---------------------------------------------------------------------------------
  subroutine transform(x,y,z,u,v,w,cosines,translations, verbose)

    !**********************************************************************************************************
    ! Subroutine dummy variables 
    real(dknd), intent(in) :: cosines(1:9)
    real(dknd), intent(in) :: translations(1:3)
    real(dknd), intent(inout) :: x, y, z
    real(dknd), intent(inout) :: u, v, w
    logical, intent(in) :: verbose
    !**********************************************************************************************************

    
    !**********************************************************************************************************
    ! Local variables      
    real(dknd) :: xxx_initial,yyy_initial,zzz_initial
    real(dknd) :: uuu_initial,vvv_initial,www_initial
    !**********************************************************************************************************    
    
    xxx_initial = x
    yyy_initial = y
    zzz_initial = z
    
    uuu_initial = u
    vvv_initial = v
    www_initial = w

    x = cosines(1)*xxx_initial + cosines(2)*yyy_initial + cosines(3)*zzz_initial + translations(1)
    y = cosines(4)*xxx_initial + cosines(5)*yyy_initial + cosines(6)*zzz_initial + translations(2)
    z = cosines(7)*xxx_initial + cosines(8)*yyy_initial + cosines(9)*zzz_initial + translations(3)


    ! TE - I am not sure if we need to rotate the direction vectors as they are random
    ! AT - a future implementation may use directional biasing hence may not be random
    !      hence included for completeness
    u = cosines(1)*uuu_initial + cosines(2)*vvv_initial + cosines(3)*www_initial
    v = cosines(4)*uuu_initial + cosines(5)*vvv_initial + cosines(6)*www_initial
    w = cosines(7)*uuu_initial + cosines(8)*vvv_initial + cosines(9)*www_initial



    if (verbose) write(*,*) 'Transformed point : ',x,y,z

  end subroutine transform
  !--------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------
  ! function to sample a CDF  
  !-------------------------------------------------------------------------------------------

  integer function sample_CDF_bin(cdf_target, cdf_array)
    implicit none

    !**********************************************************************************************************
    ! Subroutine dummy variables 
    real(dknd), intent(in) :: cdf_target
    real(dknd), intent(in) :: cdf_array(:)
    !**********************************************************************************************************

    !**********************************************************************************************************
    ! Local variables     
    integer :: cdf_L, cdf_R, cdf_M
    !**********************************************************************************************************
    
    
    !this implements a binary search algorithm.

    cdf_L = 1
    cdf_R = size(cdf_array)
    do
       cdf_M = floor((cdf_L + cdf_R)/2.0)
       if (cdf_array(cdf_M) .le. cdf_target) cdf_L = cdf_M+1 
       if (cdf_array(cdf_M) .gt. cdf_target) cdf_R = cdf_M 
       if (cdf_R .eq. cdf_L)  exit  !success
    end do

    sample_CDF_bin = cdf_R

  end function sample_CDF_bin
  ! ---------------------------------------------


  !-----------------------------------------------------------------------------------
  ! function to construct a CDF  
  !-------------------------------------------------------------------------------------------
  subroutine construct_CDF(prob_array, cdf_array)
    implicit none
    
    !**********************************************************************************************************
    ! Subroutine dummy variables 
    real(dknd), intent(in) :: prob_array(:)
    real(dknd), intent(out) :: cdf_array(:)
    !**********************************************************************************************************
    
    !**********************************************************************************************************
    ! Local variables     
    integer :: i
    !**********************************************************************************************************
    
    
    cdf_array(1)= prob_array(1)
    ! sum the array
    do i=2,size(prob_array)
       cdf_array(i) = cdf_array(i-1) + prob_array(i) 
    end do

    !renormalise
       cdf_array = cdf_array/cdf_array(size(prob_array))

  end subroutine construct_CDF

end module cR2S_modules

