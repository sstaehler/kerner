!******************************************************************************
!
!    This file is part of:
!    MC Kernel: Calculating seismic sensitivity kernels on unstructured meshes
!    Copyright (C) 2016 Simon Staehler, Martin van Driel, Ludwig Auer
!
!    You can find the latest version of the software at:
!    <https://www.github.com/tomography/mckernel>
!
!    MC Kernel is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    MC Kernel is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with MC Kernel. If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

!=========================================================================================
module type_parameter
    use global_parameters,  only : sp, dp, pi, deg2rad, verbose, lu_out, master, testing
                                   
    use source_class,       only : src_param_type
    use kernel,             only : kernelspec_type
    use receiver_class,     only : rec_param_type
    use filtering,          only : filter_type
    use commpi,             only : pbroadcast_char, pbroadcast_int, pbroadcast_dble, &
                                   pbroadcast_dble_arr, pbroadcast_log, pbarrier, &
                                   pbroadcast_char_arr, pabort 
    use simple_routines,    only : lowtrim
    implicit none    

    type parameter_type
        type(src_param_type)                 :: source
        type(rec_param_type),   allocatable  :: receiver(:)
        type(kernelspec_type),  allocatable  :: kernel(:)
        type(filter_type),      allocatable  :: filter(:, :)
        integer                              :: nrec

        real(kind=dp)                        :: allowed_error
        real(kind=dp)                        :: allowed_relative_error = 1d-10
        real(kind=dp)                        :: damp_radius = -1d0
        real(kind=dp)                        :: time_for_dumping = 0

        character(len=512)                   :: fwd_dir
        character(len=512)                   :: bwd_dir
        character(len=512)                   :: source_file
        character(len=512)                   :: stf_file = 'stf.dat'
        character(len=512)                   :: receiver_file
        character(len=512)                   :: filter_file
        character(len=512)                   :: mesh_file
        character(len=512)                   :: mesh_file_type = 'abaqus'
        character(len=512)                   :: mesh_file_vert
        character(len=512)                   :: mesh_file_face
        character(len=512)                   :: output_file = 'kerner'
        character(len=512)                   :: hetero_file
        character(len=1)                     :: component
        character(len=512)                   :: int_type
        character(len=512)                   :: int_scheme = 'PARSEVAL'
        character(len=512)                   :: dump_type 
        character(len=32)                    :: fftw_plan = 'MEASURE'
        character(len=32)                    :: strain_type_fwd
        integer                              :: nsim_fwd, nsim_bwd
        integer                              :: nkernel
        integer                              :: nfilter
        integer                              :: nelems_per_task
        integer                              :: npoints_per_step
        integer                              :: max_iter
        integer                              :: strain_buffer_size
        integer                              :: displ_buffer_size
        logical                              :: plot_wavefields      = .false.
        logical                              :: parameters_read      = .false.
        logical                              :: receiver_read        = .false.
        logical                              :: source_read          = .false.
        logical                              :: filter_read          = .false.
        logical                              :: kernel_read          = .false.
        logical                              :: detailed_convergence = .false.
        logical                              :: deconv_stf           = .false.
        logical                              :: write_smgr           = .true.
        logical                              :: quasirandom          = .true.
        logical                              :: relative_kernel      = .true.
        logical                              :: int_over_volume      = .true.
        logical                              :: int_over_background  = .false.
        logical                              :: int_over_hetero      = .false.
        logical                              :: sort_mesh_elements   = .false.
        logical                              :: mask_src_rec         = .false.
        logical                              :: create_intermediate  = .false.
        logical                              :: parallel_read        = .false.
        contains
           procedure, pass                   :: read_parameters
           procedure, pass                   :: read_receiver
           procedure, pass                   :: read_source
           procedure, pass                   :: read_kernel
           procedure, pass                   :: read_filter
    end type

contains

!-----------------------------------------------------------------------------------------
subroutine read_parameters(this, input_file_in)
   use background_model, only       : backgroundmodel_type
   class(parameter_type)           :: this
   character(len=*), intent(in), optional :: input_file_in
   character(len=256)              :: input_file
   integer                         :: lu_inparam_basic, ioerr, narg
   character(len=256)              :: line
   character(len=256)              :: keyword, keyvalue
   logical                         :: temp_logical

   call pbarrier

   if (master) then

     if (present(input_file_in)) then
         input_file = input_file_in
     else
         narg = command_argument_count()
         if (narg<1) then
             print *, 'ERROR: Argument "parameter file" is missing'
             call pabort
         end if
         call get_command_argument(1, input_file)
     end if

     write(lu_out,'(" Reading ", A, "...")') trim(input_file)
     open(newunit=lu_inparam_basic, file=trim(input_file), status='old', &
          action='read', iostat=ioerr)
     if (ioerr /= 0) then
        print *, 'ERROR: Check input file ''', trim(input_file), '''! Is it still there?' 
        call pabort
     end if

     do
        read(lu_inparam_basic, fmt='(a256)', iostat=ioerr) line
        if (ioerr < 0) exit
        if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
     
        read(line,*, iostat=ioerr) keyword, keyvalue 
        if (ioerr < 0) cycle
      
        parameter_to_read : select case(trim(keyword))
        case('ALLOWED_ERROR')
           read(keyvalue, *) this%allowed_error

        case('ALLOWED_RELATIVE_ERROR')
           read(keyvalue, *) this%allowed_relative_error

        case('FWD_DIR')
           this%fwd_dir = keyvalue

        case('BWD_DIR')
           this%bwd_dir = keyvalue

        case('SRC_FILE')
           this%source_file = keyvalue

        case('STF_FILE')
           this%stf_file = keyvalue

        case('REC_FILE')
           this%receiver_file = keyvalue

        case('FILT_FILE')
           this%filter_file = keyvalue

        case('MESH_FILE_TYPE')
           this%mesh_file_type = keyvalue

        case('MESH_FILE_ABAQUS')
           this%mesh_file = keyvalue

        case('MESH_FILE_VERTICES')
           this%mesh_file_vert = keyvalue

        case('MESH_FILE_FACETS')
           this%mesh_file_face = keyvalue

        case('NO_SORT_MESH_ELEMENTS')
           read(keyvalue, *) temp_logical
           this%sort_mesh_elements = .not. temp_logical

        case('OUT_PREFIX')
           this%output_file = keyvalue

        case('DUMP_TYPE')
           this%dump_type = keyvalue

        case('HET_FILE')
           this%hetero_file = keyvalue

        case('STRAIN_BUFFER_SIZE')
           read(keyvalue, *) this%strain_buffer_size

        case('DISPL_BUFFER_SIZE')
           read(keyvalue, *) this%displ_buffer_size

        case('MAXIMUM_ITERATIONS')
           read(keyvalue, *) this%max_iter

        case('ELEMENTS_PER_TASK')
           read(keyvalue, *) this%nelems_per_task

        case('POINTS_PER_MC_STEP')
           read(keyvalue, *) this%npoints_per_step

        case('WRITE_DETAILED_CONVERGENCE')
           read(keyvalue, *) this%detailed_convergence

        case('USE_PSEUDORANDOM_NUMBERS')
           read(keyvalue, *) temp_logical
           this%quasirandom = .not. temp_logical

        case('KERNEL_FOR_ABSOLUTE_PERTURBATIONS')
           read(keyvalue, *) temp_logical
           this%relative_kernel = .not. temp_logical

        case('NO_INT_OVER_VOLUME')
           read(keyvalue, *) temp_logical
           this%int_over_volume = .not. temp_logical

        case('INT_OVER_BACKGROUND_MODEL')
           read(keyvalue, *) this%int_over_background

        case('INT_OVER_3D_HETEROGENEITIES')
           read(keyvalue, *) this%int_over_hetero

        case('NO_DECONVOLVE_STF')
           read(keyvalue, *) temp_logical
           this%deconv_stf = .not. temp_logical

        case('FFTW_PLAN')
           read(keyvalue, *) this%fftw_plan

        case('PLOT_WAVEFIELDS')
           read(keyvalue, *) this%plot_wavefields 

        case('INT_TYPE')
           this%int_type = keyvalue

        case('WRITE_SEISMOGRAMS')
           read(keyvalue, *) this%write_smgr

        case('INTEGRATION_SCHEME')
           read(keyvalue, *) this%int_scheme

        case('MASK_SOURCE_RECEIVER')
           read(keyvalue, *) this%mask_src_rec

        case('DAMP_RADIUS_SOURCE_RECEIVER')
           read(keyvalue, *) this%damp_radius

        case('CREATE_INTERMEDIATE')
           read(keyvalue, *) this%create_intermediate

        case('INTERMEDIATE_DUMP_TIME')
           read(keyvalue, *) this%time_for_dumping 
           ! Value is given in hours
           this%time_for_dumping = this%time_for_dumping * 3600

        case('PARALLEL_READING')
           read(keyvalue, *) this%parallel_read

        case default
           print *, 'Unknown parameter', trim(keyword)
           stop


        end select parameter_to_read
        if (verbose>0) write(lu_out, "('  ', A32,' ',A)")  keyword, trim(keyvalue)
     end do
     close(lu_inparam_basic)

  else
     write(lu_out,*) 'Waiting for master to distribute inparam_basic...'

  end if

  
  ! broadcast values to other processors
  call pbroadcast_char(this%fwd_dir, 0) 
  call pbroadcast_char(this%bwd_dir, 0) 
  call pbroadcast_dble(this%allowed_error, 0)
  call pbroadcast_dble(this%allowed_relative_error, 0)
  call pbroadcast_char(this%source_file, 0) 
  call pbroadcast_char(this%receiver_file, 0)
  call pbroadcast_char(this%filter_file, 0)
  call pbroadcast_char(this%stf_file, 0)
  call pbroadcast_char(this%mesh_file_type, 0)
  call pbroadcast_log(this%sort_mesh_elements, 0)
  call pbroadcast_char(this%output_file, 0)
  call pbroadcast_char(this%stf_file, 0)
  call pbroadcast_char(this%dump_type, 0)
  call pbroadcast_int(this%strain_buffer_size, 0)
  call pbroadcast_int(this%displ_buffer_size, 0)
  call pbroadcast_int(this%max_iter, 0)
  call pbroadcast_int(this%nelems_per_task, 0)
  call pbroadcast_int(this%npoints_per_step, 0)
  call pbroadcast_log(this%detailed_convergence, 0)
  call pbroadcast_log(this%deconv_stf, 0)
  call pbroadcast_log(this%write_smgr, 0)
  call pbroadcast_log(this%quasirandom, 0)
  call pbroadcast_log(this%relative_kernel, 0)
  call pbroadcast_log(this%int_over_volume, 0)
  call pbroadcast_log(this%int_over_background, 0)
  call pbroadcast_log(this%int_over_hetero, 0)
  call pbroadcast_char(this%fftw_plan, 0)
  call pbroadcast_log(this%plot_wavefields, 0)
  call pbroadcast_char(this%int_type, 0)
  call pbroadcast_char(this%hetero_file, 0)
  call pbroadcast_char(this%int_scheme, 0)
  call pbroadcast_log(this%parallel_read, 0)

  select case(lowtrim(this%mesh_file_type))
  case('abaqus') 
    call pbroadcast_char(this%mesh_file, 0)
  case('tetrahedral')
    call pbroadcast_char(this%mesh_file_vert, 0)
    call pbroadcast_char(this%mesh_file_face, 0)
  case default
    print *, 'ERROR: Unknown MESH_FILE_TYPE in ', trim(input_file)
    print *, '       Allowed values: ''abaqus'', ''tetrahedral'''
    call pabort(do_traceback=.false.)
  end select

  write(lu_out,*)
  call flush(lu_out)

  this%parameters_read = .true.

end subroutine read_parameters
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_source(this)
   class(parameter_type)          :: this
   real(kind=dp)                  :: Mij_dyncm(6), latd, lond, depth
   character(len=16)              :: junk
   character(len=16)              :: event_name
   integer, parameter             :: lu_source=1000

   if (.not.this%parameters_read) then
       write(6,*) 'General parameters not read. Call read_parameters before read_source!'
       stop
   end if

   call pbarrier

   write(lu_out,'(A)') '***************************************************************'
   if (master) then
       write(lu_out,*) 'Reading source from file ', trim(this%source_file)
       open(unit=lu_source, file=trim(this%source_file), status='old')
       read(lu_source,*) junk
       read(lu_source,*) junk, event_name
       read(lu_source,*) junk
       read(lu_source,*) junk
       read(lu_source,*) junk, latd
       read(lu_source,*) junk, lond  
       read(lu_source,*) junk, depth
       read(lu_source,*) junk, Mij_dyncm(1) !Mrr
       read(lu_source,*) junk, Mij_dyncm(2) !Mtt
       read(lu_source,*) junk, Mij_dyncm(3) !Mpp
       read(lu_source,*) junk, Mij_dyncm(4) !Mrt
       read(lu_source,*) junk, Mij_dyncm(5) !Mrp
       read(lu_source,*) junk, Mij_dyncm(6) !Mtp
       close(lu_source)

   else
      write(lu_out,*) 'Waiting for master to distribute source...'

   end if

   write(lu_out,'(A)') '***************************************************************'

   call pbroadcast_dble(latd, 0)
   call pbroadcast_dble(lond, 0)
   call pbroadcast_dble(depth, 0)
   call pbroadcast_dble_arr(Mij_dyncm, 0)

   call this%source%init(lat=latd, lon=lond, mij=Mij_dyncm*1.E-7, depth=depth)
   call this%source%read_stf(this%stf_file)

   write(lu_out,*)
   call flush(lu_out)

   this%source_read = .true.

end subroutine read_source
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_receiver(this)
   use simple_routines, only      : to_lower 
   use kernel,          only      : tabulate_kernels
   class(parameter_type)         :: this
   integer, parameter            :: lu_receiver = 1001
   integer                       :: irec, firstkernel, lastkernel
   integer                       :: ikernel, recnkernel, recistf
   real(kind=dp)                 :: reclatd, reclond
   character(len=16)             :: recname, trash(5), model_parameter, src_name
   character(len=80)             :: fmtstring

   if (.not.this%parameters_read) then
       write(6,*) 'General parameters not read. Call read_parameters before read_receiver!'
       call pabort
   end if

   if (.not.this%source_read) then
       write(6,*) 'Source parameters not read. Call read_source before read_receiver!'
       call pabort
   end if

   this%strain_type_fwd = 'straintensor_trace'

   call pbarrier

   if (master) then
       write(lu_out,*)'Reading receivers from file ', trim(this%receiver_file)
       open(unit=lu_receiver, file=trim(this%receiver_file), status='old')
       read(lu_receiver,*) this%nrec
       read(lu_receiver,*) this%component

   else
      write(lu_out,*) 'Waiting for master to distribute receivers...'

   end if

   call pbroadcast_int(this%nrec, 0)
   call pbroadcast_char(this%component, 0)


   fmtstring = '("  Using ", I5, " receivers")'
   write(lu_out, fmtstring) this%nrec

   fmtstring = '("  Kernel on component ", A)'
   write(lu_out, fmtstring) this%component

   allocate(this%receiver(this%nrec))

   firstkernel = 1
   lastkernel  = 0
   do irec = 1, this%nrec
      if (master) then
          read(lu_receiver, *) recname, reclatd, reclond, src_name, recnkernel
          recistf = this%source%get_src_index(src_name)
      end if
      call pbroadcast_char(recname, 0)
      call pbroadcast_dble(reclatd, 0)
      call pbroadcast_dble(reclond, 0)
      call pbroadcast_int(recnkernel, 0)
      call pbroadcast_int(recistf, 0)

      lastkernel = lastkernel + recnkernel

      call this%receiver(irec)%init(name        = recname        , &
                                    lat         = reclatd        , &
                                    lon         = reclond        , &
                                    component   = this%component , &
                                    nkernel     = recnkernel     , &  
                                    firstkernel = firstkernel    , &
                                    lastkernel  = lastkernel     , &
                                    istf        = recistf         )
      firstkernel = firstkernel + recnkernel

      if (master) then
          do ikernel = 1, recnkernel
             ! Just get the model parameter
             ! If just one kernel needs the full strain tensor, it has to 
             ! be read all the file. Otherwise, we would have to create separate
             ! buffers and what else. 
             ! TODO: Improve this
             read(lu_receiver, *) trash(1:5), model_parameter
             select case(model_parameter)
             case('lam', 'vp ', 'vph', 'vpv')

             case('vs ', 'rho', 'vsh', 'vsv', 'eta', 'phi', 'xi ', 'mu ')
               this%strain_type_fwd = 'straintensor_full'
             case default
               print *, 'Unknown model parameter ', model_parameter, &
                        'in kernel ', ikernel
               stop
             end select
                  
          end do
       end if


       call this%receiver(irec)%rotate_receiver( this%source%trans_rot_mat )
       
       fmtstring = '("  Receiver ''", A8, "'', lat:  ", F8.3, ", lon: ", F8.3)'
       write(lu_out, fmtstring) trim(recname), reclatd, reclond
       fmtstring = '("                       dist: ", F8.3, ", azi: ", F8.3)'
       write(lu_out, fmtstring) this%receiver(irec)%theta / deg2rad, &
                                this%receiver(irec)%phi / deg2rad

   end do

   call pbroadcast_char(this%strain_type_fwd, 0)

   if (master) then 
       close(lu_receiver)
   end if

   this%nkernel = lastkernel

   write(lu_out,*) ' In total ', this%nkernel, ' Kernels to calculate'

   this%receiver_read = .true.

   write(lu_out,*)
   call flush(lu_out)

end subroutine read_receiver
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_kernel(this, sem_data, filter)
   use readfields, only            : semdata_type
   use filtering,  only            : filter_type
   class(parameter_type)          :: this
   type(semdata_type), intent(in), optional :: sem_data
   type(filter_type), target, intent(in), optional  :: filter(:, :)
   integer                        :: irec, nfilter, istf
   integer                        :: ikernel, ifilter
   integer                        :: lu_receiver
   real(kind=dp)                  :: timewindow(2)
   character(len=4)               :: misfit_type
   character(len=16)              :: model_parameter
   character(len=32)              :: kernelname, filtername, kernel_shortname
   character(len=80)              :: fmtstring

   if (.not.this%parameters_read) then
       write(6,*) 'General parameters not read. Call read_parameters before read_kernel!'
       call pabort
   end if

   if (.not.this%source_read) then
       write(6,*) 'Source parameters not read. Call read_source before read_kernel!'
       call pabort
   end if

   if (.not.this%receiver_read) then
       write(6,*) 'Receiver parameters not read. Call read_receiver before read_kernel!'
       call pabort
   end if

   if (.not.this%filter_read) then
       write(6,*) 'Filter parameters not read. Call read_filter before read_kernel!'
       call pabort
   end if

   call pbarrier

   if (master) then
       write(lu_out,*)'Reading kernels from file ', trim(this%receiver_file)
       open(newunit=lu_receiver, file=trim(this%receiver_file), status='old')
       read(lu_receiver,*) 
       read(lu_receiver,*) this%component

   else
      write(lu_out,*) 'Waiting for master to distribute kernels...'

   end if

   allocate(this%kernel(this%nkernel))

   if ((.not.master).or.testing) then
       nfilter = size(filter, dim=1)

       if (this%deconv_stf) then
          do ifilter = 1, nfilter
            do istf = 1, this%source%nstf
              call filter(ifilter, istf)%add_stfs(sem_data%stf_fwd,         &
                                                  sem_data%dt,              &
                                                  sem_data%amplitude_fwd,   &
                                                  istf,                     &
                                                  this%source%stf(:,istf),  &
                                                  this%source%stf_dt(istf), &
                                                  this%source%stf_shift(istf))
            end do
          end do
       end if
   else
       nfilter = 0
   end if

   do irec = 1, this%nrec
      if (master) then
          read(lu_receiver, *) 
      end if
      
      ! Start with the assumption that this receiver needs no base kernels at all
      this%receiver(irec)%needs_basekernel = .false.

      fmtstring = '("  Receiver ", A, ", has ", I3, " kernel, first and last:", 2(I5))'
      write(lu_out,fmtstring) trim(this%receiver(irec)%name),  this%receiver(irec)%nkernel,   &
                              this%receiver(irec)%firstkernel, this%receiver(irec)%lastkernel

      do ikernel = this%receiver(irec)%firstkernel, this%receiver(irec)%lastkernel

         if (master) read(lu_receiver, *) kernel_shortname, filtername, & 
                                          misfit_type, timewindow, model_parameter

         call pbroadcast_char(model_parameter, 0)
         call pbroadcast_char(kernel_shortname, 0)
         call pbroadcast_char(filtername, 0)
         call pbroadcast_char(misfit_type, 0)
         call pbroadcast_dble_arr(timewindow, 0)
         call pbroadcast_char(model_parameter, 0)

         kernelname = trim(this%receiver(irec)%name)//'_'//trim(kernel_shortname)

         if ((.not.master).or.(testing)) then
             if (.not.(present(sem_data).and.present(filter))) then
                write(*,*) 'ERROR: Only master may call read_kernel without sem_data and filter'
                call pabort()
             end if
             do ifilter = 1, nfilter
                if (trim(filtername).eq.trim(filter(ifilter, 1)%name)) exit
             end do
             if (ifilter == nfilter + 1) then
                print *, 'Could not find filter ', trim(filtername), ', which was requested'
                print *, 'by kernel ', trim(kernelname)
                print *, 'Available filters: ', [(filter(ifilter, 1)%name, ifilter = 1, nfilter)]
                call pabort
             end if

             call this%kernel(ikernel)%init(name            = kernelname                ,  &
                                            time_window     = timewindow                ,  &
                                            filter          = filter(ifilter,              &
                                                                this%receiver(irec)%istf), &
                                            misfit_type     = misfit_type               ,  &
                                            model_parameter = model_parameter           ,  &
                                            seis            = sem_data%seis(:,irec)     ,  &
                                            dt              = sem_data%dt               ,  &
                                            deconv_stf      = this%deconv_stf           ,  &
                                            timeshift_fwd   = sem_data%timeshift_fwd    ,  &
                                            write_smgr      = this%write_smgr)

             ! Update which base kernels are needed for this receiver
             this%receiver(irec)%needs_basekernel =   &
                                  this%receiver(irec)%needs_basekernel &
                             .or. this%kernel(ikernel)%needs_basekernel

             select case(model_parameter)
             case('vp', 'vph', 'vpv', 'lam')
               this%receiver(irec)%t_last_p = max(this%receiver(irec)%t_last_p, &
                                                  timewindow(2))
             case default
               this%receiver(irec)%t_last_s = max(this%receiver(irec)%t_last_s, &
                                                  timewindow(2))
             end select
         else
             this%kernel(ikernel)%name = kernelname
         end if
        
         if (.not.testing) then
           ! Communicate the model parameter index back to the master
           call pbroadcast_int(this%kernel(ikernel)%model_parameter_index, 1)
           call pbroadcast_int(this%kernel(ikernel)%hetero_parameter_index, 1)
         end if


      end do ! ikernel
 
      ! TODO: This might allow a finer control over what has to be read. The final idea:
      ! If any kernel of this receiver needs more than a lambda base kernel,
      ! the whole strain tensor has to be read in. If only lambda, the trace is enough.
      ! However, this would require some involved changes in readfields.f90. 
      if (any(this%receiver(irec)%needs_basekernel(2:6))) then
         this%receiver(irec)%strain_type = 'straintensor_full'
         !this%strain_type_fwd = 'straintensor_full'
      else
         this%receiver(irec)%strain_type = 'straintensor_trace'
      end if

   end do ! irec
   if (master) close(lu_receiver)

   this%nkernel = this%receiver(this%nrec)%lastkernel

   write(lu_out,*)
   call flush(lu_out)

   this%kernel_read = .true.

end subroutine read_kernel
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_filter(this, nomega, df, nstf)
   class(parameter_type)     :: this
   integer, intent(in), optional        :: nomega
   real(kind=dp), intent(in), optional  :: df
   integer, intent(in), optional        :: nstf
   integer                              :: lu_filter
   integer                              :: ifilter, i
   character(len=32)                    :: filtername, filtertype
   character(len=80)                    :: fmtstring
   real(kind=dp)                        :: freqs(4)

   if (.not.this%parameters_read) then
       write(6,*) 'General parameters not read. Call read_parameters before read_filter!'
       call pabort
   end if

   write(lu_out,*)'Reading filters from file ', trim(this%filter_file)
   if (master) then
       open(newunit=lu_filter, file=trim(this%filter_file), status='old')
       read(lu_filter, *) this%nfilter

   else
      write(lu_out,*) 'Waiting for master to distribute filters...'

   end if

   call pbroadcast_int(this%nfilter, 0)

   if (present(nstf)) then
     allocate(this%filter(this%nfilter, nstf))
   else
     allocate(this%filter(this%nfilter, 1))
   end if
  
   fmtstring = '("  Creating ", I5, " filters")'
   write(lu_out,fmtstring) this%nfilter

   fmtstring = '("  Creating filter ", A, " of type ", A/, "   freqs: ", 4(F8.3))'
   do ifilter = 1, this%nfilter
      if (master) read(lu_filter, *) filtername, filtertype, freqs
      call pbroadcast_char(filtername, 0)
      call pbroadcast_char(filtertype, 0)
      call pbroadcast_dble_arr(freqs, 0)

      write(lu_out,fmtstring) trim(filtername), trim(filtertype), freqs
      
      if((.not.master).or.(testing)) then
          if (.not.(present(df).and.present(nomega))) then
             write(*,*) 'ERROR: Only master may call read_filter without nomega and df'
             call pabort()
          end if
          if (present(nstf)) then
            do i = 1, nstf
              call this%filter(ifilter, i)%create(filtername, df, nomega, &
                                                  filtertype, freqs)
            end do
          else
            call this%filter(ifilter, 1)%create(filtername, df, nomega, &
                                                filtertype, freqs)
          end if
      end if

   end do

   if (master) close(lu_filter)

   write(lu_out,*)
   call flush(lu_out)

   this%filter_read = .true.

end subroutine read_filter
!-----------------------------------------------------------------------------------------

end module type_parameter
!=========================================================================================
