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
module slave_mod

  use global_parameters,           only: sp, dp, long, pi, deg2rad, &
                                         init_random_seed, myrank, lu_out
  use work_type_mod
# ifndef include_mpi
  use mpi
# endif
  implicit none

# ifdef include_mpi
  include 'mpif.h'
# endif

  type slave_result_type  
    real(kind=dp), allocatable :: kernel_values(:,:,:)
    real(kind=dp), allocatable :: kernel_variance(:,:,:)
    integer,       allocatable :: niterations(:,:)
    real(kind=dp), allocatable :: computation_time(:)
    real(kind=dp), allocatable :: model(:,:,:) !< (nmodel_parameter, nbasisfuncs_per_elem, nelements)
                                               !! Model parameters in the order as defined in backgroundmodel:
                                               !! vp, vs, rho, vph, vpv, vsh, vsv, eta, phi, xi, lam, mu
    real(kind=dp), allocatable :: hetero_model(:,:,:) !< (nmodel_parameter_hetero, nbasisfuncs_per_elem, nelements)
                                               !! dlnvp, dlnvs, dlnrho
    real(kind=sp), allocatable :: fw_field(:,:,:,:)
    real(kind=sp), allocatable :: bw_field(:,:,:,:)
    real(kind=sp), allocatable :: conv_field(:,:,:,:)
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine do_slave() 
  use global_parameters,           only: DIETAG, id_mpi, id_read_params, id_init_fft, &
                                         id_int_hetero
  use inversion_mesh,              only: inversion_mesh_data_type
  use readfields,                  only: semdata_type
  use heterogeneities,             only: hetero_type
  use type_parameter,              only: parameter_type
  use fft,                         only: rfft_type, taperandzeropad
  use clocks_mod,                  only: tick

  implicit none
  type(parameter_type)                :: parameters
  type(inversion_mesh_data_type)      :: inv_mesh
  type(semdata_type)                  :: sem_data
  type(rfft_type)                     :: fft_data
  type(hetero_type)                   :: het_model
  type(slave_result_type)             :: slave_result

  integer                             :: ndumps, ntimes, nomega
  integer                             :: ikernel, ierror
  integer(kind=long)                  :: iclockold
  integer                             :: mpistatus(MPI_STATUS_SIZE)
  real(kind=dp)                       :: df
  integer                             :: itask
  character(len=255)                  :: fmtstring
  integer                             :: nptperstep
  integer(kind=long)                  :: iclock, iclock_ref, ticks_per_sec
  real(kind=dp)                       :: waittime

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Read input files for parameters, source and receivers'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)

  iclockold = tick()

  call parameters%read_parameters()
  call parameters%read_source()
  call parameters%read_receiver()

  nptperstep = parameters%npoints_per_step

  ! Master and slave part ways here for some time. 
  ! Master reads in the inversion mesh, slaves initialize the FFT

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Initialize and open AxiSEM wavefield files'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)

  call sem_data%set_params(parameters%fwd_dir,     &
                           parameters%bwd_dir,     &
                           parameters%strain_buffer_size, & 
                           parameters%displ_buffer_size, & 
                           parameters%strain_type_fwd,    &
                           parameters%source%depth,     &
                           parameters%npoints_per_step,     &
                           parameters%parallel_read)

  call sem_data%open_files()
  call sem_data%read_meshes()

  call sem_data%load_seismogram(parameters%receiver, parameters%source)

  ndumps = sem_data%ndumps

  iclockold = tick(id=id_read_params, since=iclockold)

  if (parameters%int_over_hetero) then
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Initialize heterogeneity structure'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)

    call het_model%load_het_rtpv(parameters%hetero_file)
    call het_model%build_hetero_kdtree()

    iclockold = tick(id=id_int_hetero, since=iclockold)
  end if

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Initialize FFT'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)

  call fft_data%init(ndumps, sem_data%get_ndim(), nptperstep, sem_data%dt, &
                     fftw_plan=parameters%fftw_plan)

  ntimes = fft_data%get_ntimes()
  nomega = fft_data%get_nomega()
  df     = fft_data%get_df()
  fmtstring = '(A, I8, A, I8)'
  write(lu_out,fmtstring) '  ntimes: ',  ntimes,     '  , nfreq: ', nomega
  fmtstring = '(A, F8.3, A, F8.3, A)'
  write(lu_out,fmtstring) '  dt:     ', sem_data%dt, ' s, df:    ', df*1000, ' mHz'

  iclockold = tick(id=id_init_fft, since=iclockold)

  ! Master and slave synchronize again

  iclockold = tick()
  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Define filters'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)

  call parameters%read_filter(nomega, df, parameters%source%nstf)

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Define kernels'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)

  call parameters%read_kernel(sem_data, parameters%filter)

  iclockold = tick(id=id_read_params, since=iclockold)

  itask = 0

  call system_clock(count=iclock_ref, count_rate=ticks_per_sec)
  ! Receive new tasks from master until the DIETAG is received and then exit the loop
  receive_tasks: do 
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A)') ' Waiting for tasks from master rank'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)


    ! Receive a message from the master
    iclockold = tick()
    call MPI_Recv(wt, 1, wt%mpitype, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, ierror)
    call system_clock(count=iclock)
    waittime = real(iclock-iclock_ref, kind=dp) / real(ticks_per_sec, kind=dp)

    ! Check the tag of the received message. If no more work to do, exit loop
    ! and return to main program
    iclockold = tick(id=id_mpi, since=iclockold)
    if (mpistatus(MPI_TAG) == DIETAG) exit

    itask = itask + 1
    write(lu_out,'(A)') '***************************************************************'
    write(lu_out,'(A, I6, A, F8.3, A)') ' Received task # ', itask, ' after', waittime, &
      ' sec, going to work'
    write(lu_out,'(A)') '***************************************************************'
    call flush(lu_out)


    call inv_mesh%initialize_mesh(wt%ielement_type, wt%vertices, &
                                  wt%connectivity, wt%nbasisfuncs_per_elem)


    slave_result = slave_work(parameters, sem_data, inv_mesh, fft_data, het_model)

    wt%kernel_values(:,:,:)    = slave_result%kernel_values(:,:,:)
    wt%kernel_variance(:,:,:)  = slave_result%kernel_variance(:,:,:)
    wt%niterations(:,:)        = slave_result%niterations(:,:)
    wt%computation_time(:)     = slave_result%computation_time(:)
    wt%model(:,:,:)            = slave_result%model(:,:,:)
    wt%hetero_model(:,:,:)     = slave_result%hetero_model(:,:,:)
    wt%fw_field(:,:,:,:)       = slave_result%fw_field(:,:,:,:)
    wt%bw_field(:,:,:,:)       = slave_result%bw_field(:,:,:,:)
    wt%conv_field(:,:,:,:)     = slave_result%conv_field(:,:,:,:)

    ! Finished writing my part of the mesh
    call inv_mesh%freeme

    ! Send the result back
    call flush(lu_out)
    iclockold = tick()
    call system_clock(count=iclock_ref)
    call MPI_Send(wt, 1, wt%mpitype, 0, 0, MPI_COMM_WORLD, ierror)
    iclockold = tick(id=id_mpi, since=iclockold)

  end do receive_tasks

  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A, I6, A)') ' All work done. Did ', itask, ' tasks in total'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)


  write(lu_out,'(A)')
  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Free memory of Kernelspecs'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)
  do ikernel = 1, parameters%nkernel
    call parameters%kernel(ikernel)%freeme() 
  end do

  write(lu_out,'(A)')
  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Free memory of inversion mesh datatype'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)
  call inv_mesh%freeme

  write(lu_out,'(A)')
  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Close AxiSEM wavefield files'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)
  call sem_data%close_files()

  write(lu_out,'(A)')
  write(lu_out,'(A)') '***************************************************************'
  write(lu_out,'(A)') ' Free memory of FFT datatype'
  write(lu_out,'(A)') '***************************************************************'
  call flush(lu_out)
  call fft_data%freeme()

end subroutine do_slave
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function slave_work(parameters, sem_data, inv_mesh, fft_data, het_model) result(slave_result)

  use global_parameters,           only: sp, dp, pi, deg2rad, myrank, &
                                         id_fwd, id_bwd, id_fft, id_mc, id_filter_conv, &
                                         id_inv_mesh, id_kernel, id_init, id_int_model, &
                                         id_element, id_int_hetero

  use inversion_mesh,              only: inversion_mesh_data_type
  use readfields,                  only: semdata_type, dampen_field
  use type_parameter,              only: parameter_type
  use fft,                         only: rfft_type, taperandzeropad
  use filtering,                   only: timeshift_type
  use mc_integration,              only: integrated_type, allallconverged, allisconverged
  use kernel,                      only: calc_basekernel, calc_physical_kernels
  use background_model,            only: backgroundmodel_type, nmodel_parameters
  use heterogeneities,             only: hetero_type, nmodel_parameters_hetero         
  use clocks_mod,                  only: tick, get_clock, reset_clock
  use simple_routines,             only: mult3d_1d, cumsum_trapezoidal

  type(inversion_mesh_data_type), intent(in)    :: inv_mesh
  type(parameter_type),           intent(in)    :: parameters
  type(semdata_type),             intent(in)    :: sem_data
  type(rfft_type),                intent(in)    :: fft_data
  type(hetero_type),    intent(in), optional    :: het_model
  type(slave_result_type)                       :: slave_result

  type(timeshift_type)                :: timeshift_fwd, timeshift_bwd
  type(integrated_type), allocatable  :: int_kernel(:), int_model(:), int_hetero(:)
  type(backgroundmodel_type)          :: bg_model

  real(kind=dp),    allocatable       :: fw_field(:,:,:)
  real(kind=dp),    allocatable       :: bw_field(:,:,:)
  complex(kind=dp), allocatable       :: fw_field_fd(:,:,:)
  complex(kind=dp), allocatable       :: bw_field_fd(:,:,:)
  complex(kind=dp), allocatable       :: conv_field_fd(:,:), conv_field_fd_filt(:,:)
  real(kind=dp),    allocatable       :: conv_field(:,:)
  real(kind=dp),    allocatable       :: random_points(:,:) 
  real(kind=dp),    allocatable       :: weights(:)
  real(kind=dp),    allocatable       :: kernelvalue_basekers(:,:,:)
  real(kind=dp),    allocatable       :: kernelvalue_weighted(:,:,:) 
  real(kind=dp),    allocatable       :: kernelvalue_physical(:,:) ! this is linear combs of ...
  integer,          allocatable       :: niterations(:), kernel_list(:)

  ! This block is only needed for wavefield plotting and only 
  ! allocated, if it is enabled
  real(kind=dp),    allocatable       :: fw_field_filt(:,:,:)
  real(kind=dp),    allocatable       :: bw_field_filt(:,:,:)
  real(kind=dp),    allocatable       :: fw_field_out(:,:,:)
  real(kind=dp),    allocatable       :: bw_field_out(:,:,:)
  real(kind=dp),    allocatable       :: conv_field_out(:,:), conv_field_temp(:,:,:,:)
  complex(kind=dp), allocatable       :: fw_field_fd_filt(:,:,:)
  complex(kind=dp), allocatable       :: bw_field_fd_filt(:,:,:)

  real(kind=dp)                       :: time_in_element
  real(kind=dp)                       :: volume
  real(kind=dp)                       :: norm_fields
  real(kind=dp)                       :: planet_radius

  character(len=256)                  :: fmtstring
  integer                             :: ielement, irec, ikernel, ibasisfunc, ibasekernel
  integer                             :: nptperstep, ndumps, ntimes, nomega, nelements
  integer                             :: nbasisfuncs_per_elem, nkernel
  integer                             :: nbasekernels = 6
  integer(kind=long)                  :: iclockold, iclockold_element
  integer                             :: ndim
  integer, parameter                  :: taper_length = 10      !< This is the bare minimum. It does 
                                                                !! not produce artifacts yet, at least 
                                                                !! no obvious ones

  iclockold = tick()

  planet_radius = sem_data%fwd(1)%planet_radius

  ndim = sem_data%get_ndim()
  nptperstep = parameters%npoints_per_step
  ndumps = sem_data%ndumps
  ntimes = fft_data%get_ntimes()
  nomega = fft_data%get_nomega()
  nbasisfuncs_per_elem = inv_mesh%nbasisfuncs_per_elem
  nelements = inv_mesh%get_nelements()
  nkernel = parameters%nkernel

  allocate(slave_result%kernel_values(nkernel, nbasisfuncs_per_elem, nelements))
  allocate(slave_result%kernel_variance(nkernel, nbasisfuncs_per_elem, nelements))
  allocate(slave_result%niterations(nkernel, nelements))
  allocate(slave_result%computation_time(nelements))
  allocate(slave_result%model(nmodel_parameters, nbasisfuncs_per_elem, nelements))
  allocate(slave_result%hetero_model(nmodel_parameters_hetero, nbasisfuncs_per_elem, nelements))

  if (parameters%plot_wavefields) then
    allocate(slave_result%fw_field(ndumps, ndim, nkernel, nelements)) 
    allocate(slave_result%bw_field(ndumps, ndim, nkernel, nelements))
    allocate(slave_result%conv_field(ndumps, 1,  nkernel, nelements))

    allocate(fw_field_out(ndumps, ndim, nkernel))
    allocate(bw_field_out(ndumps, ndim, nkernel))
    allocate(conv_field_out(ndumps, nkernel))
    allocate(conv_field_temp(ndumps, nptperstep, nbasekernels, nkernel))
    allocate(fw_field_filt(ntimes, ndim, nptperstep))
    allocate(bw_field_filt(ntimes, ndim, nptperstep))
    allocate(fw_field_fd_filt(nomega, ndim, nptperstep))
    allocate(bw_field_fd_filt(nomega, ndim, nptperstep))
  else
    allocate(slave_result%fw_field(1,1,1,1))
    allocate(slave_result%bw_field(1,1,1,1))
    allocate(slave_result%conv_field(1,1,1,1))
  end if

  slave_result%kernel_values    = 0 
  slave_result%kernel_variance  = 0  
  slave_result%niterations      = 0  
  slave_result%computation_time = 0 
  slave_result%model            = 0 
  slave_result%hetero_model     = 0 
  slave_result%fw_field         = 0 
  slave_result%bw_field         = 0 
  slave_result%conv_field       = 0  


  allocate(fw_field(ndumps, ndim, nptperstep))
  allocate(bw_field(ndumps, ndim, nptperstep))

  allocate(fw_field_fd       (nomega, ndim, nptperstep))
  allocate(bw_field_fd       (nomega, ndim, nptperstep))

  allocate(conv_field        (ntimes, nptperstep))
  allocate(conv_field_fd     (nomega, nptperstep))
  allocate(conv_field_fd_filt(nomega, nptperstep))

  allocate(kernelvalue_basekers(nptperstep, parameters%nkernel, nbasekernels))
  allocate(kernelvalue_weighted(nptperstep, parameters%nkernel, nbasekernels))
  allocate(kernelvalue_physical(nptperstep, parameters%nkernel))

  volume = 0.0

  allocate(niterations(parameters%nkernel))

  if (.not.parameters%deconv_stf) then
    call timeshift_fwd%init_ts(fft_data%get_f(), sem_data%timeshift_fwd)
    call timeshift_bwd%init_ts(fft_data%get_f(), sem_data%timeshift_bwd)
  end if

  !    write(*,*) '***************************************************************'
  !    write(*,*) ' Start loop over elements'
  !    write(*,*) '***************************************************************'

  allocate(random_points(3, nptperstep))
  allocate(weights(nptperstep))
  allocate(int_kernel(nbasisfuncs_per_elem))

  if (parameters%int_over_background) then
    allocate(int_model(nbasisfuncs_per_elem))
  end if


  if (parameters%int_over_hetero) then
    allocate(int_hetero(nbasisfuncs_per_elem))
  end if

  iclockold = tick(id=id_init, since=iclockold)

  loop_elements: do ielement = 1, nelements 

    write(lu_out, '(A,I10,A)') 'Element: ', ielement, ':'

    ! Set element-specific clock to zero and start it again
    call reset_clock(id_element)
    iclockold_element = tick(id=id_element)

    ! Get volume of current element
    if (parameters%int_over_volume) then
      iclockold = tick()
      volume = inv_mesh%get_volume(ielement)
      iclockold = tick(id=id_inv_mesh)
    else
      volume = 1.d0
    end if

    ! Output arrays for wavefield plotting
    if (parameters%plot_wavefields) then
      fw_field_out = 0.0d0
      bw_field_out = 0.0d0
      conv_field_out = 0.0d0
      conv_field_temp = 0.0d0
    end if

    ! Background Model Integration 
    if (parameters%int_over_background) then
      iclockold = tick()
      int_model = integrate_1d_model(sem_data, inv_mesh, ielement)
      iclockold = tick(id=id_int_model, since=iclockold)
    end if

    ! Heterogeneity Model Integration 
    ! Calculate integrated heterogeneity structure for this element, optional
    if (parameters%int_over_hetero) then
      iclockold = tick()
      int_hetero = integrate_3d_model(het_model, inv_mesh, ielement)
      iclockold = tick(id=id_int_hetero, since=iclockold)
    end if

    ! Check, whether source is inside element. If so, do not do Monte Carlo integration
    if (inv_mesh%point_in_element(ielement, &
                                  parameters%source%get_r(planet_radius)) &
        .and.parameters%mask_src_rec ) then

        write(lu_out, '(A)') ' Contains source, not calculating kernel' 

        do ibasisfunc = 1, nbasisfuncs_per_elem          
          ! Model should be communicated even for the masked elements
          if (parameters%int_over_background) then
            slave_result%model(:, ibasisfunc, ielement) = int_model(ibasisfunc)%getintegral()
            call int_model(ibasisfunc)%freeme()
          end if
          if (parameters%int_over_hetero) then
            slave_result%hetero_model(:, ibasisfunc, ielement) = int_hetero(ibasisfunc)%getintegral()
            call int_hetero(ibasisfunc)%freeme()              
          end if
        end do

      else

        !> Sensitivity Kernel Integration <
        !  Initialize basis kernel Monte Carlo integrals for current element
        iclockold = tick()
        do ibasisfunc = 1, nbasisfuncs_per_elem
          ! Allowed absolute error is multiplied by volume to ensure that convergence
          ! is not affected by element size. It should just be a function of the kernel value
          call int_kernel(ibasisfunc)%initialize_montecarlo(parameters%nkernel,               &
                                                            volume,                           &
                                                            parameters%allowed_error          &
                                                             * volume,                        &
                                                            parameters%allowed_relative_error)
        end do

        niterations = 0

        loop_all_converged: do while (.not.allallconverged(int_kernel)) ! Beginning of Monte Carlo loop

          iclockold = tick(id=id_mc)

          ! Generate random points
          random_points = inv_mesh%generate_random_points(ielement, nptperstep, &
                                                          parameters%quasirandom)
          iclockold = tick(id=id_inv_mesh)

          ! Stop MC integration in this element after max_iter iterations
          if (any(niterations==parameters%max_iter)) exit 

          ! Load forward field
          iclockold = tick()
          fw_field = sem_data%load_fw_points( random_points, parameters%source, & 
                                              model = bg_model)
          ! Integrate time series once. The SEM simulation is run with a Gaussian source term,
          ! which means the wavefield has the form of velocities/strain rates. For the forward
          ! field, we need the displacement response to a Heaviside source term.
          fw_field = cumsum_trapezoidal(fw_field, sem_data%dt)

          iclockold = tick(id=id_fwd, since=iclockold)

          ! FFT of forward field
          call fft_data%rfft( taperandzeropad(fw_field, ntimes, taper_length), fw_field_fd )
          iclockold = tick(id=id_fft, since=iclockold)


          ! Timeshift forward field
          ! Not necessary, if STF is deconvolved
          if (.not.parameters%deconv_stf) then
            call timeshift_fwd%apply(fw_field_fd)
            iclockold = tick(id=id_filter_conv, since=iclockold)
          end if

          ! Loop over all receivers in receiver input file
          loop_receivers: do irec = 1, parameters%nrec

            if (inv_mesh%point_in_element(ielement, parameters%receiver(irec)%r) &
                .and.parameters%mask_src_rec) then

              do ikernel = parameters%receiver(irec)%firstkernel, &
                  parameters%receiver(irec)%lastkernel 

                ! Set scalar misfit base kernels to zero
                kernelvalue_basekers(:,ikernel,:) = 0.d0
                iclockold = tick(id=id_kernel, since=iclockold)

                niterations(ikernel) = niterations(ikernel) + 1
              end do ! ikernel


            elseif (travel_time_larger_Tmax(source   = parameters%source,                &
                                            receiver = parameters%receiver(irec),        &
                                            coordinates = inv_mesh%get_center(ielement), &
                                            v_max    = sem_data%get_v_max(),             &
                                            planet_radius = planet_radius)               &
                    .and. (.not.parameters%plot_wavefields)) then
              ! Test whether the minimum travel time from source to element 
              ! midpoint to receiver is smaller than the end of the latest 
              ! time window on this receiver.
              ! If wavefield plotting is switched on, no elements should be masked,
              ! since we want to see the wavefield everywhere.
              do ikernel = parameters%receiver(irec)%firstkernel, &
                  parameters%receiver(irec)%lastkernel 

                ! Set scalar misfit base kernels to zero
                kernelvalue_basekers(:,ikernel,:) = 0.d0
                iclockold = tick(id=id_kernel, since=iclockold)

                niterations(ikernel) = niterations(ikernel) + 1
              end do ! ikernel

            else
              !Receiver not in element
              
              ! List of kernels for this receiver
              if (allocated(kernel_list)) deallocate(kernel_list)
              allocate(kernel_list(parameters%receiver(irec)%nkernel))
              kernel_list = [(ikernel, ikernel = parameters%receiver(irec)%firstkernel,     &
                                                 parameters%receiver(irec)%lastkernel) ]

              ! If all kernels of this receiver are converged, cycle to next receiver
              if (allallconverged(int_kernel, kernel_list)) cycle
              iclockold = tick(id=id_mc, since=iclockold)

              ! Load backward field             
              bw_field = sem_data%load_bw_points( random_points, parameters%receiver(irec) )
              iclockold = tick(id=id_bwd, since=iclockold)

              ! FFT of backward field
              call fft_data%rfft( taperandzeropad(bw_field, ntimes, taper_length), bw_field_fd )
              iclockold = tick(id=id_fft, since=iclockold)

              ! Timeshift backward field
              ! Not necessary, if STF is deconvolved
              if (.not.parameters%deconv_stf) then
                call timeshift_bwd%apply(bw_field_fd)
                iclockold = tick(id=id_filter_conv, since=iclockold)
              end if



              ! Calculate the scalar kernel in the basic parameters ("base kernels")
              ! (lambda, mu, rho, a, b, c)  
              loop_basekernels: do ibasekernel = 1, nbasekernels
                ! Check whether any actual kernel on this receiver needs this base kernel
                if (.not.parameters%receiver(irec)%needs_basekernel(ibasekernel)) cycle

                ! Calculate waveform base kernel numero ibasekernel
                ! Mainly convolution of the right components of the fields.
                conv_field_fd = calc_basekernel(ibasekernel,                &
                                                parameters%strain_type_fwd, &
                                                parameters%strain_type_fwd, &
                                                fw_field_fd, bw_field_fd)
                iclockold = tick(id=id_filter_conv, since=iclockold)

                ! Calculate scalar base kernels for kernels of one receiver
                loop_kernels:                                        &
                 do ikernel = parameters%receiver(irec)%firstkernel, &
                              parameters%receiver(irec)%lastkernel 

                  ! Check whether kernel (ikernel) actually needs the base kernel 
                  ! (ibasekernel), otherwise cycle
                  if (.not.parameters%kernel(ikernel)%needs_basekernel(ibasekernel)) cycle

                  ! If this kernel is already converged, go to the next one
                  if (allisconverged(int_kernel, ikernel)) cycle
                  iclockold = tick(id=id_mc, since=iclockold)

                  ! Apply Filter 
                  conv_field_fd_filt = parameters%kernel(ikernel)%apply_filter(conv_field_fd)
                  iclockold = tick(id=id_filter_conv, since=iclockold)

                  ! Apply Inverse FFT
                  call fft_data%irfft(conv_field_fd_filt, conv_field)
                  iclockold = tick(id=id_fft, since=iclockold)

                  ! Calculate scalar misfit base kernels from convolved time traces
                  kernelvalue_basekers(:,ikernel,ibasekernel) =         &
                    parameters%kernel(ikernel)%calc_misfit_kernel(conv_field, parameters%int_scheme)


                  ! If the wavefields are to be plotted, it needs to be filtered and iFTd here.
                  if (parameters%plot_wavefields) then
                    conv_field_temp(:, :, ibasekernel, ikernel) =  & 
                      conv_field_temp(:, :, ibasekernel, ikernel) + conv_field(1:ndumps, :)
                  end if

                  iclockold = tick(id=id_kernel, since=iclockold)

                  niterations(ikernel) = niterations(ikernel) + 1
                end do loop_kernels

              end do loop_basekernels

              if (parameters%plot_wavefields) then
                ! Since the fwd and bwd wavefield do not depend on the basekernel type,
                ! this can and should be done outside of the basekernel loop
                do ikernel = parameters%receiver(irec)%firstkernel, &
                              parameters%receiver(irec)%lastkernel 

                  !Test of planar wave with lambda=1000km and f = 1/50s
                  !do idump = 1, ndumps
                  !  fw_field_out(idump,:,ikernel) = &
                  !    sum(cos(2*pi* (random_points(1,:)/1e6 + idump*sem_data%dt / (50.d0)) ))
                  !end do

                  fw_field_fd_filt = parameters%kernel(ikernel)%filter%apply(fw_field_fd, kind='fwd')
                  call fft_data%irfft(fw_field_fd_filt, fw_field_filt)

                  bw_field_fd_filt = parameters%kernel(ikernel)%filter%apply(bw_field_fd, kind='bwd')
                  call fft_data%irfft(bw_field_fd_filt, bw_field_filt)

                  ! Sum over all MC points in this iteration
                  ! TODO: Introduce weighting of MC point contributions
                  fw_field_out(:, :, ikernel) =  & 
                    fw_field_out(:, :, ikernel) + sum(fw_field_filt(1:ndumps, :, :), 3)
                  bw_field_out(:, :, ikernel) =  & 
                    bw_field_out(:, :, ikernel) + sum(bw_field_filt(1:ndumps, :, :), 3) 

                  ! Compute the waveform kernels for the actual physical parameters of interest
                  conv_field_out(:, ikernel) = sum(calc_physical_kernels( &
                    parameters%kernel(ikernel)%model_parameter,       &
                    conv_field_temp(:, :, :, ikernel),                &
                    bg_model = bg_model,                              &
                    relative_kernel = parameters%relative_kernel), 2)

                  ! Multiply the waveform kernel with the seismogram for this receiver
                  ! to mask the "coda sensitivity" in regions, where almost no energy arrives
                  ! This results in very strange waveform kernels, so only use it, if
                  ! you know what you're doing.
                  ! conv_field_out(:, ikernel) = conv_field_out(:, ikernel) * &
                  !   parameters%kernel(ikernel)%seis * parameters%kernel(ikernel)%normalization


                end do
              end if

              ! Dampen kernel values, if the points are close to receiver or source
              call dampen_field(kernelvalue_basekers, random_points,             &
                                parameters%receiver(irec)%get_r(planet_radius),  &
                                parameters%damp_radius)
              call dampen_field(kernelvalue_basekers, random_points,             &
                                parameters%source%get_r(planet_radius),          &
                                parameters%damp_radius)

              iclockold = tick(id=id_mc, since=iclockold)

              ! Deallocate list of kernels for this receiver, since the next receiver 
              ! will have another number of kernels
              deallocate(kernel_list)

            end if ! Receiver is inside element?

        end do loop_receivers

        iclockold = tick(id=id_mc, since=iclockold)

        ! Check for convergence           
        ! For volumetric basis functions, this loop has just one iteration,
        ! for vertex-based basis functions, nvertices_per_elem iterations
        do ibasisfunc = 1, nbasisfuncs_per_elem
          iclockold = tick()

          weights = inv_mesh%weights(ielement, ibasisfunc, random_points)                 
          kernelvalue_weighted = mult3d_1d(kernelvalue_basekers, weights)

          ! Compute weighted base kernels
          do ikernel = 1, parameters%nkernel

            ! Compute the kernels for the actual physical parameters of interest
            kernelvalue_physical(:, ikernel) = calc_physical_kernels(        &
              parameters%kernel(ikernel)%model_parameter, &
              kernelvalue_weighted(:, ikernel, :),        &
              bg_model = bg_model,                        &
              relative_kernel = parameters%relative_kernel)

          end do

          ! Check if the summed kernels for each basis function have converged
          iclockold = tick(id=id_kernel, since=iclockold)
          call int_kernel(ibasisfunc)%check_montecarlo_integral(kernelvalue_physical)
          iclockold = tick(id=id_mc, since=iclockold)
        end do        

        ! Print convergence info
        if (parameters%detailed_convergence) then
          write(fmtstring,"('(I6, I6, ES10.3, A, I0.5, A, I0.5, ', I4, '(ES11.3), A, ', I4, '(ES10.3))')") &
            parameters%nkernel, parameters%nkernel 
          write(lu_out,fmtstring) myrank, ielement, inv_mesh%get_volume(ielement), ' Converged? ', &
                                  int_kernel(1)%countconverged(), '/', parameters%nkernel,             &
                                  int_kernel(1)%getintegral(), ' +- ', sqrt(int_kernel(1)%getvariance())
          call flush(lu_out)
        end if

      end do loop_all_converged ! End of Monte Carlo loop

      ! Print number of iterations
      if (any(niterations(:)==parameters%max_iter)) then
        fmtstring = "(' Max number of iterations reached. ', I6 , ' kernels were not converged.')"
        write(lu_out, fmtstring) parameters%nkernel - int_kernel(1)%countconverged()
      else
        fmtstring = "(' All kernels converged after', I8, ' iterations')"
        write(lu_out, fmtstring) maxval(niterations(:))
      end if

      ! Get time spent in this element
      iclockold_element = tick(id=id_element, since = iclockold_element)
      call get_clock(id = id_element, total_time = time_in_element) 


      ! Fill slave_result object, which is later passed over to master
      slave_result%computation_time(ielement)    = time_in_element 
      slave_result%niterations(:,ielement)       = niterations(:)
      do ibasisfunc = 1, nbasisfuncs_per_elem          
        slave_result%kernel_values(:, ibasisfunc, ielement)   = int_kernel(ibasisfunc)%getintegral()
        slave_result%kernel_variance(:, ibasisfunc, ielement) = int_kernel(ibasisfunc)%getvariance()

        if (parameters%int_over_background) then
          slave_result%model(:, ibasisfunc, ielement)         = int_model(ibasisfunc)%getintegral()
          call int_model(ibasisfunc)%freeme()
        end if
        if (parameters%int_over_hetero) then
          slave_result%hetero_model(:, ibasisfunc, ielement)  = int_hetero(ibasisfunc)%getintegral()
          call int_hetero(ibasisfunc)%freeme()              
        end if

        if (parameters%plot_wavefields) then
          ! Wavefield plotting is only allowed for volumetric mode, where nbasisfuncs_per_elem==1
          ! Therefore, this block is only called once.
          do ikernel = 1, parameters%nkernel
            norm_fields = real(niterations(ikernel) * nptperstep, kind=dp)
            slave_result%fw_field(:, :, ikernel, ielement)   &
              = real(fw_field_out(:, :, ikernel) / norm_fields, kind=sp)
            slave_result%bw_field(:, :, ikernel, ielement)   &
              = real(bw_field_out(:, :, ikernel) / norm_fields, kind=sp)
            slave_result%conv_field(:, 1, ikernel, ielement) &
              = real(conv_field_out(:, ikernel)  / norm_fields, kind=sp)
          end do ! ikernel
        end if
        call int_kernel(ibasisfunc)%freeme()
      end do ! ibasisfunc
    end if 

  end do loop_elements ! ielement

  call timeshift_fwd%freeme()
  call timeshift_bwd%freeme()
end function slave_work
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Do the integration/projection of the (1D) background model on the base functions
function integrate_1d_model(sem_data, inv_mesh, ielement) result(int_model)
  use global_parameters,                     only: sp, dp
  use readfields,                            only: semdata_type
  use inversion_mesh,                        only: inversion_mesh_data_type
  use background_model,                      only: backgroundmodel_type, nmodel_parameters
  use mc_integration,                        only: integrated_type, allallconverged, allisconverged

  type(semdata_type),   intent(in)              :: sem_data
  type(inversion_mesh_data_type), intent(in)    :: inv_mesh
  integer,              intent(in)              :: ielement
  type(integrated_type), allocatable            :: int_model(:)

  real(kind=dp),    allocatable                 :: random_points(:,:) 
  real(kind=dp),    allocatable                 :: model_random_points(:,:) 
  type(backgroundmodel_type)                    :: coeffs_random_points
  integer                                       :: ibasisfunc, istep_model
  integer                                       :: nbasisfuncs_per_elem
  integer, parameter                            :: nptperstep_model = 100 !< Points per MC iteration for 
                                                                          !! Integration of model parameters
                                                                          !! Larger than for kernels, since
                                                                          !! evaluation is very cheap
  integer, parameter                            :: max_iterations = 9999  !< The model integration should 
                                                                          !! not take forever

  nbasisfuncs_per_elem = inv_mesh%nbasisfuncs_per_elem
  allocate(int_model(nbasisfuncs_per_elem))
  allocate(random_points(3, nptperstep_model))
  allocate(model_random_points(nmodel_parameters, nptperstep_model))


  !  Calculate integrated model parameters for this element
  write(lu_out,'(A)', advance='no') ' Integrate model parameters...' 
  call flush(lu_out)
  istep_model = 0
  do ibasisfunc = 1, nbasisfuncs_per_elem
     call int_model(ibasisfunc)%initialize_montecarlo(nfuncs = nmodel_parameters,   & 
                                                      volume = 1d0,                 & 
                                                      allowed_error = 1d-2,         &
                                                      allowed_relerror = 1d-2)
  end do
  do while (.not.allallconverged(int_model))
     random_points = inv_mesh%generate_random_points( ielement, nptperstep_model, .true.)
                                                      
     do ibasisfunc = 1, nbasisfuncs_per_elem
        coeffs_random_points = sem_data%load_model_coeffs(random_points)
        model_random_points  = coeffs_random_points%weight(inv_mesh%weights(ielement,      &  
                                                                            ibasisfunc,    & 
                                                                            random_points) )

        call int_model(ibasisfunc)%check_montecarlo_integral(transpose(model_random_points))
     end do
     istep_model = istep_model + 1
     if (istep_model.ge.max_iterations) exit 
  end do
  write(lu_out, '(A,I5,A)') ' done after ', istep_model, ' steps.'
  call flush(lu_out)

end function integrate_1d_model
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Do the integration/projection of a (3D) heterogeneity model on the base functions
function integrate_3d_model(het_model, inv_mesh, ielement) result(int_hetero)
  use global_parameters,                     only: sp, dp
  use inversion_mesh,                        only: inversion_mesh_data_type
  use heterogeneities,                       only: hetero_type, nmodel_parameters_hetero         
  use mc_integration,                        only: integrated_type, allallconverged

  type(hetero_type),    intent(in)              :: het_model
  type(inversion_mesh_data_type), intent(in)    :: inv_mesh
  integer,              intent(in)              :: ielement
  type(integrated_type), allocatable            :: int_hetero(:)

  real(kind=dp),    allocatable                 :: random_points(:,:) 
  real(kind=dp),    allocatable                 :: hetero_random_points(:,:) 
  integer                                       :: ibasisfunc, istep_model
  integer                                       :: nbasisfuncs_per_elem
  integer, parameter                            :: nptperstep_model = 10  !< Points per MC iteration for 
                                                                          !! Integration of model parameters
                                                                          !! Larger than for kernels, since
                                                                          !! evaluation is very cheap

  nbasisfuncs_per_elem = inv_mesh%nbasisfuncs_per_elem
  allocate(int_hetero(nbasisfuncs_per_elem))
  allocate(random_points(3,nptperstep_model))
  allocate(hetero_random_points(nmodel_parameters_hetero,nptperstep_model))


  write(lu_out,'(A)', advance='no') ' Integrate heterogeneities... '
  call flush(lu_out)
  istep_model = 0
  do ibasisfunc = 1, nbasisfuncs_per_elem
     call int_hetero(ibasisfunc)%initialize_montecarlo(nfuncs = nmodel_parameters_hetero,   & 
                                                       volume = 1d0,                 & 
                                                       allowed_error = 1d-2,         &
                                                       allowed_relerror = 1d-1)
  end do
  do while (.not.allallconverged(int_hetero))
     random_points = inv_mesh%generate_random_points( ielement, nptperstep_model, .true.)
          
     do ibasisfunc = 1, nbasisfuncs_per_elem
        ! load weighted random points from heterogeneity structure
        hetero_random_points= het_model%load_model_coeffs(random_points,                  &
                                                          inv_mesh%weights(ielement,      &  
                                                                           ibasisfunc,    & 
                                                                           random_points) )
        call int_hetero(ibasisfunc)%check_montecarlo_integral(transpose(hetero_random_points))
     end do
     istep_model = istep_model + 1
     if (istep_model.ge.9999) exit ! The model integration should not take forever
  end do
  write(lu_out, '(A,I4,A)') ' done after ', istep_model, ' steps.'
  call flush(lu_out)

end function integrate_3d_model
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
logical function travel_time_larger_Tmax(source, receiver, coordinates, v_max, &
                                         planet_radius)
  use source_class,                 only : src_param_type
  use receiver_class,               only : rec_param_type
  type(src_param_type)                  :: source
  type(rec_param_type)                  :: receiver
  real(kind=dp), intent(in)             :: coordinates(3)
  real(kind=sp), intent(in)             :: v_max(2)
  real(kind=dp), intent(in)             :: planet_radius
  
  real(kind=dp)                         :: distance
  real(kind=dp)                         :: t_max_p, t_max_s

  distance = norm2(source%get_r(planet_radius) - coordinates) + &
             norm2(receiver%get_r(planet_radius) - coordinates)

  t_max_p = distance / v_max(1)
  t_max_s = distance / v_max(2)
  
  travel_time_larger_Tmax = (t_max_p > receiver%t_last_p .and. &
                             t_max_s > receiver%t_last_s)

end function travel_time_larger_Tmax
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
