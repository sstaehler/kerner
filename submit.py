#!/usr/bin/env python
# coding: utf-8
"""
Submit a MC kerner job locally or into a HPC queue
Author: Simon Stähler, LMU München
"""

import os
import argparse
import shutil
import glob
import datetime
import subprocess
import math
from numpy import isclose
from netCDF4 import Dataset


def nextpow2(i):
    """
    Find 2^n that is equal to or greater than.
    """
    n = 1
    while n < i:
        n *= 2
    return n


def estimate_memory():
    # Estimate memory usage which cannot be controlled by input parameters

    memory_mesh = 4 * (2*npoints_fwd + 8*nelems_fwd + 25*nelems_fwd) * 2
    # print 'Mesh size in memory: %f MB' % (memory_mesh/(2**20))

    nomega = nextpow2(ndumps_fwd)
    if full_strain:
        ndim = 6
    else:
        ndim = 1

    memory_fft = (8 * (nomega + ndumps_fwd) * ndim *
                  float(params['ELEMENTS_PER_TASK']))
    # print 'Memory for FFT types: %f MB' % (memory_fft/(2**20))

    memory_fields = (8 * float(params['ELEMENTS_PER_TASK']) * ndim *
                     (5*nomega + 2*ndumps_fwd))
    # print 'Memory for Wavefields: %f MB' % (memory_fields/(2**20))

    '''
    The following memory requirements were determined using the Massif heap
    profiler from the Valgrind package on two Ubuntu 14.04 machines. It might
    be that the memory requirements of the HDF5 library are different on other
    architectures. Memory usage of HDF5 library seems to be pretty constant
    around 120 MB.
    '''
    memory_hdf5 = 150 * 2**20
    # print 'Memory for HDF5/NetCDF4 library: %f MB' % (memory_hdf5/(2**20))

    # Memory usage of KD-Trees is roughly 80 Byte per mesh point
    memory_kdtree = 80 * npoints_fwd
    # print 'Memory for KD-Trees: %f MB' % (memory_kdtree/(2**20))

    memory_total = (memory_mesh + memory_fft + memory_hdf5 +
                    memory_kdtree + memory_fields)
    return memory_total


def auto_buffer_size(memory_available):
    if full_strain:
        ndim = 6
    else:
        ndim = 1

    memory_for_buffers = (memory_available - estimate_memory())*0.9

    if merged_db:
        size_one_strain_element = (8.0 *         # 8 Byte per number
                                   25 *          # Number of GLL pts per elem
                                   ndumps_fwd *  # number of time samples
                                   ndim *        # number of strain dimensions
                                   6)            # 6 files (4 fwd, 2 bwd)
        size_one_disp_element = (4.0 *         # 4 Byte per number
                                 25 *          # Number of GLL points per elem
                                 ndumps_fwd *  # number of time samples
                                 15)           # 15 disp. dimensions

    else:
        size_one_strain_element = (4.0 *  # 8 Byte per number
                                   25 *  # Number of GLL points per elem
                                   ndumps_fwd *  # number of time samples
                                   ndim *        # number of strain dimensions
                                   6)            # 6 files (4 fwd, 2 bwd)
        size_one_disp_element = (4.0 *         # 4 Byte per number
                                 3 *           # 3 dimensions
                                 ndumps_fwd *  # number of time samples
                                 6)            # 6 files (4 fwd, 2 bwd)

    # Rule: Strain buffer gets 60% of the available memory, displ. buffer 40%
    size_strain_buffer = int(memory_for_buffers * 0.6 /
                             size_one_strain_element)
    size_disp_buffer = int(memory_for_buffers * 0.4 /
                           size_one_disp_element)

    # memory_buffers_strain = size_one_strain_element * size_strain_buffer
    # print 'Strain buffer size: %f MB' % (memory_buffers_strain/(2**20))
    # memory_buffers_disp = size_one_disp_element *  size_disp_buffer
    # print 'Displ. buffer size: %f MB' % (memory_buffers_disp/(2**20))

    if memory_for_buffers < 0:
        raise ValueError('Not enough memory for buffers')

    return str(size_strain_buffer), str(size_disp_buffer)


def read_receiver_dat(rec_file):
    with open(rec_file) as f_rec:
        nkernel_total = 0
        fullstrain_kernel = False
        # Read number of receivers
        str_line = f_rec.readline()
        nrec = int(str_line.split()[0])
        print('Number of receivers: %d' % nrec)
        # Read seismogram component
        _ = f_rec.readline()
        # seis_cmp = str_line.split()[0]
        for irec in range(0, nrec):
            str_line = f_rec.readline()
            # rec_name = str_line.split()[0]
            # rec_lat = float(str_line.split()[1])
            # rec_lon = float(str_line.split()[2])
            nkernel = int(str_line.split()[4])
            # print 'Receiver: %s, coordinates: (%f, %f), %d kernels' % \
            #          (rec_name, rec_lat, rec_lon, nkernel)
            for ikernel in range(0, nkernel):
                str_line = f_rec.readline()
                kernel_name = str_line.split()[0]
                # filter_name = str_line.split()[1]
                # misfit_name = str_line.split()[2]
                time_window_start = float(str_line.split()[3])
                time_window_stop = float(str_line.split()[4])
                # Check whether kernel time windows run into taper at the end of
                # the wavefield time series.
                if time_window_stop > dt_fwd * ndumps_fwd * 0.95:
                    errmsg = 'Time window for kernel %s, (%7.1fs, %7.1fs) ' + \
                             'exceeds safe length (0.95 * sim. len. = %7.1f s)'
                    raise RuntimeError(errmsg %
                                       (kernel_name,
                                        time_window_start, time_window_stop,
                                        0.95 * dt_fwd * ndumps_fwd))

                model_param = str_line.split()[5]
                if model_param in ('vs', 'rho', 'vsh', 'vsv',
                                   'eta', 'phi', 'xi', 'mu'):
                    fullstrain_kernel = True
                elif model_param not in ('lam', 'vp', 'vph', 'vpv'):
                    raise RuntimeError('Unknown model parameter %s in %s' %
                                       (model_param, rec_file))

            nkernel_total += nkernel

        print('Number of kernels: %d' % nkernel_total)

    return nkernel_total, nrec, fullstrain_kernel


def define_arguments():
    helptext = 'Create MC Kernel input file and submit job.'
    formatter_class = argparse.RawTextHelpFormatter
    parser = argparse.ArgumentParser(description=helptext,
                                     formatter_class=formatter_class)

    helptext = "Job directory name. \n" + \
               "If this argument is an absolute path (starting with /),\n" + \
               "this path will be created and used. \n" + \
               "If not, a run directory with this name will be created\n" + \
               "in the 'RUNS_DIRECTORY' directory set in\n" + \
               "make_mc_kernel.macros."
    parser.add_argument('job_name', help=helptext)

    helptext = "Input file to use. It will overwrite default values, \n" + \
               "but will be overwritten by any argument to this function."
    parser.add_argument('-i', '--input_file', help=helptext)

    helptext = "Number of slaves to use. Default is local number \n" + \
               "of CPUs - 1 (for the master).                     \n" + \
               "Number of slaves to use. If --queue==SuperMUC, it \n" + \
               "will be rounded up to a multiple of 16 (thin island) \n" + \
               "or 40 (fat island)"
    parser.add_argument('-n', '--nslaves', type=int,
                        default=ncpu - 1,
                        metavar='N',
                        help=helptext)

    helptext = "Description of run, which is saved in \n" + \
               "JOB_NAME/README.run. \n"
    parser.add_argument('-m', '--message', metavar='JOB_DESCRIPTION_MESSAGE',
                        help=helptext)

    helptext = "Plot wavefields and waveform kernels in addition to the\n" + \
               "normal misfit kernels. Be aware that this feature\n" + \
               "requires a large amount of memory."
    parser.add_argument('--plot_wavefields', default=False, action='store_true',
                        help=helptext)

    helptext = "Amount of memory available per task in MB.\n" + \
               "This number is used to determine the size of the read\n" + \
               "buffers. If the required memory is larger than this\n" + \
               "number, an error is thrown."
    parser.add_argument('-a', '--available_memory', type=int,
                        help=helptext)

    helptext = "Queue to use. Default is local, which starts a job\n" + \
               "with MPIRUN"
    parser.add_argument('-q', '--queue',
                        choices=['SuperMUC', 'background', 'foreground',
                                 'monch', 'torque'],
                        default='background',
                        help=helptext)

    ############################################################################
    # Queue options
    ############################################################################
    hpc_queue = parser.add_argument_group('Options specific to HPC queues')
    hpc_queue.add_argument('--wall_time', type=int, metavar='WALLTIME_IN_H',
                           default=10,
                           help='Walltime in hours')
    hpc_queue.add_argument('--mail_address',
                           help='Mail address for HPC notifications')
    hpc_queue.add_argument('--job_class', choices=['fat', 'thin'],
                           default='thin',
                           help='Job class on SuperMUC')
    hpc_queue.add_argument('--tasks_per_node', type=int,
                           help='Tasks per node on SuperMUC')
    hpc_queue.add_argument('--parallel_reading', default=False,
                           action='store_true',
                           help='Use parallel NetCDF4 for reading.')

    ############################################################################
    # AxiSEM wavefield directories
    ############################################################################

    axisem_dirs = parser.add_argument_group('AxiSEM run directories')
    axisem_dirs.add_argument('--fwd_dir', default='./wavefield/fwd/',
                             help='Path to AxiSEM forward run')
    axisem_dirs.add_argument('--bwd_dir', default='./wavefield/bwd/',
                             help='Path to AxiSEM backward run')

    ############################################################################
    # input files
    ############################################################################

    input_files = parser.add_argument_group('Required input files')
    input_files.add_argument('--src_file', default='CMTSOLUTION',
                             help='Path to source file in CMTSOLUTION format')
    input_files.add_argument('--rec_file', default='receiver.dat',
                             help='Path to receiver and kernel file')
    input_files.add_argument('--filt_file', default='filters.dat',
                             help='Path to filter file')
    input_files.add_argument('--stf_file', default='stf_20s.dat',
                             help='Path to Source Time Function file')

    ############################################################################
    # Mesh file-related options
    ############################################################################

    mesh_files = parser.add_argument_group('Inversion mesh')
    helptext = "Select the mesh file type. Allowed values are\n" + \
               "abaqus     : .inp file, can be generated with\n" + \
               "             Qubit and other codes .Can contain\n" + \
               "             various geometries and multiple\n" + \
               "             sub-objects.\n" + \
               "             Supported geometries (so far):\n" + \
               "             tetrahedra, triangles.\n" + \
               "             Set file name in MESH_FILE_ABAQUS.\n" + \
               "tetrahedral: tetrahedral mesh in two separate files\n" + \
               "             1. coordinates of the vertices\n" + \
               "                (MESH_FILE_VERTICES)\n" + \
               "             2. the connectivity of the facets of\n" + \
               "                tetrahedrons (MESH_FILE_FACETS)"
    mesh_files.add_argument('--mesh_file_type', default='abaqus',
                            choices=['abaqus', 'tetrahedral'],
                            help=helptext)
    mesh_files.add_argument('--mesh_file_abaqus',
                            default='Meshes/flat_triangles.inp',
                            help='Path to Abaqus mesh file')
    helptext = "Path to Vertices file:\n" + \
               "(only if --mesh_file_type=tetrahedral)"
    mesh_files.add_argument('--mesh_file_vertices',
                            default='tests/vertices.TEST',
                            help=helptext)
    helptext = "Path to Facets file:\n" + \
               "(only if --mesh_file_type=tetrahedral)"
    mesh_files.add_argument('--mesh_file_facets',
                            default='tests/facets.TEST',
                            help=helptext)

    ############################################################################
    # Kernel-related options
    ############################################################################

    kernel_options = parser.add_argument_group('Kernel calculation options')
    helptext = "Calculate kernels for absolute values (e.g. Vp) instead\n" + \
               "of relative perturbations (dVp) with respect to the \n" + \
               "background model"

    kernel_options.add_argument('--kernel_for_absolute_perturbations',
                                action="store_true", default=False,
                                help=helptext)

    helptext = "On which base functions are the kernels defined?\n" + \
               "volumetric (default): Each voxel is a base function\n" + \
               "                      (Boschi & Auer)\n" + \
               "onvertices:           Each vertex has a set of non-\n" + \
               "                      orthogonal base functions defined\n" + \
               "                      on it (Nolet & Sigloch)"
    kernel_options.add_argument('--int_type',
                                choices=['volumetric', 'onvertices'],
                                default='volumetric', help=helptext)

    helptext = "For plotting reasons one may wish to skip the \n" + \
               "integration over cell-volume. This makes the \n" + \
               "values independent of cell size \n" + \
               "Resulting kernels bear the unit [s/m^3]"
    kernel_options.add_argument('--no_int_over_volume', action="store_true",
                                default=False, help=helptext)

    ############################################################################
    # Monte Carlo-related options
    ############################################################################

    mc_options = parser.add_argument_group('Monte Carlo options')
    helptext = "Number of points on which the kernel is evaluated per\n" + \
               "Monte Carlo step. Default value is 4."
    mc_options.add_argument('--points_per_mc_step', type=int, default=4,
                            help=helptext)

    helptext = "Maximum number of Monte Carlo iterations. Allows to \n" + \
               "skip evaluation of single problematic cells. \n" + \
               "Default value is 1E6"
    mc_options.add_argument('--maximum_iterations', type=int, default=1000000,
                            help=helptext)

    helptext = "Allowed absolute error before Monte Carlo integration \n" + \
               "is considered to be converged. When calculating this \n" + \
               "value, the volume is not considered, no matter whether \n" + \
               "--no_int_over_volume is set or not."
    mc_options.add_argument('--allowed_error', type=float, default=1e-4,
                            help=helptext)

    helptext = "Allowed relative error before Monte Carlo integration \n" + \
               "in one cell is considered to be converged. \n" + \
               "Default value is 1e-2"
    mc_options.add_argument('--allowed_relative_error', type=float,
                            default=1e-2, help=helptext)

    helptext = "Use pseudorandom numbers instead of quasirandom"
    mc_options.add_argument('--use_pseudorandom_numbers', action="store_true",
                            default=False, help=helptext)

    ############################################################################
    # Debugging-related options
    ############################################################################

    debug_options = parser.add_argument_group('Debugging options')
    helptext = "This activates the optional (linearity test) to \n" + \
               "integrate relative kernels over model perturbations in \n" + \
               "percent, to assess how well our kernels predict measured \n" + \
               "traveltime perturbations for the same model. This only \n" + \
               "makes sense when not calculating absolute kernels. "
    debug_options.add_argument('--int_over_3d_heterogeneities',
                               action="store_true", default=False,
                               help=helptext)

    helptext = "Path to heterogeneity file"
    debug_options.add_argument('--het_file', default='tests/savani.rtpv',
                               help=helptext)

    helptext = "Integrate the kernel over the background model. \n" + \
               "Classically, this was assumed to result in the travel \n" + \
               "time of a phase. This assumption is highly dubious for\n" + \
               "wavefield-derived kernels. For legacy reasons, we can\n" + \
               "still leave it in.  Adds a version of the background\n" + \
               "model interpolated on the inversion mesh to the \n" + \
               "output file."
    debug_options.add_argument('--int_over_background_model',
                               action="store_true", default=False,
                               help=helptext)

    helptext = "Every slave writes out the values of all the kernels \n" + \
               "and their respective estimated errors into his \n" + \
               "OUTPUT_???? file after each MC step. This can lead to \n" + \
               "huge ASCII files (>1GB) with inane line lengths \n" + \
               "(approx. 20 x nkernel).  However, it might be \n" + \
               "interesting to study the convergence behaviour."
    debug_options.add_argument('--write_detailed_convergence',
                               action="store_true", default=False,
                               help=helptext)

    helptext = "Do not deconvolve the Source Time Function and \n" + \
               "reconvolve with the one set in --stf_file, but just \n" + \
               "timeshift the wavefields."
    debug_options.add_argument('--no_deconvolve_stf', action="store_true",
                               default=False, help=helptext)

    helptext = "Integration scheme to calculate scalar kernels from \n" + \
               "seismograms and waveforms.  \n" + \
               "parseval (default): Integration in frequency domain, \n" + \
               "                    using the Parseval theorem.  \n" + \
               "trapezoidal:        Integration in time domain using \n" + \
               "                    the trapezoidal rule."
    debug_options.add_argument('--integration_scheme',
                               choices=['parseval', 'trapezoidal'],
                               default='parseval', help=helptext)

    ############################################################################
    # Output-related options
    ############################################################################

    output_options = parser.add_argument_group('Output options')
    helptext = "Output format when dumping kernels and wavefields. \n" + \
               "Choose between xdmf, Yale-style csr binary format \n" + \
               "(compressed sparse row) and ascii. Yet, the allowed \n" + \
               "error below is assumed as the truncation threshold in \n" + \
               "csr and ascii storage"
    output_options.add_argument('--dump_type', choices=['xdmf', 'ascii', 'csr'],
                                default='xdmf', help=helptext)

    helptext = "Write out Seismograms (raw full trace, filtered full \n" + \
               "trace and cut trace) into run_dir/SEISMOGRAMS. Produces \n" + \
               "three files per kernel.  Disable to avoid congesting \n" + \
               "your file system."
    output_options.add_argument('--write_seismograms', default=True,
                                help=helptext)

    helptext = "Prefix of output file names.  \n" + \
               "Kernel files are called $OUTPUT_FILE_kernel.xdmf \n" + \
               "Wavefield movies are called $OUTPUT_FILE_wavefield.xdmf"
    output_options.add_argument('--out_prefix', default='kerner',
                                help=helptext)

    ############################################################################
    # Performance-related options
    ############################################################################

    performance_options = parser.add_argument_group(
        'Performance-related options')
    helptext = "Size of buffer for strain values. Since the strain \n" + \
               "has to be calculated from the displacement stored in \n" + \
               "the AxiSEM files, increasing this buffer size saves \n" + \
               "IO access and CPU time."
    performance_options.add_argument('--strain_buffer_size', type=int,
                                     default=1000, help=helptext)

    helptext = "Size of buffer for displacement values. Displacement \n" + \
               "values are use to calculate strain later. Having a \n" + \
               "separate buffer here allows to save some IO accesses."
    performance_options.add_argument('--displ_buffer_size', type=int,
                                     default=100, help=helptext)

    helptext = "Number of elements per Slave task. \n" + \
               "A larger value allows to the Slave to have more \n" + \
               "contiguous parts of the earth to work on, smaller \n" + \
               "values improve load balancing. It should be chosen such \n" + \
               "that each slave gets at least 50-100 tasks to work on."
    performance_options.add_argument('--elements_per_task', type=int,
                                     default=100, help=helptext)

    helptext = "Do not sort the mesh elements. Just for debugging."
    performance_options.add_argument('--no_sort_mesh_elements',
                                     action="store_true", default=False,
                                     help=helptext)

    helptext = "Create a file with intermediate results. Probably \n" + \
               "useful, if you have reason to expect the job to be \n" + \
               "cancelled.  Can inhibit performance significantly for \n" + \
               "large numbers of kernels and large inversion grids"
    performance_options.add_argument('--create_intermediate', default=False,
                                     help=helptext)

    helptext = "Start creation of intermediate files only after N hours\n" + \
               "Set this to value close to walltime on HPC systems"
    performance_options.add_argument('--intermediate_dump_time', type=float,
                                     default=0, help=helptext)

    helptext = "Mask the source and the receiver element and set the \n" + \
               "kernel to zero in each. A rough way to avoid spending \n" + \
               "hours until convergence in these two elements in reached."
    performance_options.add_argument('--mask_source_receiver',
                                     action="store_true", default=False,
                                     help=helptext)

    helptext = "Dampen the kernel in a radius around source and receiver.\n" + \
               "If a negative value is chosen, damping is \n" + \
               "switched off (default)."
    performance_options.add_argument('--damp_radius_source_receiver',
                                     type=float, default=-100.E3,
                                     help=helptext)

    helptext = "FFTW Planning to use\n" + \
               " Options:\n" + \
               " ESTIMATE:   Use heuristic to find best FFT plan\n" + \
               " MEASURE:    Compute several test FFTs to find best plan \n" + \
               "             (default)\n" + \
               " PATIENT:    Compute a lot of test FFTs to find best plan\n" + \
               " EXHAUSTIVE: Compute an awful amount of test FFTs to \n" + \
               "             find best plan\n" + \
               " This option did not prove to be very useful on most systems."
    performance_options.add_argument('--fftw_plan', default='MEASURE',
                                     choices=['ESTIMATE', 'MEASURE',
                                              'PATIENT', 'EXHAUSTIVE'],
                                     help=helptext)
    return parser

# Try to get the number of CPUS. If psutils is not insalled, use a default
# value of 2.
try:
    import psutil
    ncpu = psutil.cpu_count()
except ImportError:
    ncpu = 2

parser = define_arguments()

args = parser.parse_args()

# Parse input file, if one was given
if args.input_file:
    with open(args.input_file) as f:
        args_input_file = {}
        for line in f:
            # Skip comment and empty lines
            if line[0] != '#' and line.strip() != '':
                (key, val) = line.split()
                args_input_file[key] = val

# Merge input variables from input_file, arguments and default values
params = {}
# Loop over all possible arguments
for key, value in vars(args).items():
    if key not in ('nslaves', 'job_name', 'queue', 'available_memory'):
        # If an input file is selected, get values from there by default
        if args.input_file:
            if value == parser.get_default(key):
                key_out = key.upper()
                try:
                    value_out = str(args_input_file[key.upper()]).strip("'")
                except KeyError:
                    # If value is not set in input file
                    value_out = str(value)
            else:
                # Unless values were explicitly given
                key_out = key.upper()
                value_out = str(value)
        else:
            # In all other cases, take default values
            key_out = key.upper()
            value_out = str(value)

        params[key_out] = value_out

# Check for AxiSEM wavefield files to get mesh size
fwd_path = os.path.join(os.path.realpath(params['FWD_DIR']),
                        'MZZ', 'Data', 'ordered_output.nc4')
bwd_path = os.path.join(os.path.realpath(params['BWD_DIR']),
                        'PZ', 'Data', 'ordered_output.nc4')

fwd_path_merged = os.path.join(os.path.realpath(params['FWD_DIR']),
                               'merged_instaseis_db.nc4')
bwd_path_merged = os.path.join(os.path.realpath(params['BWD_DIR']),
                               'merged_instaseis_db.nc4')

if os.path.exists(fwd_path):
    nc_fwd = Dataset(fwd_path)
    merged_db = False
elif os.path.exists(fwd_path_merged):
    nc_fwd = Dataset(fwd_path_merged)
    merged_db = True
else:
    errmsg = 'Could not find a wavefield file in the fwd_dir %s\n%s' % \
             (params['FWD_DIR'], fwd_path)
    raise IOError(errmsg)

npoints_fwd = getattr(nc_fwd, "npoints")
nelems_fwd = getattr(nc_fwd, "nelem_kwf_global")
ndumps_fwd = getattr(nc_fwd, "number of strain dumps")
dt_fwd = getattr(nc_fwd, "strain dump sampling rate in sec")
period_fwd = getattr(nc_fwd, "dominant source period")
nc_fwd.close()

if os.path.exists(bwd_path):
    nc_bwd = Dataset(bwd_path)
elif os.path.exists(bwd_path_merged):
    nc_bwd = Dataset(bwd_path_merged)
else:
    errmsg = 'Could not find a wavefield file in the bwd_dir %s' % \
             params['BWD_DIR']
    raise IOError(errmsg)

npoints_bwd = getattr(nc_bwd, "npoints")
nelems_bwd = getattr(nc_bwd, "nelem_kwf_global")
ndumps_bwd = getattr(nc_bwd, "number of strain dumps")
dt_bwd = getattr(nc_bwd, "strain dump sampling rate in sec")
period_bwd = getattr(nc_bwd, "dominant source period")
nc_bwd.close()

# Sanity check, whether fwd and bwd mesh have the same sizes and the same number
# of wavefield time steps.
if not isclose(dt_fwd, dt_bwd):
    errmsg = 'The two mesh files do not have the same dt:' + \
        'Forward wavefield:  %9.6f sec\n' % dt_fwd + \
        'Backward wavefield: %9.6f sec\n' % dt_bwd
    raise ValueError(errmsg)

if not isclose(period_fwd, period_bwd):
    errmsg = 'The two mesh files do not have the same period:' + \
        'Forward wavefield:  %5.2f sec\n' % period_fwd + \
        'Backward wavefield: %5.2f sec\n' % period_bwd
    raise ValueError(errmsg)

if npoints_fwd != npoints_bwd or \
   nelems_fwd != nelems_bwd or \
   ndumps_fwd != ndumps_bwd:

    raise RuntimeError('Forward and backward run did not use' +
                       'the same parameters')

# Read receiver file and get number of receivers, kernels and whether the
# full strain has to be read for any kernel (increases the memory footprint
# of the buffers by factor of 6)
nrec, nkernel, full_strain = read_receiver_dat(params['REC_FILE'])

# Get mpirun and runs_directory from make_mckernel.macros
with open('make_mc_kernel.macros') as f:
    for line in f:
        if line.strip() != '':
            key = line.split()[0]
            if key == 'MPIRUN':
                mpirun_cmd = line.split()[2]
            elif key == 'RUNS_DIRECTORY':
                runs_directory = line.split()[2]

# Create run_dir
# Check whether absolute or relative path is given
if os.path.isabs(args.job_name):
    run_dir = args.job_name
else:
    run_dir = os.path.join(runs_directory, args.job_name)

if os.path.exists(run_dir):
    raise RuntimeError('Run directory \n %s \n already exists' % run_dir)

os.mkdir(run_dir)


# Define buffer sizes based on available memory
if args.available_memory:
    params['STRAIN_BUFFER_SIZE'], params['DISPL_BUFFER_SIZE'] = \
        auto_buffer_size(args.available_memory*2**20)


params_out = {}

# Copy necessary files to rundir
for key, value in params.items():

    if key == 'SRC_FILE':
        src_file_name = os.path.split(value)[1]
        shutil.copy(value, os.path.join(run_dir, src_file_name))
        params_out[key] = src_file_name

    elif key == 'REC_FILE':
        rec_file_name = os.path.split(value)[1]
        shutil.copy(value, os.path.join(run_dir, rec_file_name))
        params_out[key] = rec_file_name

    elif key == 'FILT_FILE':
        filt_file_name = os.path.split(value)[1]
        shutil.copy(value, os.path.join(run_dir, filt_file_name))
        params_out[key] = filt_file_name

    elif key == 'STF_FILE':
        stf_file_name = os.path.split(value)[1]
        shutil.copy(value, os.path.join(run_dir, stf_file_name))
        params_out[key] = stf_file_name

    elif key == 'MESH_FILE_TYPE':
        params_out[key] = value
        if value == 'abaqus':
            # get file name only
            mesh_file_name = os.path.split(params['MESH_FILE_ABAQUS'])[1]
            shutil.copy(params['MESH_FILE_ABAQUS'],
                        os.path.join(run_dir, mesh_file_name))
            params_out['MESH_FILE_ABAQUS'] = mesh_file_name
        elif value == 'tetrahedral':
            mesh_vertices_name = os.path.split(params['MESH_FILE_VERTICES'])[1]
            mesh_facets_name = os.path.split(params['MESH_FILE_FACETS'])[1]
            shutil.copy(params['MESH_FILE_VERTICES'],
                        os.path.join(run_dir, mesh_vertices_name))
            shutil.copy(params['MESH_FILE_FACETS'],
                        os.path.join(run_dir, mesh_facets_name))
            params_out['MESH_FILE_VERTICES'] = mesh_vertices_name
            params_out['MESH_FILE_FACETS'] = mesh_facets_name
        else:
            raise KeyError('Unknown mesh file type: %s' % args.mesh_file_type)

    elif key == 'INT_OVER_3D_HETEROGENEITIES':
        params_out[key] = value
        het_file_name = os.path.split(params['HET_FILE'])[1]
        shutil.copy(params['HET_FILE'],
                    os.path.join(run_dir, het_file_name))
        params_out['HET_FILE'] = het_file_name

    elif key in ('FWD_DIR', 'BWD_DIR'):
        # Set mesh dir to absolute path
        params_out[key] = os.path.realpath(value)

    elif key not in ('NSLAVES', 'JOB_NAME', 'MESH_FILE_ABAQUS', 'WALL_TIME',
                     'MESH_FILE_VERTICES', 'MESH_FILE_FACETS', 'HET_FILE',
                     'INPUT_FILE', 'MESSAGE', 'AVAILABLE_MEMORY',
                     'MAIL_ADDRESS', 'JOB_CLASS', 'TASKS_PER_NODE'):
        # These variables should not be written into the MC Kernel input
        # file and are only used by the submit script.
        params_out[key] = value

# Open editor window to write run descriptor
out_readme = 'readme_temp.txt'
f_readme = open(out_readme, 'w')
current_time = str(datetime.datetime.now())
f_readme.write('MC KERNEL run for %d CPUs, started on %s\n' % (args.nslaves,
                                                               current_time))
f_readme.write('  by user ''%s'' on ''%s''\n' % (os.environ.get('USER'),
                                                 os.environ.get('HOSTNAME')))
f_readme.write('  Minimum AxiSEM period: %4.1f\n' % period_fwd)
if args.message:
    f_readme.write(args.message)
f_readme.close()

# Move README file to rundir
shutil.move(out_readme, os.path.join(run_dir, 'README.run'))

# Create directory for seismogram output
os.mkdir(os.path.join(run_dir, 'Seismograms'))

# Create directory for filter output
os.mkdir(os.path.join(run_dir, 'Filters'))

# Create input file for run
out_input_file = os.path.join(run_dir, 'inparam')
with open(out_input_file, 'w') as f_out:
    for key, value in params_out.items():
        if value.find('/') == -1:
            f_out.write('%s  %s\n' % (key, value))
        else:
            f_out.write('%s "%s"\n' % (key, value))

# Make MC kernel code
subprocess.check_call('make -sj', shell=True)

# Copy code files into run_dir, tar it and delete it.
# A bit clumsy, but ensures that the internal path is Code/*.f90 etc.
code_dir = os.path.join(run_dir, 'Code')
archive_name = os.path.join(run_dir, 'Code')
os.mkdir(code_dir)
for f90_file in glob.glob('./src/*.90'):
    shutil.copy(f90_file, code_dir)
shutil.copy('Makefile', code_dir)
shutil.copy('make_mc_kernel.macros', code_dir)
shutil.make_archive(archive_name, 'gztar', code_dir)
shutil.rmtree(code_dir)

shutil.copy('./bin/mc_kernel', run_dir)

if args.queue == 'background':
    # Change dir and submit
    os.chdir(run_dir)

    cmd_string = \
        'nohup %s -n %d -quiet ./mc_kernel inparam 2>&1 &> OUTPUT_0000 &'
    run_cmd = cmd_string % (mpirun_cmd, args.nslaves + 1)
    print('Starting local job in %s' % run_dir)
    print('Check %s/OUTPUT_0000 for progress' % run_dir)
    subprocess.call(run_cmd, shell=True)

if args.queue == 'foreground':
    # Change dir and submit
    os.chdir(run_dir)

    cmd_string = \
        '%s -n %d ./mc_kernel inparam'
    run_cmd = cmd_string % (mpirun_cmd, args.nslaves + 1)
    print('Starting local job in %s' % run_dir)
    subprocess.check_call(run_cmd, shell=True)
    
if args.queue == 'torque':
    # Create a .pbs submission script for archer
    job_script = os.path.join(run_dir, 'submit_job.pbs')
    if not args.tasks_per_node:
        tasks_per_node = 24
    else:
        tasks_per_node = args.tasks_per_node
    nodes = math.ceil((args.nslaves)/tasks_per_node)
    with open(job_script, 'w') as f:
        text_out = "#! /bin/bash --login \n"
        text_out += "#PBS -l select=%s \n" % str(int(nodes))
        text_out += "#PBS -l walltime=%s:00:00 \n" % str(args.wall_time)
        text_out += "#PBS -A n03-ox2 \n"
        text_out += "cd %s \n" % run_dir
        text_out += "module unload PrgEnv-cray/5.2.82 \n"
        text_out += "module load PrgEnv-gnu \n"
        text_out += "module load fftw \n"
        text_out += "module load cray-netcdf-hdf5parallel/4.4.1.1 \n"
        text_out += "module load atp \n"
        text_out += "module load doxygen/1.8.4 \n"
        text_out += "export ATP_ENABLED=1 \n"
        text_out += "ulimit -s unlimited \n"
        text_out += "aprun -n %s %s/mc_kernel inparam \n" % (str(args.nslaves), run_dir)
        f.write(text_out)
    os.chdir(run_dir)
    subprocess.call(['qsub', job_script])    
    
elif args.queue == 'SuperMUC':
    # Create a LoadLeveler job script for SuperMUC

    # Master gets his own node, since his memory requirements can become quite
    # huge for big meshes and a large number of kernels and cannot be changed.
    job_script = os.path.join(run_dir, 'job.cmd')
    if args.job_class == 'fat':
        if not args.tasks_per_node:
            tasks_per_node = 40
        else:
            tasks_per_node = args.tasks_per_node

        nodes = math.ceil((args.nslaves)/tasks_per_node) + 1

        job_class = 'big'

        if args.available_memory > 6000.*(40./tasks_per_node):
            raise IOError('Fat island has only 6GB RAM per node')

    elif args.job_class == 'thin':
        if not args.tasks_per_node:
            tasks_per_node = 28
        else:
            tasks_per_node = args.tasks_per_node

        nodes = math.ceil((args.nslaves)/tasks_per_node) + 1

        if nodes > 20:
            job_class = 'general'
        else:
            job_class = 'micro'

        if args.available_memory > 2000.*(28./tasks_per_node):
            raise IOError('Thin island has only 2.0GB RAM per node')

    with open(job_script, 'w') as f:
        text_out = "# Job file automatically created by submit.py on %s\n" % \
                     str(datetime.datetime.now())
        text_out += "#@ output = job_$(jobid).out \n"
        text_out += "#@ error = job_$(jobid).err\n"
        text_out += "#@ job_type = parallel \n"
        text_out += "#@ network.MPI = sn_all,not_shared,us \n"
        text_out += "#@ notification=always \n"
        text_out += "#@ notify_user = staehler@geophysik.uni-muenchen.de \n"
        text_out += "#@ energy_policy_tag = MCKernel\n"
        text_out += "#@ minimize_time_to_solution = yes     \n"
        text_out += "#@ class = %s\n" % job_class
        text_out += "#@ tasks_per_node = %d\n" % tasks_per_node
        text_out += "#@ first_node_tasks=1\n"
        text_out += "#@ node = %d\n" % nodes
        text_out += "#@ wall_clock_limit = %d:00:00\n" % args.wall_time
        text_out += "#@ job_name = %s\n" % args.job_name
        text_out += "#@ initialdir = %s\n" % os.path.realpath(run_dir)
        text_out += "#@ queue \n"
        text_out += ". /etc/profile \n"
        text_out += ". /etc/profile.d/modules.sh \n"
        text_out += "module load netcdf/mpi \n"
        text_out += "module load fftw \n"
        text_out += "export MP_TIMEOUT=3600\n"
        text_out += "export MP_PULSE=0\n"
        text_out += "poe ./mc_kernel inparam 2>&1 &> OUTPUT_0000\n"
        f.write(text_out)
    print('Submitting to SuperMUC loadleveler queue')
    subprocess.call(['llsubmit', job_script])

elif args.queue == 'monch':
    with open(os.path.join(run_dir, 'sbatch.sh'), 'w') as f:
        text_out = "#!/bin/bash -l\n"
        text_out += "#SBATCH --ntasks=%d\n" % (args.nslaves + 1)
        text_out += "#SBATCH --ntasks-per-node=%d\n" % int(args.nslaves/20.)
        text_out += "#SBATCH --nodes=%d\n" % (int((args.nslaves+1) /
                                              int(args.nslaves/20.)))
        text_out += "#SBATCH --mem-per-cpu=%d\n" % int(args.available_memory)
        text_out += "#SBATCH --time=%d:00:00\n" % args.wall_time
        text_out += "#SBATCH --partition=_compute\n"
        text_out += "#SBATCH --job-name=%s\n" % args.job_name
        text_out += "#SBATCH --output=mc_kernel_out.o\n"
        text_out += "#SBATCH --error=mc_kernel_err.o\n"
        text_out += "echo The current job ID is $SLURM_JOB_ID\n"
        text_out += "echo Running on $SLURM_JOB_NUM_NODES nodes\n"
        text_out += "echo Using $SLURM_NTASKS_PER_NODE tasks per node\n"
        text_out += "echo A total of $SLURM_NTASKS tasks is used\n"
        text_out += "mpirun -n %d ./mc_kernel inparam_basic &> OUTPUT_0000\n" %\
                    (args.nslaves + 1)
        f.write(text_out)
    os.chdir(run_dir)
    run_cmd = 'sbatch sbatch.sh'
    print(run_cmd)
    subprocess.call(run_cmd, shell=True)
