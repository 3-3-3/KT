import time
import pandas as pd
import numpy as np
import argparse
import json


from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic import utils
from cosmic.evolve import Evolve

import schwimmbad
from schwimmbad import MultiPool, MPIPool

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_commandline():
    """Parse the arguments given on the command-line.
    """
    # Parse any inifile specification
    # We make this parser with add_help=False so that
    # it doesn't parse -h and print help.
    conf_parser = argparse.ArgumentParser(
        description=__doc__, # printed with -h/--help
        # Don't mess with format of description
        formatter_class=argparse.RawDescriptionHelpFormatter,
        # Turn off help, so we print all options in response to -h
        add_help=False
        )
    conf_parser.add_argument("--inifile",
                            help="Name of ini file of params",
                            metavar="FILE",)
    args, remaining_argv = conf_parser.parse_known_args()

    defaults = {}
    if not (args.inifile is None and (('-h' in remaining_argv) or ('--help' in remaining_argv))):
        BSEDict, seed_int, filters, convergence, sampling = utils.parse_inifile(args.inifile)
        defaults.update(sampling)
        defaults.update(filters)
        defaults.update(convergence)
        defaults.update({'seed' : seed_int})
        defaults.update({'inifile' : args.inifile})

    # Parse rest of arguments
    # Don't suppress add_help here so it will handle -h
    parser = argparse.ArgumentParser(
        # Inherit options from config_parser
        parents=[conf_parser]
        )
    parser.set_defaults(**defaults)
    parser.add_argument("--final-kstar1",
                        help="Specify the final condition of kstar1 "
                        ", you want systems to end at for your samples",
                        required=True, type=int, nargs='+')
    parser.add_argument("--final-kstar2",
                        help="Specify the final condition of kstar2, you want "
                        "systems to end at for your samples",
                        required=True, type=int, nargs='+')
    parser.add_argument("--Niter",
                        help="Number of iterations of binaries "
                        "to try, will check ever Nstep for convergence",
                        type=int, default=10000000)
    parser.add_argument("--Nstep",
                        help="Number of binaries to try before checking for "
                        "convergence, it will check ever Nstep binaries until "
                        "it reach Niter binaries", type=int, default=10000)
    parser.add_argument("--binary_state", nargs='+', type=int)
    parser.add_argument("--sampling_method")
    parser.add_argument("--primary_model", help="Chooses the initial primary mass function from: salpeter55, kroupa93, kroupa01", type=str)
    parser.add_argument("--binfrac_model", help="Chooses the binary fraction model from: a float between [0,1] and vanHaaften", type=float)
    parser.add_argument("--ecc_model", help="Chooses the initial eccentricity distribution model from: thermal, uniform, and sana12", type=str)
    parser.add_argument("--porb_model", help="Chooses the initial orbital period distribution model from: log_uniform and sana12", type=str)
    parser.add_argument("--SF_start", help="Sets the time in the past when star formation initiates in Myr", type=float)
    parser.add_argument("--SF_duration", help="Sets the duration of constant star formation in Myr", type=float)
    parser.add_argument("--metallicity", type=float)
    parser.add_argument("--total_mass", help="Gives the total mass to sample in creating population", type=float)
    parser.add_argument("--convergence_params", nargs='+', help="specifies the list of parameters for which you"
                                                                " would like to track the distribution shapes for convergence")
    parser.add_argument("--convergence_limits", type=json.loads, help="dictionary that can contain limits for convergence params")
    parser.add_argument("--pop_select", help="Used in combination with the specified final_kstar1 and final_kstar2 values"
                                             " to select the subpopulation of interest from the evolved population")
    parser.add_argument("--match", type=float, help="provides the tolerance for the convergence calculation")
    parser.add_argument("--apply_convergence_limits", type=str2bool, nargs='?',
                        const=True, default=False, help="filters the evolved binary population to contain"
                                                        " only the binaries that satsify the convergence limits")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Run in Verbose Mode")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-n", "--nproc",
                        help="number of processors", type=int, default=1)
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")

    args = parser.parse_args(remaining_argv)

    if len(args.final_kstar1) > 2 or len(args.final_kstar2) > 2:
        raise parser.error('final kstar1 and final kstar2 '
                           'must be either a single value or '
                           'a range between two values.')

    if (len(args.final_kstar1) == 2):
        if (args.final_kstar1[0] >= args.final_kstar1[1]):
            raise parser.error('Range provided for final-kstar1 invalid')

    if (len(args.final_kstar2) == 2):
        if (args.final_kstar2[0] >= args.final_kstar2[1]):
            raise parser.error('Range provided for final-kstar2 invalid')

    if (len(args.final_kstar2) == 1) and (len(args.final_kstar1) == 1):
        if (args.final_kstar2 > args.final_kstar1):
            raise parser.error('final-kstar1 must be greater than or equal to '
                               'final-kstar2.')

    return args

if __name__ == '__main__':
    print('parsing arguments')
    args = parse_commandline()
    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.nproc)
    if isinstance(pool, MPIPool):
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        nproc = len(pool.workers)
    else:
        nproc = args.nproc
    #Parse inifile
    BSEDict, seed_int, filters, convergence, sampling = utils.parse_inifile(args.inifile)

    # we now overwrite the inifile values with what was specified from the command line
    # (which could mean not overwriting anything at all because they are populated
    # by default from the inifile).
    for argument in vars(args):
        if argument in filters.keys():
            if filters[argument] != getattr(args, argument):
                warnings.warn("You are overriding the inifile value of {0}={1} "
                              "with {0}={2} from the commandline".format(argument, filters[argument], getattr(args, argument)))
                filters[argument] = getattr(args, argument)

        if argument in convergence.keys():
            if convergence[argument] != getattr(args, argument):
                warnings.warn("You are overriding the inifile value of {0}={1} "
                              "with {0}={2} from the commandline".format(argument, convergence[argument], getattr(args, argument)))
                convergence[argument] = getattr(args, argument)

        if argument in sampling.keys():
            if (sampling[argument] == "independent") or (getattr(args, argument) == "independent"):
               for model in ["primary_model", "porb_model", "ecc_model", "binfrac_model"]:
                    if (model not in sampling.keys()) and not (getattr(args, model)):
                        raise ValueError("You have selected the {0} sampler "
                                         "but not specified a model for {1} "
                                         "in the inifile or command line".format(sampling[argument], model))
            if sampling[argument] != getattr(args, argument):
                warnings.warn("You are overriding the inifile value of {0}={1} "
                              "with {0}={2} from the commandline".format(argument, sampling[argument], getattr(args, argument)))
                sampling[argument] = getattr(args, argument)

        if argument == 'seed':
            if getattr(args, argument) != seed_int:
                warnings.warn("You are overriding the inifile value of {0}={1} "
                              "with {0}={2} from the commandline".format(argument, seed_int, getattr(args, argument)))
                seed_int = getattr(args, argument)

    # Check that the values in BSEDict, filters, and convergence are valid
    utils.error_check(BSEDict, filters, convergence, sampling)

    if seed_int != 0:
        np.random.seed(seed_int)
    else:
        np.random.seed(0)

    # Set up final_kstar1 and final_kstar2 strings for saved data files
    if len(args.final_kstar1) == 2:
        kstar1_range = np.arange(args.final_kstar1[0], args.final_kstar1[1]+1)
        kstar1_range_string = str(int(args.final_kstar1[0]))+'_'+str(int(args.final_kstar1[1]))
    else:
        kstar1_range = args.final_kstar1
        kstar1_range_string = str(int(args.final_kstar1[0]))

    if len(args.final_kstar2) == 2:
        kstar2_range = np.arange(args.final_kstar2[0], args.final_kstar2[1]+1)
        kstar2_range_string = str(int(args.final_kstar2[0]))+'_'+str(int(args.final_kstar2[1]))
    else:
        kstar2_range = args.final_kstar2
        kstar2_range_string = str(int(args.final_kstar2[0]))

    # Open the hdf5 file to store the fixed population data
    try:
        dat_store = pd.HDFStore('dat_kstar1_{0}_kstar2_{1}_SFstart_{2}_SFduration_{3}_metallicity_{4}.h5'.format(kstar1_range_string, kstar2_range_string, sampling['SF_start'], sampling['SF_duration'], sampling['metallicity']))
        log_file = open('log_kstar1_{0}_kstar2_{1}_SFstart_{2}_SFduration_{3}_metallicity_{4}.txt'.format(kstar1_range_string, kstar2_range_string, sampling['SF_start'], sampling['SF_duration'], sampling['metallicity']), 'a')
        log_file.write('There are already: '+str(conv_save.shape[0])+' '+kstar1_range_string+'_'+kstar2_range_string+' binaries evolved\n')
        log_file.write('\n')
        total_mass_singles = np.max(pd.read_hdf(dat_store, 'mass_singles'))[0]
        total_mass_binaries = np.max(pd.read_hdf(dat_store, 'mass_binaries'))[0]
        total_mass_stars = np.max(pd.read_hdf(dat_store, 'mass_stars'))[0]
        total_n_singles = np.max(pd.read_hdf(dat_store, 'n_singles'))[0]
        total_n_binaries = np.max(pd.read_hdf(dat_store, 'n_binaries'))[0]
        total_n_stars = np.max(pd.read_hdf(dat_store, 'n_stars'))[0]
        idx = int(np.max(pd.read_hdf(dat_store, 'idx'))[0])
    except:
        dat_store = pd.HDFStore('dat_kstar1_{0}_kstar2_{1}_SFstart_{2}_SFduration_{3}_metallicity_{4}.h5'.format(kstar1_range_string, kstar2_range_string, sampling['SF_start'], sampling['SF_duration'], sampling['metallicity']))
        total_mass_singles = 0
        total_mass_binaries = 0
        total_mass_stars = 0
        total_n_singles = 0
        total_n_binaries = 0
        total_n_stars = 0
        idx = 0
        log_file = open('log_kstar1_{0}_kstar2_{1}_SFstart_{2}_SFduration_{3}_metallicity_{4}.txt'.format(kstar1_range_string, kstar2_range_string, sampling['SF_start'], sampling['SF_duration'], sampling['metallicity']), 'w')


    # save configuration settings to output file
    configuration_settings = {'BSEDict' : BSEDict, 'filters' : filters,
                              'convergence' : convergence, 'sampling' : sampling}

    for k, v in configuration_settings.items():
        for k1, v1 in v.items():
            dat_store.put('config/{0}/{1}/'.format(k, k1), pd.Series(v1))
    dat_store.put('config/rand_seed/', pd.Series(seed_int))

    Nstep = idx - np.mod(idx, args.Nstep)

    #Main loop: keep simulating until we have sampled the amount of mass we would like
    print('Running main loop')
    t = time.time()

    while total_mass_stars < args.total_mass:
        # Set random seed such that each iteration gets a unique, determinable seed
        rand_seed = seed_int + Nstep
        np.random.seed(rand_seed)

        # Select the initial binary sample method from user input
        if sampling['sampling_method'] == 'independent':
            if hasattr(args,'qmin'):
                init_samp_list = InitialBinaryTable.sampler(format_ = sampling['sampling_method'],
                                                            final_kstar1 = kstar1_range,
                                                            final_kstar2 = kstar2_range,
                                                            binfrac_model = args.binfrac_model,
                                                            primary_model = args.primary_model,
                                                            ecc_model = args.ecc_model,
                                                            porb_model = args.porb_model,
                                                            keep_singles = args.keep_singles,
                                                            SF_start = sampling['SF_start'],
                                                            SF_duration = sampling['SF_duration'],
                                                            met = sampling['metallicity'],
                                                            size = args.Nstep,
                                                            qmin = args.qmin,
                                                            params = args.inifile)
            elif hasattr(args,'m2_min'):
                init_samp_list = InitialBinaryTable.sampler(format_ = sampling['sampling_method'],
                                                            final_kstar1 = kstar1_range,
                                                            final_kstar2 = kstar2_range,
                                                            binfrac_model = args.binfrac_model,
                                                            primary_model = args.primary_model,
                                                            ecc_model = args.ecc_model,
                                                            porb_model = args.porb_model,
                                                            keep_singles = args.keep_singles,
                                                            SF_start = sampling['SF_start'],
                                                            SF_duration = sampling['SF_duration'],
                                                            met = sampling['metallicity'],
                                                            size = args.Nstep,
                                                            m2_min = args.m2_min,
                                                            params = args.inifile,
                                                            )
            else:
                raise ValueError("You must specify either qmin or m2_min in the",
                                 " inifile if you are using the independent sampler")
            IBT, mass_singles, mass_binaries, n_singles, n_binaries = init_samp_list

        if sampling['sampling_method'] == 'multidim':
            init_samp_list = InitialBinaryTable.sampler(format_ = sampling['sampling_method'],
                                                        final_kstar1 = kstar1_range,
                                                        final_kstar2 = kstar2_range,
                                                        keep_singles = args.keep_singles,
                                                        rand_seed = rand_seed,
                                                        nproc = args.nproc,
                                                        SF_start = sampling['SF_start'],
                                                        SF_duration = sampling['SF_duration'],
                                                        met = sampling['metallicity'],
                                                        size = args.Nstep,
                                                        pool=pool)
            IBT, mass_singles, mass_binaries, n_singles, n_binaries = init_samp_list

        # Log the total sampled mass from the initial binary sample
        # for future Galactic occurence rate calculation
        total_mass_singles += mass_singles
        total_mass_binaries += mass_binaries
        total_mass_stars += mass_singles + mass_binaries
        total_n_singles += n_singles
        total_n_binaries += n_binaries
        total_n_stars += n_singles + 2*n_binaries

        # Now that we have all these initial conditions
        # let's create an Evolve class and evolve these systems

        # check what kind of time resolution for the bcm array the user specified

        # assume none
        dtp = IBT['tphysf'].values

        # check
        if isinstance(filters['timestep_conditions'], str):
            dtp_inifile = filters['timestep_conditions'].split('=')[-1]
            try:
                dtp = float(dtp_inifile)
            except:
                pass
            filters['timestep_conditions'] = []

        # Create a pool
        bpp, bcm, initCond, kick_info = Evolve.evolve(initialbinarytable=IBT,
                                                      pool=pool,
                                                      BSEDict=BSEDict,
                                                      idx=idx,
                                                      dtp=dtp,
                                                      timestep_conditions=filters['timestep_conditions'],)

        # get any nans and pull them out for now
        nans = np.isnan(bpp.sep)
        if nans.any():
            nan_bin_nums = np.unique(bpp[nans]["bin_num"].values)
            initCond_nan = initCond.loc[initCond.bin_num.isin(nan_bin_nums)]
            dat_store.append("nan_initC", initCond_nan)
            log_file.write(f"There are {len(nan_bin_nums)} NaNs stored in the datfile with key: 'nan_initC'")
            log_file.write(f"These NaNs likely arise because you have pts1 = 0.001, try running with pts1 = 0.01")

            bcm = bcm.loc[~bcm.bin_num.isin(nan_bin_nums)]
            bpp = bpp.loc[~bpp.bin_num.isin(nan_bin_nums)]
            initCond = initCond.loc[~initCond.bin_num.isin(nan_bin_nums)]
            kick_info = kick_info.loc[~kick_info.bin_num.isin(nan_bin_nums)]

        # Keep track of the index
        idx = int(bcm.bin_num.max()+1)

        # If dtp is not set, filter out first timestep in bcm
        if np.all(dtp == IBT['tphysf'].values):
            bcm = bcm.loc[bcm['tphys'].isin(dtp)]

        #Now, we need to filter bcm and bpp so that we are only keeping DWD LISA sources
        porb_lim = 2 / (10**-5) * 1 / (60 * 60 * 24)
        query = f'(kstar_1.isin(@kstar1_range)) & (kstar_2.isin(@kstar2_range)) & (porb < @porb_lim)'
        bcm = bcm.query(query)
        bin_nums = bcm.bin_num.unique()

        bpp = bpp.loc[bpp.bin_num.isin(bin_nums)]
        initCond = initCond.loc[initCond.bin_num.isin(bin_nums)]
        kick_info = kick_info.loc[kick_info.bin_num.isin(bin_nums)]


        #Save everything!
        # write the data and the logs!
        mass_list = [total_mass_singles, total_mass_binaries, total_mass_stars]
        n_list = [total_n_singles, total_n_binaries, total_n_stars]
        m_keys = ["mass_singles", "mass_binaries", "mass_stars"]
        n_keys = ["n_singles", "n_binaries", "n_stars"]
        for m_write, m_key, n_write, n_key in zip(mass_list, m_keys, n_list, n_keys):
            # save the total_sampled_mass so far
            dat_store.append(m_key, pd.DataFrame([m_write]))
            dat_store.append(n_key, pd.DataFrame([n_write]))

        # Save the bcm dataframe
        dat_store.append("bcm", bcm)

        # Save the bpp dataframe
        dat_store.append("bpp", bpp)

        # Save the initial binaries
        # ensure that the index corresponds to bin_num
        dat_store.append("initCond", initCond.set_index("bin_num", drop=False))

        # Save the kick_info dataframe
        dat_store.append("kick_info", kick_info)

        # Save the index
        dat_store.append("idx", pd.DataFrame([idx]))

        #How fast are we going?
        ellapsed_time = time.time() - t
        estimated_time_remaining = args.total_mass / total_mass_stars * ellapsed_time - ellapsed_time
        print(f'[*] {round(total_mass_stars,2)} mass sampled of {round(args.total_mass,2)} in {round(ellapsed_time,2)}. Estimated remaining time: {round(estimated_time_remaining,2)}')
        log_file.write(f'[*] {round(total_mass_stars,2)} mass sampled of {round(args.total_mass,2)} in {round(ellapsed_time,2)}. Estimated remaining time: {round(estimated_time_remaining,2)}')

        Nstep += 1

    # Close the data storage file
    dat_store.close()

    log_file.write('All done friend!')
    log_file.close()
