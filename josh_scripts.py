# Some scripts that are helpful for me
import math
import general_purpose_python_modules.reproduce_system as rs
from types import SimpleNamespace
from astropy import units
from stellar_evolution.library_interface import MESAInterpolator
#MESAInterpolator.set_quantity_lower_limit('iconv', 1e-5)
from stellar_evolution.manager import StellarEvolutionManager
import numpy
from orbital_evolution.transformations import phase_lag
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from multiprocessing_util import setup_process,setup_process_map
import matplotlib.pyplot as plt
#from orbital_evolution.binary import Binary
#from orbital_evolution.star_interface import EvolvingStar
from POET.solver import poet_solver

from multiprocessing import Manager
thelock = Manager().Lock()

import logging

def masses(mttot,mrat):
    """Returns the masses of the two bodies given the total mass and mass ratio"""
    m1 = mttot/(1+mrat)
    m2 = mttot - m1
    return m1,m2

def e(sinC,cosC):
    c=sinC
    d=cosC
    return c/math.sin(math.atan(c/d))

def age(t):
    return (10**t) / (10**9)

def retry_system():
    print("Retrying system...")

    disk_dissipation_age = 0.01 * units.Gyr
    lgQ_inertial_boost = 0.0
    lgQ_inertial_sharpness = 10.0
    primary_disk_lock_period = 5.0
    primary_wind_strength = 0.17
    primary_wind_saturation = 2.45
    primary_core_envelope_coupling_timescale = 0.01# * units.Gyr
    secondary_disk_lock_period = None#5.0
    secondary_wind_strength = 0.17
    secondary_wind_saturation = 2.45
    secondary_core_envelope_coupling_timescale = 0.01# * units.Gyr

    # Variables
    initial_porb = 21.72109769815944 * units.day
    lgQ_min = 5#5.718532564697535
    lgQ_break_period = 2.9578562443466514 * units.day
    lgQ_powerlaw = -4#-4.44462284528867
    initial_eccentricity = 'solve'#0.8
    age = 0.3#0.2747635194868908
    feh = 0.04#0.03627395537563404
    orbital_period = 7.964402489325128
    primary_mass = 0.7#0.6686791407972886
    secondary_mass = 0.6#0.6951645709715627
    cmd_primary_radius = 0.61#0.609697911391692
    cmd_secondary_radius = 0.62#00039277087531

    solve = True#False

    # parameters = SimpleNamespace(
    #               disk_dissipation_age = 0.01 * units.Gyr,
    #               ,
    #               primary_wind_strength = 0.17,
    #               primary_wind_saturation = 2.45,
    #               primary_core_envelope_coupling_timescale = 0.01 * units.Gyr,
    #               initial_eccentricity = 0.8,
    #               secondary_wind_strength = 0.17,
    #               secondary_wind_saturation = 2.45,
    #               secondary_core_envelope_coupling_timescale = 0.01 * units.Gyr,
    #               ,
    #               precision = 1e-05)

    # From parameters
    disk_period = primary_disk_lock_period * units.day
    secondary_is_star = True

    star_dissipation = dict(
        spin_frequency_breaks=None,
        spin_frequency_powers=numpy.array([0.0]),
        reference_phase_lag=phase_lag(
            lgQ_min+max(lgQ_inertial_boost,0.0)
        ),
        inertial_mode_enhancement=10.0**(lgQ_inertial_boost),
        inertial_mode_sharpness=lgQ_inertial_sharpness
    )
    break_frequency = (2.0 * numpy.pi
                        /
                        lgQ_break_period.to_value(units.day))
    powerlaw = lgQ_powerlaw

    if powerlaw > 0:
        star_dissipation['tidal_frequency_breaks'] = numpy.array([
            2.0 * numpy.pi / 50.0,
            break_frequency
        ])
        star_dissipation['tidal_frequency_powers'] = numpy.array([
            0.0,
            powerlaw,
            0.0
        ])
        star_dissipation['reference_phase_lag'] *= numpy.power(
            star_dissipation['tidal_frequency_breaks'][0]
            /
            star_dissipation['tidal_frequency_breaks'][1],
            powerlaw
        )
    else:
        star_dissipation['tidal_frequency_breaks'] = numpy.array([
            break_frequency
        ])
        star_dissipation['tidal_frequency_powers'] = numpy.array([
            0.0,
            powerlaw
        ])

    dissipation = dict(
                        primary = star_dissipation,#dict(
                        #     #spin_frequency_breaks = None,
                        #     #spin_frequency_powers = array([0.]),
                        #     #reference_phase_lag = 2.81578188e-07,
                        #     #inertial_mode_enhancement = 1.,
                        #     #inertial_mode_sharpness = 10.,
                        #     #tidal_frequency_breaks = array([0.38596437]),
                        #     #tidal_frequency_powers = array([ 0.        , -1.99262346])
                        #     lgQ_min = 5.43559054457393,
                        #     lgQ_inertial_boost = 0.0,
                        #     lgQ_inertial_sharpness = 10.0,
                        #     lgQ_break_period = 18.323358198954484,
                        #     lgQ_powerlaw = 2.7489352220413537
                        # ),
                        secondary = dict(
                            # spin_frequency_breaks = None,
                            # spin_frequency_powers = array([0.]),
                            # reference_phase_lag = 2.81578188e-07,
                            # inertial_mode_enhancement = 1.,
                            # inertial_mode_sharpness = 10.,
                            # tidal_frequency_breaks = array([0.38596437]),
                            # tidal_frequency_powers = array([ 0.        , -1.99262346])
                            )
                        )
    system = SimpleNamespace(
                        orbital_period= orbital_period * units.day,
                        age= age * units.Gyr,
                        eccentricity= 0.2561160543396648,#eccentricity,
                        primary_mass= primary_mass * units.solMass,
                        secondary_mass= secondary_mass * units.solMass,
                        feh= feh,
                        Rprimary= cmd_primary_radius * units.solRad,
                        Rsecondary= cmd_secondary_radius * units.solRad
        )

    dissipation['secondary'] = dissipation['primary']
    print('time for interpolatorer')
    interpolator = StellarEvolutionManager('/home/vortebo/ctime/poet/stellar_evolution_interpolators').get_interpolator_by_name(
            'default'
        )
    print('interpolatorer overed')
    print(interpolator)

    print('loading')
    orbital_evolution_library.prepare_eccentricity_expansion(
        '/home/vortebo/eccentricity_expansion_coef_O400.sqlite'.encode('ascii'),
        1e-4,
        True,
        True
    )
    print('ready')

    #Call
    print('run')
    gary = rs.find_evolution(system = system,
        interpolator = interpolator,
        dissipation = dissipation,
        initial_porb = initial_porb,
        initial_eccentricity = initial_eccentricity,
        #initial_obliquity = 3.0,
        disk_period=disk_period,
        disk_dissipation_age=disk_dissipation_age,
        primary_wind_strength=primary_wind_strength,
        primary_wind_saturation=primary_wind_saturation,
        primary_core_envelope_coupling_timescale=primary_core_envelope_coupling_timescale * units.Gyr,
        secondary_wind_strength=secondary_wind_strength,
        secondary_wind_saturation=secondary_wind_saturation,
        secondary_core_envelope_coupling_timescale=secondary_core_envelope_coupling_timescale * units.Gyr,
        secondary_disk_period=secondary_disk_lock_period,
        #orbital_period_tolerance=1e-6,
        #eccentricity_tolerance=1e-6,
        #obliquity_tolerance=1e-6,
        #period_search_factor=2.0,
        #scaled_period_guess=1.0,
        #eccentricity_upper_limit=0.8,
        solve=solve,
        #max_iterations=49,
        secondary_is_star=secondary_is_star,
        precision = 1e-5,
        eccentricity_expansion_fname = '/home/vortebo/eccentricity_expansion_coef_O400.sqlite'.encode('ascii')
    )
    print(gary)

def get_dissipation(lgQ_min,lgQ_inertial_boost,lgQ_inertial_sharpness,lgQ_break_period,lgQ_powerlaw,age_breaks=None):
    star_dissipation = dict(
        spin_frequency_breaks=None,
        spin_frequency_powers=numpy.array([0.0]),
        age_breaks=age_breaks,
        reference_phase_lags=phase_lag(
        #reference_phase_lag=phase_lag(
            lgQ_min+max(lgQ_inertial_boost,0.0)
        ),
        inertial_mode_enhancement=10.0**(lgQ_inertial_boost),
        inertial_mode_sharpness=lgQ_inertial_sharpness
    )
    break_frequency = (2.0 * numpy.pi
                        /
                        lgQ_break_period.to_value(units.day))
    powerlaw = lgQ_powerlaw
    if powerlaw > 0:
        star_dissipation['tidal_frequency_breaks'] = numpy.array([
            2.0 * numpy.pi / 50.0,
            break_frequency
        ])
        star_dissipation['tidal_frequency_powers'] = numpy.array([
            0.0,
            powerlaw,
            0.0
        ])
        star_dissipation['reference_phase_lags'] *= numpy.power(
        #star_dissipation['reference_phase_lag'] *= numpy.power(
            star_dissipation['tidal_frequency_breaks'][0]
            /
            star_dissipation['tidal_frequency_breaks'][1],
            powerlaw
        )
    else:
        star_dissipation['tidal_frequency_breaks'] = numpy.array([
            break_frequency
        ])
        star_dissipation['tidal_frequency_powers'] = numpy.array([
            0.0,
            powerlaw
        ])

    star_dissipation['reference_phase_lags'] = numpy.array([
            star_dissipation['reference_phase_lags']
        ])

    return star_dissipation

def plot_evolution(evolution,name):

    primary_iconv_keys = [i for i in range(len(evolution.primary_iconv)) if not math.isnan(evolution.primary_iconv[i])]
    secondary_iconv_keys = [i for i in range(len(evolution.secondary_iconv)) if not math.isnan(evolution.secondary_iconv[i])]

    primary_env_omg = 2 * math.pi * numpy.array(evolution.primary_iconv) / numpy.array(evolution.primary_envelope_angmom)
    primary_core_omg = 2 * math.pi * numpy.array(evolution.primary_irad) / numpy.array(evolution.primary_core_angmom)
    secondary_env_omg = 2 * math.pi * numpy.array(evolution.secondary_iconv) / numpy.array(evolution.secondary_envelope_angmom)
    secondary_core_omg = 2 * math.pi * numpy.array(evolution.secondary_irad) / numpy.array(evolution.secondary_core_angmom)

    #name = '/home/vortebo/ctime/josh_output/'+name
    
    plt.plot(evolution.age,evolution.orbital_period)
    plt.savefig(name+'_age_porb.png')
    plt.clf()
    plt.plot(evolution.age,evolution.eccentricity)
    plt.savefig(name+'_age_eccentricity.png')
    plt.clf()
    plt.plot(evolution.age,evolution.primary_radius)
    plt.savefig(name+'_age_primary_radius.png')
    plt.clf()
    plt.plot(age,evolution.primary_lum)
    plt.savefig(name+'_age_primary_lum.png')
    plt.clf()
    plt.plot(age,evolution.primary_iconv)
    plt.plot(age,evolution.primary_irad)
    added = numpy.array(evolution.primary_iconv[primary_iconv_keys[0]:]) + numpy.array(evolution.primary_irad[primary_iconv_keys[0]:])
    plt.plot(age[primary_iconv_keys[0]:],added)
    plt.savefig(name+'_age_primary_i.png')
    plt.clf()
    plt.plot(age,evolution.secondary_radius)
    plt.savefig(name+'_age_secondary_radius.png')
    plt.clf()
    plt.plot(age,evolution.secondary_lum)
    plt.savefig(name+'_age_secondary_lum.png')
    plt.clf()
    plt.plot(age,evolution.secondary_iconv)
    plt.plot(age,evolution.secondary_irad)
    added = numpy.array(evolution.secondary_iconv[secondary_iconv_keys[0]:]) + numpy.array(evolution.secondary_irad[secondary_iconv_keys[0]:])
    plt.plot(age[secondary_iconv_keys[0]:],added)
    plt.savefig(name+'_age_secondary_i.png')
    plt.clf()
    plt.plot(age,primary_env_omg)
    plt.plot(age,primary_core_omg)
    plt.plot(evolution.age,evolution.orbital_period)
    plt.savefig(name+'_age_primary_omg.png')
    plt.clf()
    plt.plot(age,secondary_env_omg)
    plt.plot(age,secondary_core_omg)
    plt.plot(evolution.age,evolution.orbital_period)
    plt.savefig(name+'_age_secondary_omg.png')
    plt.clf()

def test_age_split():

    # tests would be something like:
    # 1. Run with PMS dissipation and no MS dissipation
    # 2. Run with MS dissipation and no PMS dissipation
    # and then hopefully we will see the orbit evolving only in the relevant regime

    # and one more very useful test occurs to me:
    # run an evolution with no age breaks and then one where dissipation just turns off after some time
    # ideally you will see the second one following the first exactly up to the age when it turns off and then only stellar spin should change
    
    disk_dissipation_age = 0.01 * units.Gyr
    primary_disk_lock_period = 5.0
    primary_wind_strength = 0.17
    primary_wind_saturation = 2.45
    primary_core_envelope_coupling_timescale = 0.01 * units.Gyr
    secondary_disk_lock_period = None#5.0
    secondary_wind_strength = 0.17
    secondary_wind_saturation = 2.45
    secondary_core_envelope_coupling_timescale = 0.01 * units.Gyr
    # From parameters
    disk_period = primary_disk_lock_period * units.day
    secondary_is_star = True
    # Variables
    initial_porb = 50 * units.day
    initial_eccentricity = 0.8
    age = 5.89
    feh = 0.007
    orbital_period = 27.3
    primary_mass = 0.94
    secondary_mass = 0.60
    cmd_primary_radius = 0.609697911391692
    cmd_secondary_radius = 0.6200039277087531
    # This is just true
    solve = False
    system = SimpleNamespace(
                        orbital_period= orbital_period * units.day,
                        age= age * units.Gyr,
                        eccentricity= 0.44,#eccentricity,
                        primary_mass= primary_mass * units.solMass,
                        secondary_mass= secondary_mass * units.solMass,
                        feh= feh,
                        Rprimary= cmd_primary_radius * units.solRad,
                        Rsecondary= cmd_secondary_radius * units.solRad
        )
    
    print('time for interpolatorer')
    interpolator = StellarEvolutionManager('/home/vortebo/ctime/poet/stellar_evolution_interpolators').get_interpolator_by_name(
            'default'
        )
    print('interpolatorer overed')
    print('loading')
    orbital_evolution_library.prepare_eccentricity_expansion(
        '/home/vortebo/eccentricity_expansion_coef_O400.sqlite'.encode('ascii'),
        #'previouslyOn/makeSomeLessGoodDats/dontcompare.db'.encode('ascii'),
        1e-4,
        True,
        True
    )
    print('ready')

    # Test 1: No dissipation at any time
    lgQ_inertial_boost = 0.0
    lgQ_inertial_sharpness = 10.0
    lgQ_min = numpy.array([20.0,20.0])
    lgQ_break_period = 2.9578562443466514 * units.day
    lgQ_powerlaw = -4.44462284528867
    age_breaks = numpy.array([1])
    star_dissipation = get_dissipation(lgQ_min,lgQ_inertial_boost,lgQ_inertial_sharpness,lgQ_break_period,lgQ_powerlaw,age_breaks)
    dissipation = dict(
                        primary = star_dissipation,
                        secondary = star_dissipation
                        )

    parameters_1 = dict(
        system = system,
        interpolator = interpolator,
        dissipation = dissipation,
        initial_porb = initial_porb,
        initial_eccentricity = initial_eccentricity,
        disk_period=disk_period,
        disk_dissipation_age=disk_dissipation_age,
        primary_wind_strength=primary_wind_strength,
        primary_wind_saturation=primary_wind_saturation,
        primary_core_envelope_coupling_timescale=primary_core_envelope_coupling_timescale,
        secondary_wind_strength=secondary_wind_strength,
        secondary_wind_saturation=secondary_wind_saturation,
        secondary_core_envelope_coupling_timescale=secondary_core_envelope_coupling_timescale,
        secondary_disk_period=secondary_disk_lock_period,
        solve=solve,
        secondary_is_star=secondary_is_star,
        precision = 1e-5,
        eccentricity_expansion_fname = '/home/vortebo/eccentricity_expansion_coef_O400.sqlite'.encode('ascii')#,
        #plotname = 'testbreak1'
    )

    # Test 2: Same dissipation for all time
    lgQ_min = numpy.array([5.0,5.0])
    age_breaks = numpy.array([1])
    star_dissipation = get_dissipation(lgQ_min,lgQ_inertial_boost,lgQ_inertial_sharpness,lgQ_break_period,lgQ_powerlaw,age_breaks)
    dissipation = dict(
                        primary = star_dissipation,
                        secondary = star_dissipation
                        )
    parameters_2 = parameters_1
    parameters_2['dissipation'] = dissipation
    #parameters_2['plotname'] = 'testbreak2'

    # Test 3: Dissipation just disappears after some time
    lgQ_min = numpy.array([5.0,20.0])
    age_breaks = numpy.array([2.3])
    star_dissipation = get_dissipation(lgQ_min,lgQ_inertial_boost,lgQ_inertial_sharpness,lgQ_break_period,lgQ_powerlaw,age_breaks)
    dissipation = dict(
                        primary = star_dissipation,
                        secondary = star_dissipation
                        )
    parameters_3 = parameters_1
    parameters_3['dissipation'] = dissipation
    #parameters_3['plotname'] = 'testbreak3'

    # Test 4: PMS but no MS dissipation
    lgQ_min = numpy.array([5.0,20.0])
    age_breaks = numpy.array([0.001])
    star_dissipation = get_dissipation(lgQ_min,lgQ_inertial_boost,lgQ_inertial_sharpness,lgQ_break_period,lgQ_powerlaw,age_breaks)
    dissipation = dict(
                        primary = star_dissipation,
                        secondary = star_dissipation
                        )
    parameters_4 = parameters_1
    parameters_4['dissipation'] = dissipation
    #parameters_4['plotname'] = 'testbreak4'

    # Test 5: MS but no PMS dissipation
    lgQ_min = numpy.array([20.0,5.0])
    star_dissipation = get_dissipation(lgQ_min,lgQ_inertial_boost,lgQ_inertial_sharpness,lgQ_break_period,lgQ_powerlaw,age_breaks)
    dissipation = dict(
                        primary = star_dissipation,
                        secondary = star_dissipation
                        )
    parameters_5 = parameters_1
    parameters_5['dissipation'] = dissipation
    #parameters_5['plotname'] = 'testbreak5'

    lgQ_min = 5.0#numpy.array([5.0])
    age_breaks = None
    star_dissipation = get_dissipation(lgQ_min,lgQ_inertial_boost,lgQ_inertial_sharpness,lgQ_break_period,lgQ_powerlaw,age_breaks)
    dissipation = dict(
                        primary = star_dissipation,
                        secondary = star_dissipation
                        )
    parameters_6 = parameters_1
    parameters_6['dissipation'] = dissipation

    aa_1 = dict(
        parameters = parameters_1,
        plotname = 'testbreak1'
    )
    aa_2 = dict(
        parameters = parameters_2,
        plotname = 'testbreak2'
    )
    aa_3 = dict(
        parameters = parameters_3,
        plotname = 'testbreak3'
    )
    aa_4 = dict(
        parameters = parameters_4,
        plotname = 'testbreak4'
    )
    aa_5 = dict(
        parameters = parameters_5,
        plotname = 'testbreak5'
    )
    aa_6 = dict(
        parameters = parameters_6,
        plotname = 'testbreak6'
    )
    #load the smaller database plz

    # Use multiprocessing to run all the tests at once
    #from multiprocessing import Pool

    #with Pool(5) as p:
    #    p.map(run_evolution,[aa_1,aa_2,aa_3,aa_4,aa_5])
    print(aa_6['parameters']['dissipation']['primary']['reference_phase_lag'])
    run_evolution(aa_6)

    print('Success?')
    #pool = Pool(processes=5)
    #pool.map(run_evolution, [parameters_1,parameters_2,parameters_3,parameters_4,parameters_5])
    #pool.close()

def run_evolution(parameters):
    evolution = rs.find_evolution(**parameters['parameters'])[0]
    plot_evolution(evolution,parameters['plotname'])

def read_fits():

    import astropy.io.fits as fits
    filename = "/home/vortebo/ctime/failed_solutions/solution_2023-11-09_16-03-22.fits"

    # Open the file and print out all the data
    with fits.open(filename) as data:
        print(data.info())
        for i in range(len(data)):
            print("Header: ", i)
            print(data[i].header)
            print("Data: ", i)
            if data[i].data is not None:
                for n in range(len(data[i].data)):
                    print(data[i].data[n])

def emergency_plotter():
    import astropy.io.fits as fits
    filename = "/home/vortebo/ctime/failed_solutions/solution_2023-11-09_16-20-46.fits"

    age = []
    porb = []
    eccentricity = []
    primary_radius = []
    primary_lum = []
    primary_iconv = []
    primary_irad = []
    secondary_radius = []
    secondary_lum = []
    secondary_iconv = []
    secondary_irad = []
    primary_L_env = []
    primary_L_core = []
    secondary_L_env = []
    secondary_L_core = []

    # Open the file and print out all the data
    with fits.open(filename) as data:
        for n in range(len(data[1].data)):
            out = data[1].data[n]
            age.append(out[0])
            porb.append(out[1])
            eccentricity.append(out[2])
            primary_radius.append(out[3])
            primary_lum.append(out[4])
            primary_iconv.append(out[5])
            primary_irad.append(out[6])
            secondary_radius.append(out[7])
            secondary_lum.append(out[8])
            secondary_iconv.append(out[9])
            secondary_irad.append(out[10])
            primary_L_env.append(out[11])
            primary_L_core.append(out[12])
            secondary_L_env.append(out[13])
            secondary_L_core.append(out[14])
    
    primary_iconv_keys = [i for i in range(len(primary_iconv)) if not math.isnan(primary_iconv[i])]
    secondary_iconv_keys = [i for i in range(len(secondary_iconv)) if not math.isnan(secondary_iconv[i])]

    primary_env_omg = 2 * math.pi * numpy.array(primary_iconv) / numpy.array(primary_L_env)
    primary_core_omg = 2 * math.pi * numpy.array(primary_irad) / numpy.array(primary_L_core)
    secondary_env_omg = 2 * math.pi * numpy.array(secondary_iconv) / numpy.array(secondary_L_env)
    secondary_core_omg = 2 * math.pi * numpy.array(secondary_irad) / numpy.array(secondary_L_core)

    plt.plot(age,porb)
    plt.savefig('age_porb.png')
    #plt.legend('porb')
    plt.clf()
    plt.plot(age,eccentricity)
    plt.savefig('age_eccentricity.png')
    plt.clf()
    plt.plot(age,primary_radius)
    plt.savefig('age_primary_radius.png')
    plt.clf()
    plt.plot(age,primary_lum)
    plt.savefig('age_primary_lum.png')
    plt.clf()
    plt.plot(age,primary_iconv)
    plt.plot(age,primary_irad)
    added = numpy.array(primary_iconv[primary_iconv_keys[0]:]) + numpy.array(primary_irad[primary_iconv_keys[0]:])
    plt.plot(age[primary_iconv_keys[0]:],added)
    #plt.legend('primary_iconv','primary_irad','primary_iconv+irad')
    plt.savefig('age_primary_i.png')
    plt.clf()
    plt.plot(age,secondary_radius)
    plt.savefig('age_secondary_radius.png')
    plt.clf()
    plt.plot(age,secondary_lum)
    plt.savefig('age_secondary_lum.png')
    plt.clf()
    plt.plot(age,secondary_iconv)
    plt.plot(age,secondary_irad)
    added = numpy.array(secondary_iconv[secondary_iconv_keys[0]:]) + numpy.array(secondary_irad[secondary_iconv_keys[0]:])
    plt.plot(age[secondary_iconv_keys[0]:],added)
    #plt.legend('secondary_iconv','secondary_irad','secondary_iconv+irad')
    plt.savefig('age_secondary_i.png')
    plt.clf()
    plt.plot(age,primary_env_omg)
    plt.plot(age,primary_core_omg)
    plt.plot(age,porb)
    #plt.legend('primary_env_omg','primary_core_omg')
    plt.savefig('age_primary_omg.png')
    plt.clf()
    plt.plot(age,secondary_env_omg)
    plt.plot(age,secondary_core_omg)
    plt.plot(age,porb)
    #plt.legend('secondary_env_omg','secondary_core_omg')
    plt.savefig('age_secondary_omg.png')
    plt.clf()

def trainnn():
    path = '/home/vortebo/ctime/ayeye'
    params = {
                "type": 'blank',
                "epochs": 30,
                "batch_size": 100,
                "verbose": 2,
                "retrain": False,
                "threshold": 20,
                "path_to_store": path,
                "version": None,
                "features": None
            }
    systems = ['12356914','3348093']

    # for system in systems:
    #     params['version'] = system
    #     # Abc
    #     params['type'] = '1d_period'
    #     params['features'] = [True, True, True, True, True, True, True, True, True, True]
    #     #X_test= numpy.array([5.0, 10.0, 0.5, 2.0, 0.06, 69.0, 1.0, 0.7, 1.0, 0.8])
    #     X_test= numpy.array([5.833126006809987, 2.0091408178138224, 4.442220825763339, 1.8952175865652539, -0.13081969223378914, 7.964403170216089, 0.6197747472062328, 0.647487183032701, 0.6041092737036013, 0.6292563347957507])
    #     poet = poet_solver.POET_IC_Solver(**params)
    #     print(poet.fit_evaluate(X_test=X_test))
    #     # Easy as one two three
    #     params['type'] = '2d_period'
    #     params['features'].append(True)
    #     #X_test = numpy.append(X_test, 0.44)
    #     X_test = numpy.append(X_test, 0.43438837036704603)
    #     poet = poet_solver.POET_IC_Solver(**params)
    #     print(poet.fit_evaluate(X_test=X_test))
    #     # As simple as do re mi
    #     params['type'] = '2d_eccentricity'
    #     poet = poet_solver.POET_IC_Solver(**params)
    #     print(poet.fit_evaluate(X_test=X_test))

    params['path_to_store'] = '/home/vortebo/ctime/ai_test'
    params['version'] = '3348093'
    params['features'] = [True, True, True, True, True, True, True, True, True, True]
    X_test= numpy.array([5.833126006809987, 2.0091408178138224, 4.442220825763339, 1.8952175865652539, -0.13081969223378914, 7.964403170216089, 0.6197747472062328, 0.647487183032701, 0.6041092737036013, 0.6292563347957507])
    params['type'] = '2d_period'
    params['features'].append(True)
    X_test = numpy.append(X_test, 0.43438837036704603)
    poet = poet_solver.POET_IC_Solver(**params)
    print(poet.fit_evaluate(X_test=X_test))
    params['type'] = '2d_eccentricity'
    poet = poet_solver.POET_IC_Solver(**params)
    print(poet.fit_evaluate(X_test=X_test))

    print('Training montage complete.')

def test_numberz():
    path = '/home/vortebo/ctime/ayeye'
    params = {
                "type": 'blank',
                "epochs": 30,
                "batch_size": 100,
                "verbose": 2,
                "retrain": False,
                "threshold": 20,
                "path_to_store": path,
                "version": None,
                "features": None
            }
    systems = ['12356914','3348093']
    print('lez go')

    for system in systems:
        params['version'] = system
        # Abc
        params['type'] = '1d_period'
        params['features'] = [True, True, True, True, True, True, True, True, True, True]
        X_test= numpy.array([5.0, 10.0, 0.5, 2.0, 0.06, 69.0, 1.0, 0.7, 1.0, 0.8])
        poet = poet_solver.POET_IC_Solver(**params)
        print(poet.data_length())
        # Easy as one two three
        params['type'] = '2d_period'
        params['features'].append(True)
        X_test = numpy.append(X_test, 0.44)
        poet = poet_solver.POET_IC_Solver(**params)
        print(poet.data_length())
        # As simple as do re mi
        params['type'] = '2d_eccentricity'
        poet = poet_solver.POET_IC_Solver(**params)
        print(poet.data_length())

    print('Training montage complete.')

def f_for_i(e_i):

    disk_dissipation_age = 0.01 * units.Gyr
    lgQ_inertial_boost = 0.0
    lgQ_inertial_sharpness = 10.0
    primary_disk_lock_period = 5.0
    primary_wind_strength = 0.17
    primary_wind_saturation = 2.45
    primary_core_envelope_coupling_timescale = 0.01# * units.Gyr
    secondary_disk_lock_period = None
    secondary_wind_strength = 0.17
    secondary_wind_saturation = 2.45
    secondary_core_envelope_coupling_timescale = 0.01# * units.Gyr

    # Variables
    #initial_porb = 21.72109769815944 * units.day
    #lgQ_min = 5.1132194465309215
    #lgQ_break_period = 6.314887627048559 * units.day
    #lgQ_powerlaw = -4.336916616556421
    #initial_eccentricity = 'solve'#e_i
    #age = 2.1739862978291256
    #feh = -0.21087014875588855
    #orbital_period = 7.9644043522307575
    #primary_mass = 0.626606245275752
    #secondary_mass = 0.6640150390960553
    #cmd_primary_radius = 0.5928346964698288
    #cmd_secondary_radius = 0.6372604662783549
    ##
    # #12356914
    # lgQ_min = 6.125745092854476
    # lgQ_break_period = 10.203311339280315 * units.day
    # lgQ_powerlaw = 0.1615887253615469
    # initial_eccentricity = 'solve'#e_i
    # age = 3.072206810727815
    # feh = -0.3586404143328353
    # orbital_period = 27.308229700438552
    # primary_mass = 1.0026067368598708
    # secondary_mass = 0.6143792971318769
    # cmd_primary_radius = 1.0293764108490684
    # cmd_secondary_radius = 0.5967444520267795
    #3348093
    # lgQ_min = 6.448394485856541
    # lgQ_break_period = 2.5389791409236775 * units.day
    # lgQ_powerlaw = -0.2537823511284367
    # initial_eccentricity = 'solve'#e_i
    # age = 2.7948625956361774
    # feh = -0.11905706044362162
    # orbital_period = 7.964403037488727
    # primary_mass = 0.6280276003650294
    # secondary_mass = 0.6559427536896645
    # cmd_primary_radius = 0.6081552567189626
    # cmd_secondary_radius = 0.618546680625624
    #11403216
    # lgQ_min = 6.730834951374322
    # lgQ_break_period = 0.9524689224317893 * units.day
    # lgQ_powerlaw = 0.6928867615956023
    # initial_eccentricity = 'solve'#e_i
    # age = 4.992647080579888
    # feh = 0.16115155325450198
    # orbital_period = 4.0532537066910015
    # primary_mass = 1.1454357634717733
    # secondary_mass = 0.47249084407776804
    # cmd_primary_radius = 1.3876509776305228
    # cmd_secondary_radius = 0.4639166819029762
    #9881258
    # lgQ_min = 7.8691057445674994
    # lgQ_break_period = 0.8961115868785646 * units.day
    # lgQ_powerlaw = -2.0557821441166246
    # initial_eccentricity = e_i
    # age = 5.488213741069469
    # feh = -0.14419623497766104
    # orbital_period = 4.05702435351295
    # primary_mass = 1.1540213398413717
    # secondary_mass = 1.0038498319940319
    # cmd_primary_radius = 1.5865313989714425
    # cmd_secondary_radius = 1.3231944551036534
    #3834364
    lgQ_min = 6.5426211182269665
    lgQ_break_period = 0.7268577941281817 * units.day
    lgQ_powerlaw = 3.5804698651669025
    initial_eccentricity = 'solve'#e_i
    age = 1.270129528874416
    feh = 0.050033604581647
    orbital_period = 2.9084545603603127
    primary_mass = 1.092516035842445
    secondary_mass = 0.5244726999066902
    cmd_primary_radius = 1.081037432242053
    cmd_secondary_radius = 0.5522193325893356

    carepackage = dict()
    carepackage['lgQ_min']=lgQ_min
    carepackage['lgQ_break_period']=lgQ_break_period
    carepackage['lgQ_powerlaw']=lgQ_powerlaw
    carepackage['system_name']='12356914'#'9881258'
    carepackage['lock']=thelock
    carepackage['path']= '/home/vortebo/ctime/ayeye'

    solve = True

    # From parameters
    disk_period = primary_disk_lock_period * units.day
    secondary_is_star = True

    system = SimpleNamespace(
                        orbital_period= orbital_period * units.day,
                        age= age * units.Gyr,
                        eccentricity= e_i,#0.2561160543396648,#eccentricity,
                        primary_mass= primary_mass * units.solMass,
                        secondary_mass= secondary_mass * units.solMass,
                        feh= feh,
                        Rprimary= cmd_primary_radius * units.solRad,
                        Rsecondary= cmd_secondary_radius * units.solRad
        )

    star_dissipation = get_dissipation(lgQ_min,lgQ_inertial_boost,lgQ_inertial_sharpness,lgQ_break_period,lgQ_powerlaw)
    dissipation = dict(
                        primary = star_dissipation,
                        secondary = star_dissipation
                        )
    print('time for interpolatorer')
    interpolator = StellarEvolutionManager('/home/vortebo/ctime/poet/stellar_evolution_interpolators').get_interpolator_by_name(
            'default'
        )
    print('interpolatorer overed')
    print(interpolator)

    print('loading')
    orbital_evolution_library.prepare_eccentricity_expansion(
        '/home/vortebo/eccentricity_expansion_coef_O400.sqlite'.encode('ascii'),
        1e-4,
        True,
        True
    )
    print('ready')

    #Call
    print('run')
    result = rs.find_evolution(system = system,
        interpolator = interpolator,
        dissipation = dissipation,
        #initial_porb = initial_porb,
        initial_eccentricity = initial_eccentricity,
        #initial_obliquity = 3.0,
        disk_period=disk_period,
        disk_dissipation_age=disk_dissipation_age,
        primary_wind_strength=primary_wind_strength,
        primary_wind_saturation=primary_wind_saturation,
        primary_core_envelope_coupling_timescale=primary_core_envelope_coupling_timescale * units.Gyr,
        secondary_wind_strength=secondary_wind_strength,
        secondary_wind_saturation=secondary_wind_saturation,
        secondary_core_envelope_coupling_timescale=secondary_core_envelope_coupling_timescale * units.Gyr,
        secondary_disk_period=secondary_disk_lock_period,
        #orbital_period_tolerance=1e-6,
        #eccentricity_tolerance=1e-6,
        #obliquity_tolerance=1e-6,
        #period_search_factor=2.0,
        #scaled_period_guess=1.0,
        #eccentricity_upper_limit=0.8,
        solve=solve,
        #max_iterations=49,
        secondary_is_star=secondary_is_star,
        carepackage = carepackage,#None,
        precision = 1e-5,
        eccentricity_expansion_fname = '/home/vortebo/eccentricity_expansion_coef_O400.sqlite'.encode('ascii')
    )
    print('THIS IS IMPORTANT. THE EHAT PRIME IS AND ECC_I ARE ',result[1])
    return result[1][1]#result[0].eccentricity[-1]

def plot_ef_v_ei(systemname):

    #a=0.7009632432260313
    #input = numpy.linspace(a-.05,a+.05,11)
    #input = numpy.linspace(0.1,0.70,11)#numpy.linspace(0.1,0.64,10)#[0.141]#numpy.linspace(0.1,0.51,11)
    input = numpy.linspace(0.05470077796297116,0.1680700293901145,11)
    output = []

    from multiprocessing import Pool

    with Pool(
                4,
                initializer=setup_process_map,
                initargs=[
                    dict(
                        fname_datetime_format='%Y%m%d%H%M%S',
                        system=systemname,
                        std_out_err_fname='josh_output_3/{task}/{system}_{now}_{pid:d}.outerr',
                        logging_fname='josh_output_3/{task}/{system}_{now}_{pid:d}.log',
                        logging_verbosity='debug',
                        logging_message_format='%(levelname)s %(asctime)s %(name)s: %(message)s | %(pathname)s.%(funcName)s:%(lineno)d'
                        #logging_datetime_format=config.logging_datetime_format
                    )
                ]
            ) as p:
        output = p.map(f_for_i,input)
    
    plt.plot(output,input)
    plt.title('Final eccentricity vs initial eccentricity')
    plt.xlabel('Initial eccentricity')
    plt.ylabel('Final eccentricity')
    plt.savefig('ef_v_ei_3348093.png')
    #print(f_for_i(a))
    print('ef: ',input)
    print('ei: ',output)


def test_data_struct():

    expandabandband = []

    def subfunc():
        expandabandband.append(1)
    
    for i in range(10):
        subfunc()
    
    print(expandabandband)

def command_line_ml_test(lgQ_min = 6.730834951374322,
                         lgQ_break_period = 0.9524689224317893 * units.day,
                         lgQ_powerlaw = 0.6928867615956023,
                         age = 4.992647080579888,
                         feh = 0.16115155325450198,
                         orbital_period = 4.0532537066910015,
                         primary_mass = 1.1454357634717733,
                         secondary_mass = 0.47249084407776804,
                         cmd_primary_radius = 1.3876509776305228,
                         cmd_secondary_radius = 0.4639166819029762,
                         is_1d = True,
                         final_eccentricity = 0.2561160543396648,
                         initial_porb = 10.0 * units.day,
                         initial_eccentricity = 0.8):

    disk_dissipation_age = 0.01 * units.Gyr
    lgQ_inertial_boost = 0.0
    lgQ_inertial_sharpness = 10.0
    primary_disk_lock_period = 5.0
    primary_wind_strength = 0.17
    primary_wind_saturation = 2.45
    primary_core_envelope_coupling_timescale = 0.01# * units.Gyr
    secondary_disk_lock_period = None
    secondary_wind_strength = 0.17
    secondary_wind_saturation = 2.45
    secondary_core_envelope_coupling_timescale = 0.01# * units.Gyr

    if is_1d:
        initial_eccentricity = 0.8
    #else:
    #    initial_eccentricity = 'solve'

    solve = False#True

    # From parameters
    disk_period = primary_disk_lock_period * units.day
    secondary_is_star = True
    system = SimpleNamespace(
                        orbital_period= orbital_period * units.day,
                        age= age * units.Gyr,
                        eccentricity= final_eccentricity,
                        primary_mass= primary_mass * units.solMass,
                        secondary_mass= secondary_mass * units.solMass,
                        feh= feh,
                        Rprimary= cmd_primary_radius * units.solRad,
                        Rsecondary= cmd_secondary_radius * units.solRad
        )
    star_dissipation = get_dissipation(lgQ_min,lgQ_inertial_boost,lgQ_inertial_sharpness,lgQ_break_period,lgQ_powerlaw)
    dissipation = dict(
                        primary = star_dissipation,
                        secondary = star_dissipation
                        )
    interpolator = StellarEvolutionManager('/home/vortebo/ctime/poet/stellar_evolution_interpolators').get_interpolator_by_name(
            'default'
        )
    orbital_evolution_library.prepare_eccentricity_expansion(
        '/home/vortebo/eccentricity_expansion_coef_O400.sqlite'.encode('ascii'),
        1e-4,
        True,
        True
    )

    #Call
    result = rs.find_evolution(system = system,
        interpolator = interpolator,
        dissipation = dissipation,
        initial_porb = initial_porb,
        initial_eccentricity = initial_eccentricity,
        disk_period=disk_period,
        disk_dissipation_age=disk_dissipation_age,
        primary_wind_strength=primary_wind_strength,
        primary_wind_saturation=primary_wind_saturation,
        primary_core_envelope_coupling_timescale=primary_core_envelope_coupling_timescale * units.Gyr,
        secondary_wind_strength=secondary_wind_strength,
        secondary_wind_saturation=secondary_wind_saturation,
        secondary_core_envelope_coupling_timescale=secondary_core_envelope_coupling_timescale * units.Gyr,
        secondary_disk_period=secondary_disk_lock_period,
        solve=solve,
        secondary_is_star=secondary_is_star,
        carepackage = None,
        precision = 1e-5,
        eccentricity_expansion_fname = '/home/vortebo/eccentricity_expansion_coef_O400.sqlite'.encode('ascii')
    )
    print('Initial period: ', result[0].orbital_period[0])
    print('Initial eccentricity: ', result[0].eccentricity[0])
    print('Final period: ', result[0].orbital_period[-1])
    print('Final eccentricity: ', result[0].eccentricity[-1])
    return 0

def whystardissipationcrashwtf(lgQ_min = 6.730834951374322,
                         lgQ_break_period = 0.9524689224317893 * units.day,
                         lgQ_powerlaw = 0.6928867615956023,
                         age = 4.992647080579888,
                         feh = 0.16115155325450198,
                         orbital_period = 4.0532537066910015,
                         primary_mass = 1.1454357634717733,
                         secondary_mass = 0.47249084407776804,
                         cmd_primary_radius = 1.3876509776305228,
                         cmd_secondary_radius = 0.4639166819029762,
                         is_1d = True,
                         final_eccentricity = 0.2561160543396648,
                         initial_porb = 10.0 * units.day,
                         initial_eccentricity = 0.8,
                         angmom_a = 0.0,
                         angmom_b = 0.0):
    
    logger=logging.getLogger(__name__)

    disk_dissipation_age = 0.01 * units.Gyr
    lgQ_inertial_boost = 0.0
    lgQ_inertial_sharpness = 10.0
    primary_disk_lock_period = 5.0
    primary_wind_strength = 0.17
    primary_wind_saturation = 2.45
    primary_core_envelope_coupling_timescale = 0.01 * units.Gyr
    secondary_disk_lock_period = None
    secondary_wind_strength = 0.17
    secondary_wind_saturation = 2.45
    secondary_core_envelope_coupling_timescale = 0.01 * units.Gyr
    
    disk_period = primary_disk_lock_period * units.day
    secondary_is_star = True
    system = SimpleNamespace(
                        orbital_period= orbital_period * units.day,
                        age= age * units.Gyr,
                        eccentricity= final_eccentricity,
                        primary_mass= primary_mass * units.solMass,
                        secondary_mass= secondary_mass * units.solMass,
                        feh= feh,
                        Rprimary= cmd_primary_radius * units.solRad,
                        Rsecondary= cmd_secondary_radius * units.solRad
        )
    star_dissipation = get_dissipation(lgQ_min,lgQ_inertial_boost,lgQ_inertial_sharpness,lgQ_break_period,lgQ_powerlaw)
    #star_dissipation['spin_frequency_breaks'] = numpy.array([None])
    dissipation = dict(
                        primary = star_dissipation,
                        secondary = star_dissipation
                        )
    MESAInterpolator.set_quantity_lower_limit('iconv', 1)#e-5)
    interpolator = StellarEvolutionManager('/home/vortebo/ctime/poet/stellar_evolution_interpolators').get_interpolator_by_name(
            'default'
        )

    logger.debug('Now loading: eccentricity expansion coefficients.')
    import os
    orbital_evolution_library.prepare_eccentricity_expansion(
        os.path.expanduser(
            #'~/projects/git/poet/eccentricity_expansion_coef_O400.sqlite'
            '~/eccentricity_expansion_coef_O400.sqlite'
        ).encode('ascii'),
        #config.eccentricity_expansion_coefficients.encode('ascii'),
        #'/home/vortebo/eccentricity_expansion_coef_O400.sqlite',
        1e-4,
        True,
        True
    )
    logger.debug('Eccentricity expansion coefficients loaded.')
    
    secondary_disk_period=None
    orbital_period_tolerance=1e-6
    period_search_factor=2.0
    scaled_period_guess=1.0
    extra_evolve_args = dict()
    extra_evolve_args['max_time_step'] = 1e-3
    extra_evolve_args['precision'] = 1e-5

    from general_purpose_python_modules.solve_for_initial_values import \
        InitialValueFinder

    value_finder = InitialValueFinder(
        system=system,
        interpolator=interpolator,
        #False positive
        #pylint: disable=no-member
        current_age=system.age,
        #pylint: enable=no-member
        disk_period=disk_period,
        initial_obliquity=0.0,
        disk_dissipation_age=disk_dissipation_age,
        dissipation=dissipation,
        secondary_star=secondary_is_star,
        primary_wind_strength=primary_wind_strength,
        primary_wind_saturation=primary_wind_saturation,
        primary_core_envelope_coupling_timescale=(
            primary_core_envelope_coupling_timescale
        ),
        secondary_wind_strength=secondary_wind_strength,
        secondary_wind_saturation=secondary_wind_saturation,
        secondary_core_envelope_coupling_timescale=(
            secondary_core_envelope_coupling_timescale
        ),
        secondary_disk_period=secondary_disk_period,
        orbital_period_tolerance=orbital_period_tolerance,
        period_search_factor=period_search_factor,
        scaled_period_guess=scaled_period_guess,
        **extra_evolve_args
    )
    logger.debug('Okay guys let\'s do this LEEEEROOOOOY JENKINS')
    initial_secondary_angmom = [angmom_a,angmom_b]#value_finder.get_secondary_initial_angmom()
    logger.debug('At least I have chicken')

    evolution = value_finder.try_system([initial_porb,initial_eccentricity,0.0],initial_secondary_angmom)
    logger.debug('I cry')
    #porb_found, ecc_found = evolution.orbital_period[-1], evolution.eccentricity[-1]
    plt.plot(evolution.age,evolution.primary_iconv,label='Primary')
    #plt.scatter(evolution.age,evolution.secondary_iconv,label='Secondary')
    plt.savefig('/home/vortebo/ctime/zzzzzz_iconv.png')

if __name__ == '__main__':
    print("starting")
    systemname='6521542'
    setup_process(
                    fname_datetime_format='%Y%m%d%H%M%S',
                    system=systemname,
                    std_out_err_fname='josh_output_9/{task}/{system}_{now}_{pid:d}.outerr',
                    logging_fname='josh_output_9/{task}/{system}_{now}_{pid:d}.log',
                    logging_verbosity='debug',
                    logging_message_format='%(levelname)s %(asctime)s %(name)s: %(message)s | %(pathname)s.%(funcName)s:%(lineno)d'#,
                    #logging_datetime_format=config.logging_datetime_format
                  )

    # whystardissipationcrashwtf(
    #     #4947726 (not fixed by iconv)
    #     # lgQ_min = 7.185163961316851,
    #     # lgQ_break_period = 1.2105510858211672 * units.day,
    #     # lgQ_powerlaw = -1.8675025503172198,
    #     # age = 4.537071400653452,
    #     # feh = -0.1476147128380225,
    #     # orbital_period = 4.726086038265389,
    #     # primary_mass = 1.0349127223165384,
    #     # secondary_mass = 0.5186184463883605,
    #     # cmd_primary_radius = 1.114127831417474,
    #     # cmd_secondary_radius = 0.5075721251761187,
    #     # is_1d = True,
    #     # final_eccentricity = 0.0,
    #     # initial_porb = 22.834226476460383,
    #     # initial_eccentricity = 0.8,
    #     # angmom_a = 0.10428821085789115,
    #     # angmom_b = 0.00019044827608288606
    #     #6521542 (yes fixed? :O)
    #     lgQ_min = 9.222211210635884,
    #     lgQ_break_period = 0.8597019976850838 * units.day,
    #     lgQ_powerlaw = 4.012989378277558,
    #     age = 1.5789897540930298,
    #     feh = -0.8460936639710941,
    #     orbital_period = 4.425754132228538,
    #     primary_mass = 1.0861970436931017,
    #     secondary_mass = 1.043207183779346,
    #     cmd_primary_radius = 1.0887086343631127,
    #     cmd_secondary_radius = 1.043463796711456,
    #     is_1d = True,
    #     final_eccentricity = 0.0,
    #     initial_porb = 13.277262396685614,
    #     initial_eccentricity = 0.8,
    #     angmom_a = 0.020857808807702508,
    #     angmom_b = 0.40876059990796304
    # )
    #retry_system()
    #read_fits()
    #emergency_plotter()
    #test_age_split()
    trainnn()
    #test_numberz()
    #plot_ef_v_ei(systemname)
    #print(f_for_i(0.35))#0.8))

        # print('Running ml test')
        # # lgQ_min = 5.05931329805408
        # # lgQ_break_period = 36.2321034083127 * units.day
        # # lgQ_powerlaw = -0.184784404574091
        # # agee = 0.344586723279352
        # # feh = 0.0740881169972182
        # # orbital_period = 8.46440234619176
        # # primary_mass = 0.654553462304074
        # # secondary_mass = 0.684806628768325
        # # cmd_primary_radius = 0.603691378563715
        # # cmd_secondary_radius = 0.621470624061517
        # final_eccentricity = 0.2561160543396648
        # # initial_porb = 11.0753937717073 * units.day
        # # initial_eccentricity = 0.5
        # lgQ_min	=	11.70853731		
        # lgQ_break_period	=	1.514387516	*	units.day
        # lgQ_powerlaw	=	-3.111613814		
        # agee	=	1.210009354		
        # feh	=	-0.076536918		
        # orbital_period	=	23.89312732		
        # primary_mass	=	0.641352039		
        # secondary_mass	=	0.669182491		
        # cmd_primary_radius	=	0.593650971		
        # cmd_secondary_radius	=	0.634637896		
        # initial_porb	=	23.89320909	*	units.day
        # initial_eccentricity	=	0.5		

        # print('1d')
        # #command_line_ml_test(lgQ_min,lgQ_break_period,lgQ_powerlaw,agee,feh,orbital_period,primary_mass,secondary_mass,cmd_primary_radius,cmd_secondary_radius,True,final_eccentricity)
        # #print('2d')
        # #final_eccentricity = 0.320081949219306
        # command_line_ml_test(lgQ_min,lgQ_break_period,lgQ_powerlaw,agee,feh,orbital_period,primary_mass,secondary_mass,cmd_primary_radius,cmd_secondary_radius,True,final_eccentricity,initial_porb,initial_eccentricity)

    print("Done!")