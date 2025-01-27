import math
#import general_purpose_python_modules.reproduce_system as rs
from types import SimpleNamespace
from astropy import units
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

import scipy
from general_purpose_python_modules.solve_for_initial_values import \
    InitialValueFinder
from josh_scripts import get_dissipation

import random
import pandas as pd
import os

import time

                #BEGIN THE CHANGES
                # Save the most recent point
                # For all points before that
                #     Discard if they are too close or too far from the most recent point
                # If no points remain
                #     Run two evolutions, one with small offset in initial period and one with small offset in initial eccentricity
                #     Use the results for the two other points
                # Else
                #     For all acceptable points before the most recent point i
                #         For all acceptable points after point i (j)
                #             If the resulting triangle area is acceptable
                #                 Use the results for the two other points
                #                 Break out of the loop
                #     If we don't have two identified points
                #         Grab most recent acceptable point
                #         Find another point perpendicular to the line between the most recent acceptable point and the most recent point (using initial values)
                #         Run an evolution for that point
                #         Use the results for the two other points
                # Proceed with existing least square stuff
                #END THE CHANGES
                # For debug, I should just have the algorithm available in a separate script and this part should save the points used when it gets here

def just_make_triangles():
    def triangle_area(A, B):
        return (1/2) * (A[0]*(B[1] - B[2] ) + A[1]*(B[2] - B[0] ) + A[2]*(B[0] - B[1]))
    # Not part of logic, gotta define some things
    porb_i = 4.809120437793041
    ecc_i = 0.04046079279696418
    obliq_i = 0.0
    thetype = '2d'
    ecc_found = 0.004150996553256088
    porb_found = 4.803506701241148
    #f
    laststeps = [[0.5, 0.273460024369773, 7.7269703200100786, 11.24430035579992], [0.5, 0.273460024369773, 7.7269703200100786, 11.24430035579992], [0.5, 0.273460024369773, 7.7269703200100786, 11.24430035579992], [0.5, 0.27346002406468156, 7.726970314942611, 11.244300352025645], [0.5000000074505806, 0.27346002699381683, 7.726970325039979, 11.244300523353049], [0.28153736625739895, 0.16756969838466734, 7.306108907458704, 8.1341801060198], [0.04494260170457748, 0.028454119863657864, 7.143764242451655, 7.232194072502805], [0.038789505211379166, 0.025168909322323543, 7.211088906571491, 7.292258679494884], [0.02477744918044396, 0.01640431204831755, 7.265000378803102, 7.337068347324254], [0.014284490031355707, 0.009578959447735914, 7.306060873347519, 7.372844875129992], [0.01471757346200407, 0.00986456610390006, 7.304416087396882, 7.371383301752334]]
    
    print('Beginning logic...')

    delta_p = 0.01*porb_found   #TODO: proper logger statements
    alpha = 0.5
    tri_max = 0.25
    tri_min = 1e-8
    max_dist = 0.1
    min_dist = 0.001
    second_point = []
    third_point = []
    acceptable_points = []
    print('There are ',len(laststeps),' point(s) in laststeps.')
    # For all points before the most recent point
    for i in range(len(laststeps)):
        # Get the distance between the most recent point and the current point in ei,pf space
        distance = numpy.sqrt((alpha*(ecc_i-laststeps[i][0]))**2 + (porb_found-laststeps[i][2])**2)
        print('Distance between ',laststeps[i],' and ',[ecc_i,ecc_found,porb_found,porb_i],' is ',distance)
        # If the point is neither too close to nor too far from the most recent point
        if max_dist > distance > min_dist:
            # Keep it
            acceptable_points.append(laststeps[i])
    print('Removed ',len(laststeps)-len(acceptable_points),' point(s).')
    # For all acceptable points before the most recent point i
    print('Acceptable points remain after filtering.')
    acceptable_points.reverse()
    acceptable_points = [[0.035277868166641996, 0.003465712275945747, 4.806879300851904, 4.808762509811369], [0.01302718714256744, 0.0011532766650002608, 4.820087669828977, 4.8117123761741505], [0.1016332591081317, 0.022201382384644384, 4.809475466147676, 4.891072617240518]]
    print('Acceptable points are ',acceptable_points)
    breakloop = False
    i = 0
    while not breakloop and i < len(acceptable_points):
        # For all acceptable points after point i (j)
        for j in range(i+1,len(acceptable_points)):
            # If the resulting triangle area (in ef vs ei) is acceptable
            print('Checking triangle area for ',acceptable_points[i],' and ',acceptable_points[j])
            tri_area = triangle_area([ecc_i,acceptable_points[i][0],acceptable_points[j][0]],[porb_found,acceptable_points[i][2],acceptable_points[j][2]])
            print('Triangle area is ',tri_area)

            A = [
                    [ecc_i,porb_found,1],
                    [acceptable_points[i][0],acceptable_points[i][2],1],
                    [acceptable_points[j][0],acceptable_points[j][2],1]
                ]
            B = [ecc_found,acceptable_points[i][1],acceptable_points[j][1]]
            A = numpy.matrix(A)
            B = numpy.matrix(B)
            fit = scipy.linalg.lstsq(A,B.T)[0][0]
            print('fit is ',fit)
        i += 1
    print('The end')

def new_logic():
    def evolve_system(initial_conditions,initial_secondary_angmom,thetype):
        print('evolve_system called with initial_conditions ',initial_conditions,' and initial_secondary_angmom ',initial_secondary_angmom,' and thetype ',thetype)
        try:
            evolution = value_finder.try_system(initial_conditions,initial_secondary_angmom,
                                                thetype,
                                                None)
            return evolution
        except AssertionError as err:
            print('AssertionError: ',err)
            return numpy.nan
        except:
            print('Something unknown went wrong while trying to evolve the system with the given parameters: ',repr(initial_conditions))
            raise
    
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
    # Variables
    lgQ_min = 6.448394485856541
    lgQ_break_period = 2.5389791409236775 * units.day
    lgQ_powerlaw = -0.2537823511284367
    initial_eccentricity = 'solve'#e_i
    age = 2.7948625956361774
    feh = -0.11905706044362162
    orbital_period = 7.964403037488727
    primary_mass = 0.6280276003650294
    secondary_mass = 0.6559427536896645
    cmd_primary_radius = 0.6081552567189626
    cmd_secondary_radius = 0.618546680625624
    solve = True
    disk_period = primary_disk_lock_period * units.day
    secondary_is_star = True
    system = SimpleNamespace(
                        orbital_period= orbital_period * units.day,
                        age= age * units.Gyr,
                        eccentricity= 0.2561288517812237,#e_i,
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
    print('Loading interpolator and preparing eccentricity expansion...')
    interpolator = StellarEvolutionManager('/home/vortebo/ctime/poet/stellar_evolution_interpolators').get_interpolator_by_name(
            'default'
        )
    orbital_evolution_library.prepare_eccentricity_expansion(
        '/home/vortebo/eccentricity_expansion_coef_O400.sqlite'.encode('ascii'),
        1e-4,
        True,
        True
    )
    extra_evolve_args = {}
    extra_evolve_args['max_time_step'] = 1e-3
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
        secondary_disk_period=secondary_disk_lock_period,
        orbital_period_tolerance=1e-6,
        period_search_factor=2.0,
        scaled_period_guess=1.0,
        precision = 1e-5,
        **extra_evolve_args
    )
    initial_secondary_angmom = numpy.array(value_finder.get_secondary_initial_angmom())

    # Not part of logic, gotta define some things
    porb_i = 16.478945137893465
    ecc_i = 0.6890106969773595
    obliq_i = 0.0
    thetype = '2d'
    ecc_found = 0.46899999872998926
    porb_found = 7.964403007017858
    #f
    laststeps = [[0.5, 0.18402829029713652, 5.11894776012314, 8.447648006477143], [0.5, 0.18402829029713652, 5.11894776012314, 8.447648006477143], [0.5, 0.18402829029713652, 5.11894776012314, 8.447648006477143], [0.5, 0.18402828961776568, 5.118947754327968, 8.44764800140039], [0.5000000074505806, 0.18402829303840673, 5.1189477495335565, 8.447648132356907], [0.4090775964701719, 0.3274266552010216, 8.202596957929908, 9.452637911808354], [0.7499378991591802, 0.6536135373606679, 15.785583429650337, 25.352197382005713], [0.4090775964701719, 0.3274266537541041, 8.202596907698585, 9.45263787553826], [0.4090776025659031, 0.32742666107698204, 8.20259701666803, 9.452637996734301], [0.49652634485781216, 0.3572773085811191, 7.572309296348227, 10.10395147571947], [0.4574269156351555, 0.35155984228616377, 8.022737408473747, 9.829214938774337], [0.5687418062586453, 0.397586914367288, 7.623964932930915, 11.627709087706807], [0.5187658156301823, 0.3864101479255369, 8.011521190777056, 10.773062463415382], [0.640607403292205, 0.44365744384449196, 7.880654218207265, 14.076925366067773], [0.6841660327805495, 0.4689658160322701, 8.019649453189164, 16.233778670392336], [0.6888950502376873, 0.46692361964618406, 7.913098359301637, 16.44221895808077], [0.6890497659062806, 0.46955758459635355, 7.978200380179627, 16.48930193775666], [0.6890184090190349, 0.46900390280343457, 7.964414439553928, 16.47940304339212], [0.6890107179626889, 0.46899813597356527, 7.96435526077015, 16.478918445884684], [0.689010704078439, 0.4690000531898011, 7.964404314718983, 16.478946317815012]]
    
    print('Beginning logic...')

    delta_p = 0.01*porb_i#porb_found
    alpha = 0.5
    tri_max = 1e-5
    tri_min = 1e-7
    max_dist = 0.001
    min_dist = 0.0001
    second_point = []
    third_point = []
    acceptable_points = []
    print('There are ',len(laststeps),' point(s) in laststeps.')
    # For all points before the most recent point
    for i in range(len(laststeps)):
        # Get the distance between the most recent point and the current point
        distance = numpy.sqrt((alpha*(ecc_found-laststeps[i][1]))**2 + (porb_found-laststeps[i][2])**2) #TODO: this isn't ei,pf
        print('Distance between ',laststeps[i],' and ',[ecc_i,ecc_found,porb_found,porb_i],' is ',distance)
        # If the point is neither too close to nor too far from the most recent point
        if max_dist > distance > min_dist:
            # Keep it
            acceptable_points.append(laststeps[i])
    print('Removed ',len(laststeps)-len(acceptable_points),' point(s).')
    # If no points remain after filtering
    if len(acceptable_points) == 0:
        # Run two evolutions, one with small offset in initial period and one with small offset in initial eccentricity
        # Use the results for the two other points
        print('No acceptable points remain after filtering.')
        point_a_period = porb_i + delta_p
        point_b_ecc = ecc_i + (alpha*delta_p)
        point_a_i = [point_a_period,ecc_i,obliq_i]
        print('point_a_i is ',point_a_i)
        evolution = evolve_system(point_a_i,initial_secondary_angmom,'2d')
        if not hasattr(evolution, 'eccentricity'):
            raise AssertionError("Something went wrong trying to find a second point.",0)
        second_point = [ecc_i,evolution.eccentricity[-1],evolution.orbital_period[-1],point_a_period]
        print('second_point is ',second_point)
        point_b_i = [porb_i,point_b_ecc,obliq_i]
        print('point_b_i is ',point_b_i)
        evolution = evolve_system(point_b_i,initial_secondary_angmom,'2d')
        if not hasattr(evolution, 'eccentricity'):
            raise AssertionError("Something went wrong trying to find a third point.",0)
        third_point = [point_b_ecc,evolution.eccentricity[-1],evolution.orbital_period[-1],porb_i]
        print('third_point is ',third_point)
    # Otherwise, points did remain, so let's work with that
    else:
        # For all acceptable points before the most recent point i
        print('Acceptable points remain after filtering.')
        acceptable_points.reverse()
        print('Acceptable points are ',acceptable_points)
        breakloop = False
        i = 0
        while not breakloop and i < len(acceptable_points):
            # For all acceptable points after point i (j)
            for j in range(i+1,len(acceptable_points)):
                # If the resulting triangle area (in ef vs ei) is acceptable
                print('Checking triangle area for ',acceptable_points[i],' and ',acceptable_points[j])
                print('Triangle area is ',(1/2) * (
                    ecc_i*(acceptable_points[i][1] - acceptable_points[j][1] )
                    +
                    acceptable_points[i][0]*(acceptable_points[j][1] - ecc_found )
                    +
                    acceptable_points[j][0]*(ecc_found - acceptable_points[i][1])
                    ))
                if (tri_max > numpy.abs((1/2) * (
                    ecc_i*(acceptable_points[i][1] - acceptable_points[j][1] )
                    +
                    acceptable_points[i][0]*(acceptable_points[j][1] - ecc_found )
                    +
                    acceptable_points[j][0]*(ecc_found - acceptable_points[i][1])
                    )) > tri_min
                ):
                    # Use the results for the two other points
                    print('Triangle area is acceptable.')
                    second_point = acceptable_points[i]
                    print('second_point is ',second_point)
                    third_point = acceptable_points[j]
                    print('third_point is ',third_point)
                    # Break out of the loop
                    breakloop = True
                    break
            i += 1
        # If we don't have two identified points
        print('second_point is ',second_point)
        print('third_point is ',third_point)
        if len(second_point) == 0 or len(third_point) == 0:
            # We're going to be working in ei,pi phase space instead of ei,pf #TODO: make sure you're doing ei,pf elsewhere
            print('Unable to find two points that avoid a degenerate solution.')
            # Grab most recent acceptable point
            second_point = acceptable_points[0]
            print('second_point is ',second_point)
            # Find another point perpendicular to the line between the most recent acceptable point and the most recent point (using initial values)
            ei_two = second_point[0]
            pi_two = second_point[3]
            distance = numpy.sqrt((alpha*(ecc_i-ei_two))**2 + (porb_i-pi_two)**2)
            print('distance is ',distance)
            #theta = numpy.arcsin(alpha*(ecc_i-ei_two)/distance)
            #print('theta is ',theta)
            print('ei_two is ',ei_two,' and ecc_i is ',ecc_i)
            print('alpha is ',alpha,)
            print('ei_two-ecc_i is ',ei_two-ecc_i)
            print('alpha*(ei_two-ecc_i) is ',alpha*(ei_two-ecc_i))
            del_e = alpha*(ei_two-ecc_i)
            del_p = pi_two-porb_i
            print('Delta e is ',del_e)
            print('Delta p is ',del_p)
            ei_three = ecc_i + del_p / alpha
            pi_three = porb_i - del_e
            print('ei_three is ',ei_three,' and pi_three is ',pi_three)
            print('New distance is ',numpy.sqrt((alpha*(ei_three-ecc_i))**2 + (pi_three-porb_i)**2))
            # Run an evolution for that point
            # TODO: Check if this is a valid point? If I've locked it to actually be perpendicular, then it should be always be a small enough
            #       difference to be fine, right?
            point_c_i = [pi_three,ei_three,obliq_i]
            evolution = evolve_system(point_c_i,initial_secondary_angmom,'2d')
            if not hasattr(evolution, 'eccentricity'):
                raise AssertionError("Something went wrong trying to find a third point.",0)
            # Use the results for the third point
            third_point = [ei_three,evolution.eccentricity[-1],evolution.orbital_period[-1],pi_three]
            print('third_point is ',third_point)
    # Final check to make sure the triangle area is acceptable
    print('second_point is ',second_point)
    print('third_point is ',third_point)
    A = [
            [ecc_i,porb_found,1],
            [second_point[0],second_point[2],1],
            [third_point[0],third_point[2],1]
        ]
    B = [ecc_found,second_point[1],third_point[1]]
    print('A is ',A)
    print('B is ',B)
    print('Triangle area is ',(1/2) * (A[0][0]*(B[1] - B[2] ) + A[1][0]*(B[2] - B[0] ) + A[2][0]*(B[0] - B[1])))
    #if not (tri_max > numpy.abs((1/2) * (A[0][0]*(B[1] - B[2] ) + A[1][0]*(B[2] - B[0] ) + A[2][0]*(B[0] - B[1]))) > tri_min):
    #    print('Unable to find two points that avoid a degenerate solution.')
    #    raise ValueError("Unable to find two points that avoid a degenerate solution.",0)

    # Perform least squares fit to find ehat_prime
    print('Performing least squares fit to find ehat_prime...')
    A = numpy.matrix(A)
    B = numpy.matrix(B)
    fit = scipy.linalg.lstsq(A,B.T)[0]
    ehat_prime = [fit[0],ecc_i]
    print('A is ',A)
    print('B is ',B)
    print('ehat_prime is ',ehat_prime)

def check_for_bad():
    path = '/home/vortebo/ctime/ai_test'
    typez='1d_period'
    version='3348093'
    
    data_df = pd.read_csv(f"/{path}/poet_output/{typez}_{version}/datasets/data.csv.gz",
                            compression='gzip')
    labels_df = pd.read_csv(f"/{path}/poet_output/{typez}_{version}/datasets/label.csv.gz",
                            compression='gzip')
    #
    # Train using specified list of features
    #
    #if self.features:
    #    data_df = data_df.loc[:, self.features]
    #
    # Convert pandas dataframe to numpy ndarray
    #
    data_df = data_df.to_numpy()
    labels_df = labels_df.to_numpy()
    to_remove = []

    # For each value in data_df, check if it's NaN
    for i in range(len(data_df)):
        #print(data_df[i])
        if numpy.isnan(data_df[i]).any():
            print('NaN found in data_df at index ',i)
            print('data_df[i] is ',data_df[i])
            print('labels_df[i] is ',labels_df[i])
            # Remember the index to remove
            to_remove.append(i)
    for i in range(len(labels_df)):
        #print(labels_df[i])
        if numpy.isnan(labels_df[i]).any():
            print('NaN found in labels_df at index ',i)
            print('data_df[i] is ',data_df[i])
            print('labels_df[i] is ',labels_df[i])
            # Remember the index to remove
            to_remove.append(i)
    
    # Remove the indices
    data_df = numpy.delete(data_df,to_remove,0)
    labels_df = numpy.delete(labels_df,to_remove,0)
    
    # Save the dataframes
    data_df = pd.DataFrame(data_df)
    labels_df = pd.DataFrame(labels_df)
    data_df.to_csv(f"/{path}/poet_output/{typez}_{version}/datasets/data.csv.gz",
                    compression='gzip',
                    index=False)
    labels_df.to_csv(f"/{path}/poet_output/{typez}_{version}/datasets/label.csv.gz",
                    compression='gzip',
                    index=False)

def join_data():
    path = '/home/vortebo/ctime/ai_test3'
    typez=['1d_period','2d_period','2d_eccentricity']
    #version=['3348093','10031409','10091257','11403216','12356914']
    version=[
                '5039441','6312521','6962018','7732791','7798259','9110346','9119652','9656543',
                '10129482','6697716','11234677','7970629','4908495','10091257','6525196','5288543',
                '8196180','5022440','4276114','5802470','4815612','3241344','7377033','3427776',
                '11403216','9881258','10935310','10031409','3973504','8957954','6521542','5816806',
                '11252617','4285087','7025851','4346875','4947726','7691527','6227560','8302455',
                '12004679','10483644','7369523','9971475','8229048','7129465','5181455','8381592',
                '7376500','8618226','3439031','6546508','10385682','5359678','9715925','7125636',
                '8580438','10965963','5263802','7960547','8746310','7362852','7128918','4753988',
                '3348093','10031409','10091257','10198109','11403216','12356914','12644769'
             ]
    
    final_version = 'joint'
    
    for thetype in typez:
        data = numpy.array([])
        labels = numpy.array([])

        for ver in version:
            if not os.path.exists(f"/{path}/poet_output/{thetype}_{ver}/datasets"):
                print(f"/{path}/poet_output/{thetype}_{ver}/datasets does not exist.")
                continue
            # Load the dataframes
            try:
                data_a = pd.read_csv(f"/{path}/poet_output/{thetype}_{ver}/datasets/data.csv")
                labels_a = pd.read_csv(f"/{path}/poet_output/{thetype}_{ver}/datasets/label.csv")
            except:
                print(f"Failed to load dataframes for {thetype}_{ver}.")
                continue

            data_a = data_a.to_numpy()
            labels_a = labels_a.to_numpy()

            if data.ndim != data_a.ndim:
                data = data_a
                labels = labels_a
            else:
                data = numpy.concatenate((data,data_a),axis=0)
                labels = numpy.concatenate((labels,labels_a),axis=0)
            
        # Save the dataframes
        data = pd.DataFrame(data)
        labels = pd.DataFrame(labels)
        # Check the directory exists. Make it if it doesn't.
        if not os.path.exists(f"/{path}/poet_output/{thetype}_{final_version}/datasets"):
            os.makedirs(f"/{path}/poet_output/{thetype}_{final_version}/datasets")

        data.to_csv(f"/{path}/poet_output/{thetype}_{final_version}/datasets/data.csv",
                        index=False)
        labels.to_csv(f"/{path}/poet_output/{thetype}_{final_version}/datasets/label.csv",
                        index=False)

def make_data():
    path = '/home/vortebo/ctime/baibaibai'
    params = {
                "type": 'blank',
                "epochs": 300,
                "batch_size": 300000,
                "verbose": 2,
                "retrain": False,
                "threshold": 20,
                "path_to_store": path,
                "version": 'Panic',
                "features": None#[True]
            }
    ai_model = poet_solver.POET_IC_Solver(**params)
    
    into = random.sample(range(500000),100000)
    into = numpy.array(into)
    out = into + 5
    #for i in range(len(into)):
    #    ai_model.store_data(X_train=numpy.matrix(into[i]), y_train=numpy.matrix(out[i]))
    #ai_model.store_data(X_train=into, y_train=out)
    
    test = [64, 11.7, 23, 12]
    for i in range(len(test)):
        print(ai_model.fit_evaluate(X_test=numpy.matrix(test[i])))

def trainnn():
    path = '/home/vortebo/ctime/ayeye'
    params = {
                "type": '1d_period',
                "epochs": 300,
                "batch_size": 100,
                "verbose": 2,
                "retrain": False,
                "threshold": 20,
                "path_to_store": path,
                "version": None,
                "features": None
            }
    systems = ['12356914','3348093']

    params['path_to_store'] = '/home/vortebo/ctime/ai_test'
    params['version'] = '3348093'
    params['features'] = [True, True, True, True, True, True, True, True, True, True]
    X_test= numpy.array([5.833126006809987, 2.0091408178138224, 4.442220825763339, 1.8952175865652539, -0.13081969223378914, 7.964403170216089, 0.6197747472062328, 0.647487183032701, 0.6041092737036013, 0.6292563347957507])
    poet = poet_solver.POET_IC_Solver(**params)
    print(poet.fit_evaluate(X_test=X_test))
    params['type'] = '2d_period'
    params['features'].append(True)
    X_test = numpy.append(X_test, 0.43438837036704603)
    poet = poet_solver.POET_IC_Solver(**params)
    print(poet.fit_evaluate(X_test=X_test))
    params['type'] = '2d_eccentricity'
    poet = poet_solver.POET_IC_Solver(**params)
    print(poet.fit_evaluate(X_test=X_test))

    print('Training montage complete.')

def try_values(lgQm,lgQb,lgQp,tf,z,pf,m1,m2,r1,r2,ef,path=None,version=None):
    params = {
                "type": '1d_period',
                "epochs": 300,
                "batch_size": 100,
                "verbose": 2,
                "retrain": False,
                "threshold": 20,
                "path_to_store": '/home/vortebo/ctime/ai_test',
                "version": '3348093',
                "features": [True, True, True, True, True, True, True, True, True, True]
            }
    
    if path is not None and version is not None:
        params['path_to_store'] = path
        params['version'] = version
    params['features'] = [True, True, True, True, True, True, True, True, False, False]
    
    #data_test = pd.read_csv(f"/{path}/poet_output/{params['type']}_{version}/datasets/data.csv")
    #labels_test = pd.read_csv(f"/{path}/poet_output/{params['type']}_{version}/datasets/label.csv")
    #data_test = data_test.to_numpy()
    #labels_test = labels_test.to_numpy()
    
    X_test= numpy.array([lgQm, lgQb, lgQp, tf, z, pf, m1, m2, r1, r2])
    poet = poet_solver.POET_IC_Solver(**params)
    for i in range(10):
        a1dpdown,a1dpup = poet.fit_evaluate(X_test=X_test)#data_test[0])
    #poet.report_model_details()
    params['type'] = '2d_period'
    params['features'].append(True)
    X_test = numpy.append(X_test, ef)
    poet = poet_solver.POET_IC_Solver(**params)
    for i in range(10):
        a2dpdown,a2dpup = poet.fit_evaluate(X_test=X_test)
    #poet.report_model_details()
    params['type'] = '2d_eccentricity'
    poet = poet_solver.POET_IC_Solver(**params)
    for i in range(10):
        a2dedown,a2deup = poet.fit_evaluate(X_test=X_test)
    #poet.report_model_details()

    result_1d = 0.5*(a1dpdown+a1dpup)
    result_2dp = 0.5*(a2dpdown+a2dpup)
    result_2de = 0.5*(a2dedown+a2deup)

    print('RESULTS')
    print('Result for 1d_period is ',result_1d)
    print('Result for 2d_period is ',result_2dp)
    print('Result for 2d_eccentricity is ',result_2de)

def validate(version='3241344'):
    types = ['1d_period']#,'2d_period','2d_eccentricity']
    #types = ['2d_period','2d_eccentricity']
    #systems = ['3348093']
    #version='3241344'#'joint'#'3348093'
    path1 = '/home/vortebo/ctime/ai_test3'
    path2 = '/home/vortebo/ctime/crycry4'
    path1 = path2
    if True:
        for typez in types:
            #for version in systems:
            #typez='1d_period'
            #version='3348093'
            #
            data_df = pd.read_csv(f"/{path1}/poet_output/{typez}_{version}/datasets/data.csv")
            labels_df = pd.read_csv(f"/{path1}/poet_output/{typez}_{version}/datasets/label.csv")
            #
            data_df = data_df.to_numpy()
            labels_df = labels_df.to_numpy()
            # 
            TRAIN = 0.85
            TEST = 0.3
            # Split the data into training and testing sets by pulling out a random 70% of the data for training
            train_indices = random.sample(range(len(data_df)),int(TRAIN*len(data_df)))
            test_indices = [i for i in range(len(data_df)) if i not in train_indices]
            data_train = data_df[train_indices]
            data_test = data_df[test_indices]
            labels_train = labels_df[train_indices]
            labels_test = labels_df[test_indices]
            # Save the dataframes
            data_train = pd.DataFrame(data_train)
            labels_train = pd.DataFrame(labels_train)
            data_test = pd.DataFrame(data_test)
            labels_test = pd.DataFrame(labels_test)

            if not os.path.exists(f"/{path2}/poet_output/{typez}_{version}/datasets"):
                os.makedirs(f"/{path2}/poet_output/{typez}_{version}/datasets")
            
            data_train.to_csv(f"/{path2}/poet_output/{typez}_{version}/datasets/data.csv",
                            index=False)
            labels_train.to_csv(f"/{path2}/poet_output/{typez}_{version}/datasets/label.csv",
                            index=False)
            data_test.to_csv(f"/{path2}/poet_output/{typez}_{version}/datasets/data_test.csv",
                            index=False)
            labels_test.to_csv(f"/{path2}/poet_output/{typez}_{version}/datasets/label_test.csv",
                            index=False)
    
    params = {
                "type": '1d_period',
                "epochs": 300,
                "batch_size": 100,
                "verbose": 2,
                "retrain": False,
                "threshold": 20,
                "path_to_store": path2,
                "version": version,#'joint',#'3348093',
                "features": None
            }
    
    for typez in types:
        params['type'] = typez
        params['features'] = [True, True, True, True, True, True, True, True, False, False, True]
        if typez == '1d_period':
            params['features'] = [True, True, True, True, True, True, True, True, False, False]
        poet = poet_solver.POET_IC_Solver(**params)
        data_test = pd.read_csv(f"/{path2}/poet_output/{typez}_{version}/datasets/data_test.csv")
        labels_test = pd.read_csv(f"/{path2}/poet_output/{typez}_{version}/datasets/label_test.csv")
        data_test = data_test.to_numpy()
        labels_test = labels_test.to_numpy()

        correct = []
        predicted = []
        difference = []
        error = []

        for i in range(len(data_test)): #10
            x_test = data_test[i]
            y_test = labels_test[i][0]
            x_poet_down, x_poet_up = poet.fit_evaluate(X_test=x_test)
            x_predict = 0.5*(x_poet_down+x_poet_up)
            x_predict = x_predict.numpy()[0][0]
            correct.append(y_test)
            predicted.append(x_predict)
            difference.append(numpy.abs(x_predict-y_test))
            error.append(numpy.abs(x_predict-y_test)*100/y_test)
        
        #print('Correct is ',correct)
        #print('Predicted is ',predicted)
        #print('Error is ',error)
        #print('Put \'em together and ',numpy.array([correct,predicted,error]))
        compiled = numpy.array([correct,predicted,difference,error]).T
        compiled = pd.DataFrame(compiled)

        compiled.to_csv(f"/{path2}/poet_output/{typez}_{version}/datasets/compiled.csv",
                        index=False)

def search_for_acceptable_systems():
    from CircularizationDissipationConstraints.source.bayesian.windemuth_et_al_util import get_summary_data,get_available_kic
    a = get_summary_data()
    b = get_available_kic()
    mass_low = 0.4
    mass_high = 1.2

    # Get the systems that are within the mass range
    # a is a pandas dataframe with relevant columns being 'maxlike_msum' and 'maxlike_mrat'
    acceptable_systems = []
    for index in b:#a.index:
        sum_of_masses = a.loc[index]['maxlike_msum']
        mass_ratio = a.loc[index]['maxlike_mrat']
        m1 = sum_of_masses/(1+mass_ratio)
        m2 = sum_of_masses - m1
        if mass_low < m1 < mass_high and mass_low < m2 < mass_high:
            acceptable_systems.append(index)
            #print('Index ',index,' has m1 of ',m1,' and m2 of ',m2)
            print(index)
    
    print('Acceptable systems are ',acceptable_systems)

if __name__ == '__main__':
    print("starting")
    # setup_process(
    #                 fname_datetime_format='%Y%m%d%H%M%S',
    #                 system='M',
    #                 std_out_err_fname='josh_output_2/{task}/{system}_{now}_{pid:d}.outerr',
    #                 logging_fname='josh_output_2/{task}/{system}_{now}_{pid:d}.log',
    #                 logging_verbosity='debug',
    #                 logging_message_format='%(levelname)s %(asctime)s %(name)s: %(message)s | %(pathname)s.%(funcName)s:%(lineno)d'#,
    #                 #logging_datetime_format=config.logging_datetime_format
    #               )

    # new_logic()
    just_make_triangles()
    #make_data()
    #trainnn()
    #check_for_bad()
    #validate(version='10091257')
    # a=time.time()
    # #validate(version='3241344')
    # #for i in range(10):
    # try_values(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,'/home/vortebo/ctime/crycrycry3','3241344')
    # b=time.time()
    # print('TIME TAKEN IS ',b-a)
    # a=time.time()
    # #validate(version='7691527')
    # #for i in range(10):
    # try_values(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,'/home/vortebo/ctime/crycrycry3','7691527')
    # b=time.time()
    # print('TIME TAKEN IS ',b-a)
    # a=time.time()
    # #validate(version='joint')
    # #for i in range(10):
    # try_values(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,'/home/vortebo/ctime/crycrycry3','joint')
    # b=time.time()
    # print('TIME TAKEN IS ',b-a)
    #join_data()
    #search_for_acceptable_systems()
    print("Done!")