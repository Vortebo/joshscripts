from types import SimpleNamespace
from astropy import units
from stellar_evolution.library_interface import MESAInterpolator
MESAInterpolator.set_quantity_lower_limit('iconv', 1e-5)
from stellar_evolution.manager import StellarEvolutionManager
import numpy
from orbital_evolution.transformations import phase_lag
from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from multiprocessing_util import setup_process

from multiprocessing import Manager
thelock = Manager().Lock()

from orbital_evolution.binary import Binary
from orbital_evolution.star_interface import EvolvingStar
#from orbital_evolution.planet_interface import LockedPlanet

import logging

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
            #2.0 * numpy.pi / 50.0
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

def runsystem(lgQ_min = 6.730834951374322,
                         lgQ_break_period = 0.9524689224317893 * units.day,
                         lgQ_powerlaw = 0.6928867615956023,
                         age = 4.992647080579888,
                         feh = 0.16115155325450198,
                         orbital_period = 4.0532537066910015,
                         primary_mass = 1.1454357634717733,
                         secondary_mass = 0.47249084407776804,
                         cmd_primary_radius = 1.3876509776305228,
                         cmd_secondary_radius = 0.4639166819029762,
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
    #MESAInterpolator.set_quantity_lower_limit('iconv', 1e-100)
    interpolator = StellarEvolutionManager('/home/vortebo/ctime/poet/stellar_evolution_interpolators').get_interpolator_by_name(
            'default'
        )

    logger.debug('Now loading: eccentricity expansion coefficients.')
    import os
    orbital_evolution_library.prepare_eccentricity_expansion(
        os.path.expanduser(
            '~/eccentricity_expansion_coef_O400.sqlite'
        ).encode('ascii'),
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
        current_age=system.age,
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
    initial_secondary_angmom = [angmom_a,angmom_b]#value_finder.get_secondary_initial_angmom()

    value_finder.try_system([initial_porb,initial_eccentricity,0.0],initial_secondary_angmom)

def runsystem2(lgQ_min = 6.730834951374322,
                         lgQ_break_period = 0.9524689224317893 * units.day,
                         lgQ_powerlaw = 0.6928867615956023,
                         age = 4.992647080579888,
                         feh = 0.16115155325450198,
                         orbital_period = 4.0532537066910015,
                         primary_mass = 1.1454357634717733,
                         secondary_mass = 0.47249084407776804,
                         cmd_primary_radius = 1.3876509776305228,
                         cmd_secondary_radius = 0.4639166819029762,
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
    #MESAInterpolator.set_quantity_lower_limit('iconv', 1e-100)
    interpolator = StellarEvolutionManager('/home/vortebo/ctime/poet/stellar_evolution_interpolators').get_interpolator_by_name(
            'default'
        )

    logger.debug('Now loading: eccentricity expansion coefficients.')
    import os
    orbital_evolution_library.prepare_eccentricity_expansion(
        os.path.expanduser(
            '~/eccentricity_expansion_coef_O400.sqlite'
        ).encode('ascii'),
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

    value_finder = SimpleNamespace(
        system=system,
        interpolator=interpolator,
        current_age=system.age,
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
        extra_evolve_args=extra_evolve_args
    )
    initial_secondary_angmom = [angmom_a,angmom_b]

    try_system(value_finder, [initial_porb,initial_eccentricity,0.0], initial_secondary_angmom)

def _create_star(mass,
                feh,
                interpolator,
                dissipation=None,
                *,
                wind_strength=0.17,
                wind_saturation_frequency=2.78,
                diff_rot_coupling_timescale=5e-3,
                interpolation_age=None):
    """
    Create the star to use in the evolution.

    Args:
        mass:    The mass of the star to create, along with astropy units.

        feh:    The [Fe/H] value of the star to create.

        interpolator:    POET stellar evolution interpolator giving the
            evolution of the star's properties.

        dissipation:    If None, no dissipation is set. Otherwise, a
            dictionary of keyword arguments to pass to
            `EvolvingStar.set_dissipation()`.

        wind_strength:    See same name argument to EvolvingStar.__init__()

        wind_saturation_frequency:    See same name argument to
            EvolvingStar.__init__()

        diff_rot_coupling_timescale:    See same name argument to
            EvolvingStar.__init__()

        interpolation_age:    The age at which to initialize the
            interpolation for the star. If `None`, the core formation age is
            used.

    Returns:
        EvolvingStar:
            The star in the system useable for calculating obital evolution.
    """

    #False positive
    #pylint: disable=no-member
    logger = logging.getLogger(__name__)
    star = EvolvingStar(
        mass=mass.to_value(units.M_sun),
        metallicity=feh,
        wind_strength=wind_strength,
        wind_saturation_frequency=wind_saturation_frequency,
        diff_rot_coupling_timescale=(
            diff_rot_coupling_timescale.to_value(units.Gyr)
        ),
        interpolator=interpolator
    )
    #pylint: enable=no-member
    star.select_interpolation_region(star.core_formation_age()
                                        if interpolation_age is None else
                                        interpolation_age)
    logger.debug('By the way, the interpolation age is %s', repr(interpolation_age))
    if dissipation is not None:
        star.set_dissipation(zone_index=0, **dissipation)
    return star

def _create_primary(system,interpolator,configuration):
    """Create the primary object in the system."""

    mprimary = getattr(system,
                        'Mprimary',
                        system.primary_mass)
    
    configuration['dissipation']['primary']
    configuration['disk_dissipation_age']
    configuration['primary_wind_strength']
    configuration['primary_wind_saturation']
    configuration['primary_core_envelope_coupling_timescale']
    
    logger = logging.getLogger(__name__)
    logger.debug('The dissipation we are setting for the primary is %s', repr(configuration['dissipation']['primary']))

    return _create_star(
        mprimary,
        system.feh,
        interpolator['primary'],
        configuration['dissipation']['primary'],
        wind_strength=configuration['primary_wind_strength'],
        wind_saturation_frequency=(
            configuration['primary_wind_saturation']
        ),
        diff_rot_coupling_timescale=(
            configuration['primary_core_envelope_coupling_timescale']
        )
    )

def _create_secondary(system,interpolator,configuration):
    """Create the two objects comprising the system to evolve."""

    msecondary = getattr(system,
                            'Msecondary',
                            system.secondary_mass)
    
    configuration['dissipation']['secondary']
    configuration['disk_dissipation_age']
    configuration['secondary_wind_strength']
    configuration['secondary_wind_saturation']
    configuration['secondary_core_envelope_coupling_timescale']

    logger = logging.getLogger(__name__)
    logger.debug('The dissipation we are setting for the secondary is %s', repr(configuration['dissipation']['secondary']))

    interpolation_age = configuration['disk_dissipation_age']
    return _create_star(
        msecondary,
        system.feh,
        interpolator['secondary'],
        configuration['dissipation']['secondary'],
        wind_strength=configuration['secondary_wind_strength'],
        wind_saturation_frequency=(
            configuration['secondary_wind_saturation']
        ),
        diff_rot_coupling_timescale=configuration[
            'secondary_core_envelope_coupling_timescale'
        ],
        interpolation_age=interpolation_age
    )

def _create_system(target_state,disk_dissipation_age,primary,
                    secondary,
                    *,
                    porb_initial,
                    initial_eccentricity,
                    initial_obliquity,
                    initial_secondary_angmom):
    """
    Create the system to evolve from the two bodies (primary & secondary).

    Args:
        primary:    The primary in the system. Usually created by calling
            :meth:`_create_star`.

        planet:    The secondary in the system. Usually created by calling
            :meth:`_create_star` or :meth:`_create_planet`.

        porb_initial:    Initial orbital period in days.

        initial_eccentricity:    The initial eccentricity of the system.

        initial_obliquity:     The initial obliquity to assume for the
            system in rad.

    Returns:
        Binary:
            The binary system ready to evolve.
    """

    #False positive
    #pylint: disable=no-member

    logger = logging.getLogger(__name__)


    binary = Binary(
        primary=primary,
        secondary=secondary,
        initial_orbital_period=porb_initial,
        initial_eccentricity=initial_eccentricity,
        initial_inclination=initial_obliquity,
        disk_lock_frequency=(2.0 * numpy.pi
                                /
                                target_state.Pdisk),
        disk_dissipation_age=disk_dissipation_age,
        secondary_formation_age=target_state.planet_formation_age
    )
    #pylint: enable=no-member
    binary.configure(age=primary.core_formation_age(),
                        semimajor=float('nan'),
                        eccentricity=float('nan'),
                        spin_angmom=numpy.array([0.0]),
                        inclination=None,
                        periapsis=None,
                        evolution_mode='LOCKED_SURFACE_SPIN')
    if isinstance(secondary, EvolvingStar):
        initial_obliquity = numpy.array([0.0])
        initial_periapsis = numpy.array([0.0])
    else:
        initial_obliquity = None
        initial_periapsis = None
    secondary.configure(
        #False positive
        #pylint: disable=no-member
        age=target_state.planet_formation_age,
        #pylint: enable=no-member
        companion_mass=primary.mass,
        #False positive
        #pylint: disable=no-member
        semimajor=binary.semimajor(porb_initial),
        #pylint: enable=no-member
        eccentricity=initial_eccentricity,
        spin_angmom=numpy.array(initial_secondary_angmom),
        inclination=initial_obliquity,
        periapsis=initial_periapsis,
        locked_surface=False,
        zero_outer_inclination=True,
        zero_outer_periapsis=True
    )
    primary.detect_stellar_wind_saturation()
    if isinstance(secondary, EvolvingStar):
        secondary.detect_stellar_wind_saturation()
    return binary

def _format_evolution(binary, interpolators, secondary_star):
    """
    Return pre-calculated evolution augmented with stellar quantities.

    Args:
        binary:    The binary used to calculate the evolution to format.

        interpolators:


    The returned evolution contains the full evolution produced by
    Binary.get_evolution(), as well as the evolution of the following
    quantities:

        * **orbital_period**: the orbital period

        * **(primary/secondary)_radius: The radius of the primary/secondary
            star.

        * **(primary/secondary)_lum: The luminosity of the primary/secondary
            star.

        * **(primary/secondary)_(iconv/irad): The convective/radiative
            zone moment of inertia of the primary/secondary star.
    """

    logger = logging.getLogger(__name__)
    logger.debug('before we get evolution')
    evolution = binary.get_evolution()
    logger.debug('after we get evolution')
    #False positive
    #pylint: disable=no-member
    evolution.orbital_period = binary.orbital_period(evolution.semimajor)
    logger.debug('we did orbital period')

    components_to_get = ['primary']
    if secondary_star:
        components_to_get.append('secondary')
    logger.debug('we are about to get the components')

    for component in components_to_get:

        if (
                len(interpolators[component].track_masses) == 1
                and
                len(interpolators[component].track_feh) == 1
        ):
            star_params = dict(
                mass=interpolators[component].track_masses[0],
                feh=interpolators[component].track_feh[0]
            )
        else:
            star_params = dict(
                mass=getattr(binary, component).mass,
                feh=getattr(binary, component).metallicity
            )

        for quantity_name in ['radius', 'lum', 'iconv', 'irad']:
            quantity = interpolators[component](
                quantity_name,
                **star_params
            )
            values = numpy.full(evolution.age.shape, numpy.nan)

            #TODO: determine age range properly (requires C/C++ code
            #modifications)
            valid_ages = numpy.logical_and(
                evolution.age > quantity.min_age * 2.0,
                evolution.age < quantity.max_age
            )
            if quantity_name in ['iconv', 'irad']:
                values[valid_ages] = getattr(
                    getattr(binary, component),
                    (
                        ('envelope' if quantity_name == 'iconv' else 'core')
                        +
                        '_inertia'
                    )
                )(evolution.age[valid_ages])
            else:
                values[valid_ages] = quantity(evolution.age[valid_ages])
            setattr(evolution,
                    component + '_' + quantity_name,
                    values)
    logger.debug('we are about to return evolution')

    return evolution

def try_system(value_finder,initial_conditions,initial_secondary_angmom):

        system = value_finder.system
        interpolator=value_finder.interpolator
        current_age=system.age
        disk_period=value_finder.disk_period
        initial_obliquity=value_finder.initial_obliquity
        disk_dissipation_age=value_finder.disk_dissipation_age
        dissipation=value_finder.dissipation
        secondary_star=value_finder.secondary_star
        primary_wind_strength=value_finder.primary_wind_strength
        primary_wind_saturation=value_finder.primary_wind_saturation
        primary_core_envelope_coupling_timescale=value_finder.primary_core_envelope_coupling_timescale
        secondary_wind_strength=value_finder.secondary_wind_strength
        secondary_wind_saturation=value_finder.secondary_wind_saturation
        secondary_core_envelope_coupling_timescale=value_finder.secondary_core_envelope_coupling_timescale
        secondary_disk_period=value_finder.secondary_disk_period
        orbital_period_tolerance=value_finder.orbital_period_tolerance
        period_search_factor=value_finder.period_search_factor
        scaled_period_guess=value_finder.scaled_period_guess

        porb = system.orbital_period
        target_state = SimpleNamespace(
            #False positive
            #pylint: disable=no-member
            age=current_age.to(units.Gyr).value,
            Porb=porb.to(units.day).value,
            Pdisk=disk_period.to(units.day).value,
            planet_formation_age=disk_dissipation_age.to(units.Gyr).value,
            evolution_max_time_step=value_finder.extra_evolve_args['max_time_step'],
            evolution_precision=value_finder.extra_evolve_args['precision']
            #pylint: enable=no-member
        )

        interpolator = dict(primary=interpolator,
                                    secondary=interpolator)
        configuration = dict(
            #False positive
            #pylint: disable=no-member
            disk_dissipation_age=disk_dissipation_age.to(units.Gyr).value,
            #pylint: enable=no-member
            orbital_period_tolerance=orbital_period_tolerance,
            dissipation=dissipation,
            initial_obliquity=initial_obliquity,
            primary_wind_strength=primary_wind_strength,
            primary_wind_saturation=primary_wind_saturation,
            primary_core_envelope_coupling_timescale=(
                primary_core_envelope_coupling_timescale
            ),
            secondary_core_envelope_coupling_timescale=(
                secondary_core_envelope_coupling_timescale
            ),
            secondary_wind_strength=secondary_wind_strength,
            secondary_wind_saturation=secondary_wind_saturation,
            secondary_disk_period=(secondary_disk_period
                                   or
                                   disk_period).to_value(units.day),
            period_search_factor=period_search_factor,
            scaled_period_guess=scaled_period_guess
        )
        
        initial_orbital_period=initial_conditions[0]
        initial_eccentricity=initial_conditions[1]
        initial_obliquity=initial_conditions[2]

        logger = logging.getLogger(__name__)

        primary = _create_primary(system,interpolator,configuration)
        secondary = _create_secondary(system,interpolator,configuration)
        if (primary.core_inertia(configuration['disk_dissipation_age']) == 0
            or
            secondary.core_inertia(configuration['disk_dissipation_age']) == 0
            ):
            logger.warning(
                'Primary or secondary core inertia is zero at current disk dissipation age: %s, ',
                repr(configuration['disk_dissipation_age'])
            )
            configuration['disk_dissipation_age'] = 0.02
        if not primary.core_inertia(configuration['disk_dissipation_age']) > 0:
            logger.error(
                'Reported primary core inertia at disk dissipation age: %s, ',
                repr(primary.core_inertia(configuration['disk_dissipation_age']))
            )
            raise ValueError("Primary core inertia is zero. Primary has not formed.",0)
        if not secondary.core_inertia(configuration['disk_dissipation_age']) > 0:
            logger.error(
                'Reported secondary core inertia at disk dissipation age: %s, ',
                repr(secondary.core_inertia(configuration['disk_dissipation_age']))
            )
            raise ValueError("Secondary core inertia is zero. Secondary has not formed.",0)

        binary=_create_system(target_state,configuration['disk_dissipation_age'],
            primary,
            secondary,
            #False positive
            #pylint: disable=no-member
            porb_initial=initial_orbital_period,
            #pylint: enable=no-member
            initial_eccentricity=initial_eccentricity,
            initial_obliquity=initial_obliquity,
            initial_secondary_angmom=initial_secondary_angmom
        )

        #if max_age is None:
        max_age = target_state.age
        #else:
        #    max_age = max_age.to(units.Gyr).value

        binary.evolve(
            max_age,
            target_state.evolution_max_time_step,
            target_state.evolution_precision,
            None,
            timeout=3600
            )

        final_state=binary.final_state()
        logger.debug('Final state age: %s, ',
                                          repr(final_state.age))
        logger.debug('Target state age: %s, ',
                                            repr(target_state.age))
        logger.debug('Initial eccentricity: %s, ',
                                            repr(initial_eccentricity))
        logger.debug('Final eccentricity: %s, ',
                                            repr(final_state.eccentricity))
        logger.debug('Initial period: %s, ',
                                            repr(initial_orbital_period))

        final_orbital_period = binary.orbital_period(final_state.semimajor)
        result = SimpleNamespace()
        setattr(result, 'orbital_period', numpy.array([final_orbital_period]))
        setattr(result, 'eccentricity', numpy.array([final_state.eccentricity]))
        
        # Clean up
        primary.delete()
        secondary.delete()
        binary.delete()

        logger.debug('Final period: %s, ',
                                            repr(result.orbital_period[-1]))
        try:
            logger.debug('ihhhh')
            assert(final_state.age==target_state.age)
        except AssertionError:
            logger.debug('uhhhhhhhhhhh')
            # Save the parameters and evolution to an astropy fits file. Parameters in header data.
            import os
            import datetime
            logger.debug('why')

            filename = 'failed_solutions'
            
            # Create the directory if it doesn't exist
            logger.debug('are')
            os.makedirs(filename, exist_ok=True)
            logger.debug('you')

            # Create the filename
            now = datetime.datetime.now()
            logger.debug('breaking')
            filename = filename + f'/solution_{now.strftime("%Y-%m-%d_%H-%M-%S")}.fits'
            logger.debug('here??????????????')

            evolution = _format_evolution(binary,
                                           interpolator,
                                           secondary_star)
            
            logger.debug('buddy??????????')
            return evolution#self._format_evolution(binary,
                            #               self.interpolator,
                            #               self.secondary_star)#.primary_iconv

            # Raise the error
            raise AssertionError(f"Final age does not match target age.")
        
        return result

if __name__ == '__main__':
    systemname='6521542'
    setup_process(
                    fname_datetime_format='%Y%m%d%H%M%S',
                    system=systemname,
                    std_out_err_fname='josh_output_8/{task}/{system}_{now}_{pid:d}.outerr',
                    logging_fname='josh_output_8/{task}/{system}_{now}_{pid:d}.log',
                    logging_verbosity='debug',
                    logging_message_format='%(levelname)s %(asctime)s %(name)s: %(message)s | %(pathname)s.%(funcName)s:%(lineno)d'
                  )

    runsystem(
        #4947726 (not fixed by iconv)
        # lgQ_min = 7.185163961316851,
        # lgQ_break_period = 1.2105510858211672 * units.day,
        # lgQ_powerlaw = -1.8675025503172198,
        # age = 4.537071400653452,
        # feh = -0.1476147128380225,
        # orbital_period = 4.726086038265389,
        # primary_mass = 1.0349127223165384,
        # secondary_mass = 0.5186184463883605,
        # cmd_primary_radius = 1.114127831417474,
        # cmd_secondary_radius = 0.5075721251761187,
        # final_eccentricity = 0.0,
        # initial_porb = 22.834226476460383,
        # initial_eccentricity = 0.8,
        # angmom_a = 0.10428821085789115,
        # angmom_b = 0.00019044827608288606
        #6521542 (yes fixed? :O)
        lgQ_min = 9.222211210635884,
        lgQ_break_period = 0.8597019976850838 * units.day,
        lgQ_powerlaw = 4.012989378277558,
        age = 1.5789897540930298,
        feh = -0.8460936639710941,
        orbital_period = 4.425754132228538,
        primary_mass = 1.0861970436931017,
        secondary_mass = 1.043207183779346,
        cmd_primary_radius = 1.0887086343631127,
        cmd_secondary_radius = 1.043463796711456,
        final_eccentricity = 0.0,
        initial_porb = 13.277262396685614,
        initial_eccentricity = 0.8,
        angmom_a = 0.020857808807702508,
        angmom_b = 0.40876059990796304
    )

    print("Done!")