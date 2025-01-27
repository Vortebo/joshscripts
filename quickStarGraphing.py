import numpy as np
import matplotlib.pyplot as plt

from orbital_evolution.star_interface import EvolvingStar
from stellar_evolution.manager import StellarEvolutionManager

age = np.linspace(0.2, 0.3, 100)

mass1 = 0.6747901609324751
metallicity = -0.11321753912301137
wind_strength = 0.17
wind_saturation_frequency = 2.45
diff_rot_coupling_timescale = 0.01
interpolator = StellarEvolutionManager('/home/vortebo/ctime/poet/stellar_evolution_interpolators').get_interpolator_by_name(
        'default'
    )
star1 = EvolvingStar(mass = mass1,
                    metallicity = metallicity,
                    wind_strength = wind_strength,
                    wind_saturation_frequency = wind_saturation_frequency,
                    diff_rot_coupling_timescale = diff_rot_coupling_timescale,
                    interpolator = interpolator)
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(age, star1.core_inertia(age))
plt.xlabel('Age')
plt.ylabel('Core Inertia')
plt.grid()
plt.subplot(2, 1, 2)
plt.plot(age, star1.envelope_inertia(age))
plt.xlabel('Age')
plt.ylabel('Envelope Inertia')
plt.grid()
plt.savefig('star1.png')
star1.delete()

mass2 = 0.7090123612219404
star2 = EvolvingStar(mass = mass2,
                    metallicity = metallicity,
                    wind_strength = wind_strength,
                    wind_saturation_frequency = wind_saturation_frequency,
                    diff_rot_coupling_timescale = diff_rot_coupling_timescale,
                    interpolator = interpolator)
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(age, star2.core_inertia(age))
plt.xlabel('Age')
plt.ylabel('Core Inertia')
plt.grid()
plt.subplot(2, 1, 2)
plt.plot(age, star2.envelope_inertia(age))
plt.xlabel('Age')
plt.ylabel('Envelope Inertia')
plt.grid()
plt.savefig('star2.png')
star2.delete()