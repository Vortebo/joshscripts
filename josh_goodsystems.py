import CircularizationDissipationConstraints.source.bayesian.windemuth_et_al_util as bob
import numpy
from astropy import units
from stellar_evolution.library_interface import \
    library as poet_stellar_evolution

gary = bob.get_available_kic()
gary.sort()
# print(gary)

# index = bob.get_summary_data().index.to_numpy()
# desired_index = [i for i in range(len(index)) if index[i] == 9468384]#'11404698'
# print(desired_index)

#######
# GOOD FOR INDIVIDUAL SYSTEMS
#######
# kic = 5622250
# wow=bob.get_summary_data().loc[kic]
# msumbob = wow['maxlike_msum']
# mratbob = wow['maxlike_mrat']
# m1bob = msumbob/(1+mratbob)
# m2bob = msumbob-m1bob
# metalbob = poet_stellar_evolution.feh_from_z(wow['posterior_z']) * units.dimensionless_unscaled
# agebob = wow['maxlike_age']
# pbob = wow['maxlike_period']
# e=(wow['maxlike_esinw']**2 + wow['maxlike_ecosw']**2)**(1/2)
# print(kic)
# print('m1',m1bob)
# print('m2',m2bob)
# print('metal',metalbob)
# print('age',agebob)
# print('P',pbob)
# print('e',e)

#gary = bob.get_available_kic()
# george = []
# for kic in gary:
#     wow=bob.get_summary_data().loc[kic]
#     msumbob = wow['maxlike_msum']
#     mratbob = wow['maxlike_mrat']
#     m1bob = msumbob/(1+mratbob)
#     m2bob = msumbob-m1bob
#     metalbob = poet_stellar_evolution.feh_from_z(wow['posterior_z']) * units.dimensionless_unscaled
#     agebob = wow['maxlike_age']
#     if (0.4 <= m1bob <= 1.2) and (0.4 <= m2bob <= 1.2) and (-1.014 <= metalbob <= 0.537):
#         george.append(kic)

# print(len(gary))
# print(len(george))
# geog=numpy.array(george)
# geog.sort()
# print(geog)


# samples=bob.get_samples(10292238)
# print(numpy.median(samples['P']))
# envelope=bob.eccentricity_envelope(numpy.median(samples['P']))
# print(envelope)
# cool=numpy.linspace(0,100,11)
# for i in cool:
#     print(i,bob.eccentricity_envelope(i))
# # wow = bob.get_summary_data().loc[10292238]
# # for item in wow.index:
# #     print(item, wow[item])


george = []
for kic in gary:
    wow=bob.get_summary_data().loc[kic]
    msumbob = wow['maxlike_msum']
    mratbob = wow['maxlike_mrat']
    m1bob = msumbob/(1+mratbob)
    m2bob = msumbob-m1bob
    metalbob = poet_stellar_evolution.feh_from_z(wow['posterior_z']) * units.dimensionless_unscaled
    agebob = wow['maxlike_age']

    samples=bob.get_samples(kic)
    envelope=bob.eccentricity_envelope(numpy.median(samples['P']))

    # if envelope > 0.7:
    #     george.remove(kic)

    if (0.4 <= m1bob <= 1.2) and (0.4 <= m2bob <= 1.2) and (-1.014 <= metalbob <= 0.537) and (envelope <= 0.7):
        george.append(kic)
    
print(len(gary))
print(len(george))
geog=numpy.array(george)
geog.sort()
print(geog)
