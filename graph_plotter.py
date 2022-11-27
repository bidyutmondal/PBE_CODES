import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import numpy as np

def plotNumberDensity():
    classSizeS = []
    numberDensityS = []

    numberDensityDataS = open('number_density_data_sfvm.txt')

    for row in numberDensityDataS:
        row = row.rstrip('\n').split(' ')
        classSizeS.append(float(row[0]))
        numberDensityS.append(float(row[1]))

    classSizeN = []
    numberDensityN = []

    numberDensityDataN = open('number_density_data_nfvm.txt')

    for row in numberDensityDataN:
        row = row.rstrip('\n').split(' ')
        classSizeN.append(float(row[0]))
        numberDensityN.append(float(row[1]))

    plt.plot(range(len(classSizeN)), numberDensityS, marker='.', markersize=4, label='SFVM')
    plt.plot(range(len(classSizeN)), numberDensityN, marker='.', markersize=4, label='NFVM')
    # plt.xticks(range(len(classSizeN)),classSizeN)
    plt.xlabel('dimensionless time')
    plt.ylabel('number density')
    plt.legend(loc='right')
    plt.ylim([-0.02, 1.02])
    plt.show()

def plotNumber():
    classSizeS = []
    numberS = []

    numberDataS = open('number_of_particles_data_sfvm.txt')

    for row in numberDataS:
        row = row.rstrip('\n').split(' ')
        classSizeS.append(float(row[0]))
        numberS.append(float(row[1]))

    classSizeN = []
    numberN = []

    numberDataN = open('number_of_particles_data_nfvm.txt')

    for row in numberDataN:
        row = row.rstrip('\n').split(' ')
        classSizeN.append(float(row[0]))
        numberN.append(float(row[1]))

    plt.plot(range(len(classSizeN)), numberN, marker='.', markersize=4, label='NFVM')
    plt.plot(range(len(classSizeN)), numberS, marker='.', markersize=4, label='SFVM')
    plt.xticks(range(len(classSizeN)))
    plt.xlabel('dimensionless size')
    plt.ylabel('number of particles ')
    plt.legend(loc='right')
    plt.show()

def plotFirstMoment():
    fisrtMomentTimeS = []
    firstMomentS = []

    firstMomentDataS = open('first_moment_data_sfvm.txt', 'r')

    for row in firstMomentDataS:
        row = row.rstrip('\n').split(' ')
        fisrtMomentTimeS.append(float(row[0]))
        firstMomentS.append(float(row[1]))

    fisrtMomentTimeN = []
    firstMomentN = []

    firstMomentDataN = open('first_moment_data_nfvm.txt', 'r')

    for row in firstMomentDataN:
        row = row.rstrip('\n').split(' ')
        fisrtMomentTimeN.append(float(row[0]))
        firstMomentN.append(float(row[1]))

    plt.plot(fisrtMomentTimeS, firstMomentS, marker='.', markersize=4, label='SFVM')
    plt.plot(fisrtMomentTimeN, firstMomentN, marker='.', markersize=4, label='NFVM')
    plt.xlabel('dimensionless time')
    plt.ylabel('normalised first moment')
    plt.legend(loc='right')
    plt.ylim([0, 1.02])
    plt.show()

def plotZerothMoment():
    zerothMomentTimeS = []
    zerothMomentS = []

    zerothMomentDataS = open('zeroth_moment_data_sfvm.txt', 'r')

    for row in zerothMomentDataS:
        row = row.rstrip('\n').split(' ')
        zerothMomentTimeS.append(float(row[0]))
        zerothMomentS.append(float(row[1]))

    zerothMomentTimeN = []
    zerothMomentN = []

    zerothMomentDataN = open('zeroth_moment_data_nfvm.txt', 'r')

    for row in zerothMomentDataN:
        row = row.rstrip('\n').split(' ')
        zerothMomentTimeN.append(float(row[0]))
        zerothMomentN.append(float(row[1]))


    plt.plot(zerothMomentTimeS, zerothMomentS, marker='.', markersize=4, label='SFVM')
    plt.plot(zerothMomentTimeN, zerothMomentN, marker='.', markersize=4, label='NFVM')
    plt.xlabel('dimensionless time')
    plt.ylabel('normalised zeroth moment')
    plt.legend(loc='right')
    plt.ylim([-0.02, 1.02])
    plt.show()

# plotNumberDensity()
plotNumber()
# plotZerothMoment()
# plotFirstMoment()