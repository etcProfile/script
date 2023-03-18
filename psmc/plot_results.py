#!/usr/bin/python

###############################################################################
# This script is used for doing the plot of the demographic history of        #
# a random-mating population from a ms command. At the same time, the script  #
# allows to plot (in the same figure) the demographic history infered by the  #
# PSMC software.                                                              #
###############################################################################

import matplotlib.pyplot as plt

# Set the values of these global variables
#==============================================================================
# The original ms command:

# Path to the output file comming from the PSMC-ZAKS
PSMC_RESULTS_JSG01 = "C:\\Personal\\R_pro\\psmc\\JSG01.psmc"
PSMC_RESULTS_PTS01 = "C:\\Personal\\R_pro\\psmc\\PTS01.psmc"
PSMC_RESULTS_ZJJ01 = "C:\\Personal\\R_pro\\psmc\\ZJJ01.psmc"
# Path to the output file comming from the PSMC-TWR
PSMC_RESULTS_AWS01 = "C:\\Personal\\R_pro\\psmc\\AWS01.psmc"
PSMC_RESULTS_CPC02 = "C:\\Personal\\R_pro\\psmc\\CPC02.psmc"
PSMC_RESULTS_GTC01 = "C:\\Personal\\R_pro\\psmc\\GTC01.psmc"
PSMC_RESULTS_XBG01 = "C:\\Personal\\R_pro\\psmc\\XBG01.psmc"
PSMC_RESULTS_YMS01 = "C:\\Personal\\R_pro\\psmc\\YMS01.psmc"
PSMC_RESULTS_ZGB01 = "C:\\Personal\\R_pro\\psmc\\ZGB01.psmc"
#Path to the output file comming from the PSMC-MSC
PSMC_RESULTS_BTJ01 = "C:\\Personal\\R_pro\\psmc\\BTJ01.psmc"
PSMC_RESULTS_JXZ01 = "C:\\Personal\\R_pro\\psmc\\JXZ01.psmc"
PSMC_RESULTS_LYS01 = "C:\\Personal\\R_pro\\psmc\\LYS01.psmc"

PSMC_RESULTS_NLB01 = "C:\\Personal\\R_pro\\psmc\\NLB01.psmc"
PSMC_RESULTS_TTS01 = "C:\\Personal\\R_pro\\psmc\\TTS01.psmc"
PSMC_RESULTS_WYL01 = "C:\\Personal\\R_pro\\psmc\\WYL01.psmc"
#Path to the output file comming from the PSMC-HSJ
PSMC_RESULTS_FSX02 = "C:\\Personal\\R_pro\\psmc\\FSX02.psmc"
PSMC_RESULTS_GZX01 = "C:\\Personal\\R_pro\\psmc\\GZX01.psmc"
PSMC_RESULTS_SNC01 = "C:\\Personal\\R_pro\\psmc\\SNC01.psmc"
# Bin size used to generate the imput of PSMC (default is 100)
BIN_SIZE = 100

# Mutation rate per base per generation
MUTATION_RATE = 4e-9

# Number of years per generation
GENERAITON_TIME = 10

# Size of the plot
X_MIN = 1.55e5
X_MAX = 2e8
Y_MIN = 1e5
Y_MAX = 9.5e5

# What plot to do
PLOT_MS = False
PLOT_PSMC_RESULTS_JSG01 = True
PLOT_PSMC_RESULTS_PTS01 = True
PLOT_PSMC_RESULTS_ZJJ01 = True
PLOT_PSMC_RESULTS_AWS01 = True
PLOT_PSMC_RESULTS_CPC02 = True
PLOT_PSMC_RESULTS_GTC01 = True
PLOT_PSMC_RESULTS_XBG01 = True
PLOT_PSMC_RESULTS_YMS01 = True
PLOT_PSMC_RESULTS_ZGB01 = True
PLOT_PSMC_RESULTS_BTJ01 = True
PLOT_PSMC_RESULTS_JXZ01 = True
PLOT_PSMC_RESULTS_LYS01 = True
PLOT_PSMC_RESULTS_NLB01 = True
PLOT_PSMC_RESULTS_TTS01 = True
PLOT_PSMC_RESULTS_WYL01 = True
PLOT_PSMC_RESULTS_FSX02 = True
PLOT_PSMC_RESULTS_GZX01 = True
PLOT_PSMC_RESULTS_SNC01 = True
#==============================================================================


def psmc2funJSG01(filename=PSMC_RESULTS_JSG01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_JSG01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funPTS01(filename=PSMC_RESULTS_PTS01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_PTS01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funZJJ01(filename=PSMC_RESULTS_ZJJ01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_ZJJ01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funAWS01(filename=PSMC_RESULTS_AWS01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_AWS01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funCPC02(filename=PSMC_RESULTS_CPC02, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_CPC02, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funGTC01(filename=PSMC_RESULTS_GTC01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_GTC01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funXBG01(filename=PSMC_RESULTS_XBG01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_XBG01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funYMS01(filename=PSMC_RESULTS_YMS01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_YMS01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funZGB01(filename=PSMC_RESULTS_ZGB01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_ZGB01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funBTJ01(filename=PSMC_RESULTS_BTJ01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_BTJ01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funJXZ01(filename=PSMC_RESULTS_JXZ01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_JXZ01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funLYS01(filename=PSMC_RESULTS_LYS01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_LYS01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funNLB01(filename=PSMC_RESULTS_NLB01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_NLB01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funTTS01(filename=PSMC_RESULTS_TTS01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_TTS01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funWYL01(filename=PSMC_RESULTS_WYL01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_WYL01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funFSX02(filename=PSMC_RESULTS_FSX02, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_FSX02, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funGZX01(filename=PSMC_RESULTS_GZX01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_GZX01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

def psmc2funSNC01(filename=PSMC_RESULTS_SNC01, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(PSMC_RESULTS_SNC01, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)


if __name__ == "__main__":
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ###ZAKS-JSG01,PTS01,ZJJ01-gold
    if PLOT_PSMC_RESULTS_JSG01:
        (estimated_times, estimated_sizes) = psmc2funJSG01(PSMC_RESULTS_JSG01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='gold', label = "ZAKS")
    
    if PLOT_PSMC_RESULTS_PTS01:
        (estimated_times, estimated_sizes) = psmc2funPTS01(PSMC_RESULTS_PTS01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='gold')
    
    if PLOT_PSMC_RESULTS_ZJJ01:
        (estimated_times, estimated_sizes) = psmc2funZJJ01(PSMC_RESULTS_ZJJ01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='gold')
    ###TWR-AWS01,CPC02,GTC01,XBG01,YMS01,ZGB01-goldenrod
    if PLOT_PSMC_RESULTS_AWS01:
        (estimated_times, estimated_sizes) = psmc2funAWS01(PSMC_RESULTS_AWS01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='goldenrod', label = "TWR")

    if PLOT_PSMC_RESULTS_CPC02:
        (estimated_times, estimated_sizes) = psmc2funCPC02(PSMC_RESULTS_ZJJ01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='goldenrod')

    if PLOT_PSMC_RESULTS_GTC01:
        (estimated_times, estimated_sizes) = psmc2funXBG01(PSMC_RESULTS_ZJJ01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='goldenrod')

    if PLOT_PSMC_RESULTS_XBG01:
        (estimated_times, estimated_sizes) = psmc2funXBG01(PSMC_RESULTS_ZJJ01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='goldenrod')

    if PLOT_PSMC_RESULTS_YMS01:
        (estimated_times, estimated_sizes) = psmc2funYMS01(PSMC_RESULTS_ZJJ01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='goldenrod')
    ###MSC-BTJ01,JXZ01,LYS01,NLB01-orangered
    if PLOT_PSMC_RESULTS_BTJ01:
        (estimated_times, estimated_sizes) = psmc2funBTJ01(PSMC_RESULTS_BTJ01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='orangered', label = "MSC")

    if PLOT_PSMC_RESULTS_JXZ01:
        (estimated_times, estimated_sizes) = psmc2funJXZ01(PSMC_RESULTS_JXZ01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='orangered')

    if PLOT_PSMC_RESULTS_LYS01:
        (estimated_times, estimated_sizes) = psmc2funLYS01(PSMC_RESULTS_LYS01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='orangered')

    if PLOT_PSMC_RESULTS_NLB01:
        (estimated_times, estimated_sizes) = psmc2funNLB01(PSMC_RESULTS_NLB01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='orangered')

    if PLOT_PSMC_RESULTS_TTS01:
        (estimated_times, estimated_sizes) = psmc2funTTS01(PSMC_RESULTS_TTS01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='salmon', label = "MEC")

    if PLOT_PSMC_RESULTS_WYL01:
        (estimated_times, estimated_sizes) = psmc2funWYL01(PSMC_RESULTS_WYL01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='salmon')

    if PLOT_PSMC_RESULTS_FSX02:
        (estimated_times, estimated_sizes) = psmc2funFSX02(PSMC_RESULTS_FSX02, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='palegreen', label = "HSJ")

    if PLOT_PSMC_RESULTS_GZX01:
        (estimated_times, estimated_sizes) = psmc2funGZX01(PSMC_RESULTS_GZX01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='palegreen')

    if PLOT_PSMC_RESULTS_SNC01:
        (estimated_times, estimated_sizes) = psmc2funGZX01(PSMC_RESULTS_SNC01, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', color='palegreen')


    ax.set_xlabel("Time in years (10 years/generation)")
    ax.set_ylabel("Effective size (x 10^4)")
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax.grid(True)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_xscale('log')
    plt.legend(loc = 'best')
    
    plt.show()

