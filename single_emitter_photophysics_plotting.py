#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:46:42 2020

@author: lisa-marianeedham
"""

def fn_intensitytrace(file_name,folder,data,time,i,x,y):
    """
    Plots photons versus time
    Inputs: data, filename and foldername which should be defined in the script
    
    """
    import numpy as np
    import matplotlib.pyplot as plt


    figure_name=file_name+'_intensity_trace'
    for a in range(0,len(data),int(len(data)/2)):
        if i==a:
            x_coord=x[i+1]
            y_coord=y[i+1]
            max_int=np.max(data[i])
            min_int=np.min(data[i])
            #norm_int = [b / max_int for b in data[i]]
            plt.figure()
            #plt.plot(time[0:len(time)-1],norm_int,'g')
            plt.plot(time[0:len(time)-1],data[i],'g')
            plt.xlim(0, 100)
            plt.ylim(min_int, (max_int+100))
            plt.xlabel('Time (s)', fontname='Arial', fontsize=12)
            plt.ylabel('Photon counts (photons)', fontname='Arial', fontsize=12)
            plt.xticks(fontname='Arial',fontsize=12)
            plt.yticks(fontname='Arial', fontsize=12)
            plt.savefig(folder+'/Figures/PDFs'+ '/' + figure_name + '_'+str(x_coord)+','+str(y_coord)+'.pdf', dpi=500)
            plt.savefig(folder+'/Figures/PNGs'+ '/' + figure_name + '_'+str(x_coord)+','+str(y_coord)+'.png', dpi=500)

    return (plt.show())



def fn_total_photon_hist(file_name,folder,total_photon):
    """
    Plots histogram of total detected photons and fits to lognormal distribution
    Inputs: data, filename and foldername which should be defined in the script
    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import lognorm
    from pylab import text
    
    n_molecules=len(total_photon)
    
    #Plot histogram
    figure_name=file_name+'_totalPhotons'
    ax = plt.subplot(111)
    num_bins = np.linspace(int(min(total_photon)), int(max(total_photon)), int(np.sqrt(len(total_photon))*3))
    ax.hist(total_photon, bins=num_bins, density=True,color='cornflowerblue',edgecolor='black')

    #Fit curve
    sigma,loc,mean = lognorm.fit(total_photon, floc=0)
    pdf = lognorm.pdf(num_bins, sigma, loc, mean) #sigma=shape, mu=np.log(scale)
    ax.plot(num_bins, pdf, 'k',linestyle='--')

    #Edit plot
    plt.xlabel('Total number of photons',fontname='Arial', fontsize=12)
    plt.ylabel('Probability density',fontname='Arial', fontsize=12)
    plt.xticks(fontname='Arial',fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    text(0.75, 0.95,'μ='+str(round(mean/10**6,2))+'$x10^6$ photons',horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,fontname='Arial', fontsize=12)
    text(0.40, 0.95,'N='+str(n_molecules),horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,fontname='Arial', fontsize=12)
    plt.savefig(folder+'/Figures/PDFs'+ '/' + figure_name + '.pdf', dpi=500)
    plt.savefig(folder+'/Figures/PNGs'+ '/' + figure_name + '.png', dpi=300)

    return (plt.show())

def fn_total_ontime_hist(file_name,folder,total_ontime):
    """
    Plots histogram of total on time and fits to lognormal distribution
    Inputs: data, filename and foldername which should be defined in the script
    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import lognorm
    from pylab import text
    
    n_molecules=len(total_ontime)
    
    #Plot histogram
    figure_name=file_name+'_totalOntime'
    ax = plt.subplot(111)
    num_bins = np.linspace(int(min(total_ontime)), int(max(total_ontime)), int(np.sqrt(len(total_ontime))*3))
    ax.hist(total_ontime, bins=num_bins, density=True,color='firebrick',edgecolor='black')

    #Fit curve
    sigma,loc,mean = lognorm.fit(total_ontime, floc=0)
    pdf = lognorm.pdf(num_bins, sigma, loc, mean) #sigma=shape, mu=np.log(scale)
    ax.plot(num_bins, pdf, 'k',linestyle='--')

    #Edit plot
    plt.xlabel('Total on time (s)', fontname='Arial', fontsize=12)
    plt.ylabel('Probability density', fontname='Arial', fontsize=12)
    plt.xticks(fontname='Arial',fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    text(0.75, 0.95,'μ='+str(round(mean,2))+' s',horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,fontname='Arial', fontsize=12)
    text(0.40, 0.95,'N='+str(n_molecules),horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,fontname='Arial', fontsize=12)
    plt.savefig(folder+'/Figures/PDFs'+ '/' + figure_name + '.pdf', dpi=500)
    plt.savefig(folder+'/Figures/PNGs'+ '/' + figure_name + '.png', dpi=500)
    
    return (plt.show())

def fn_total_blinks_hist(file_name,folder,total_blinks):
    
    """
    Plots histogram of total number of and fits to lognormal distribution
    Inputs: data, filename and foldername which should be defined in the script
    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import lognorm
    from pylab import text
    
    n_molecules=len(total_blinks)
    
    #Plot histogram
    figure_name=file_name+'_totalBlinks'
    ax = plt.subplot(111)
    num_bins = np.linspace(int(min(total_blinks)), int(max(total_blinks)), int(np.sqrt(len(total_blinks))*3))
    ax.hist(total_blinks, bins=num_bins, density=True,color='slateblue',edgecolor='black')
    
    #Fit curve
    sigma,loc,mean = lognorm.fit(total_blinks, floc=0)
    pdf = lognorm.pdf(num_bins, sigma, loc, mean) #sigma=shape, mu=np.log(scale)
    ax.plot(num_bins, pdf, 'k',linestyle='--')
    
    #Edit plot
    plt.xlabel('Number of blinks', fontname='Arial', fontsize=12)
    plt.ylabel('Probability density', fontname='Arial', fontsize=12)
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    text(0.75, 0.95,'μ='+str(round(mean,2))+' blinks',horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,fontname='Arial', fontsize=12)
    text(0.40, 0.95,'N='+str(n_molecules),horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,fontname='Arial', fontsize=12)
    plt.savefig(folder+'/Figures/PDFs'+ '/' + figure_name + '.pdf', dpi=500)
    plt.savefig(folder+'/Figures/PNGs'+ '/' + figure_name + '.png', dpi=500)
    
    return (plt.show())


def fn_photonflux_hist(file_name,folder,mean_photons_per_sec):
    
    """
    Plots histogram of total number of and fits to lognormal distribution
    Inputs: data, filename and foldername which should be defined in the script
    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import lognorm
    from pylab import text
    
    n_molecules=len(mean_photons_per_sec)
    
    #Plot photon flux
    figure_name=file_name+'_photonsPerSecond'
    ax = plt.subplot(111)
    num_bins = np.linspace(int(min(mean_photons_per_sec)), int(max(mean_photons_per_sec)), int(np.sqrt(len(mean_photons_per_sec))*4))
    ax.hist(mean_photons_per_sec, bins=num_bins, density=True, color='darkorange',edgecolor='black')
    
    #Fit curve
    sigma,loc,mean = lognorm.fit(mean_photons_per_sec, floc=0)
    pdf = lognorm.pdf(num_bins, sigma, loc, mean) #sigma=shape, mu=np.log(scale)
    ax.plot(num_bins, pdf, 'k',linestyle='--')
    
    #Edit plot
    plt.xlabel('Photon flux ($s^{-1}$)', fontname='Arial', fontsize=12)
    plt.ylabel('Probability density', fontname='Arial', fontsize=12)
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    text(0.75, 0.95,'μ='+str(round(mean,2))+' photons $s^{-1}$',horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,fontname='Arial', fontsize=12)
    text(0.40, 0.95,'N='+str(n_molecules),horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,fontname='Arial', fontsize=12)
    plt.savefig(folder+'/Figures/PDFs'+ '/' + figure_name + '.pdf', dpi=500)
    plt.savefig(folder+'/Figures/PNGs'+ '/' + figure_name + '.png', dpi=500)
    
    return (plt.show())


def fn_onstatetime_hist(file_name,folder,mean_onstatetime,distribution):
    
    """
    Plots histogram of total number of and fits to lognormal distribution
    Inputs: data, filename and foldername which should be defined in the script
    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import lognorm,norm
    from pylab import text
    
    n_molecules=len(mean_onstatetime)
    
    #Plot photon flux
    figure_name=file_name+'_onstatetime'
    ax = plt.subplot(111)
    num_bins = np.linspace(int(min(mean_onstatetime)), int(np.mean(mean_onstatetime)), int(np.sqrt(len(mean_onstatetime))*8))
    ax.hist(mean_onstatetime, bins=num_bins, density=True, color='forestgreen',edgecolor='black')
    
    #Choose distribution
    if distribution=='lognormal':
        #Fit lognormal curve
        sigma,loc,mean = lognorm.fit(mean_onstatetime, floc=0)
        pdf = lognorm.pdf(num_bins, sigma, loc, mean) #sigma=shape, mu=np.log(scale)
        ax.plot(num_bins, pdf, 'k',linestyle='--')
        
    elif distribution=='normal':
        #Fit normal curve
        mean, std = norm.fit(mean_onstatetime)
        pdf = norm.pdf(num_bins, mean, std) #sigma=shape, mu=np.log(scale)
        ax.plot(num_bins, pdf, 'k',linestyle='--')
        
        
    
    #Edit plot
    plt.xlabel('Mean on-state time (s)', fontname='Arial', fontsize=12)
    plt.ylabel('Probability density', fontname='Arial', fontsize=12)
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    text(0.75, 0.95,'μ='+str(round(mean,2))+' s',horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,fontname='Arial', fontsize=12)
    text(0.40, 0.95,'N='+str(n_molecules),horizontalalignment='center', verticalalignment='center',transform = ax.transAxes,fontname='Arial', fontsize=12)
    plt.savefig(folder+'/Figures/PDFs'+ '/' + figure_name + '.pdf', dpi=500)
    plt.savefig(folder+'/Figures/PNGs'+ '/' + figure_name + '.png', dpi=500)
    
    return (plt.show())


