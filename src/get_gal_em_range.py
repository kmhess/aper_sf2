"""
TNH, 05/2025

Script to get a range of channels with high Galactic HI
reads in the data cube of the field, measures and plots RMS noise in each channel, 
and finds range channels with high RMS + expected Galactic Emission in a certain range of frequencies,
corresponding to Galactic emission.

Run within parent directory of field directories, i.e., in mos_<>/..
The galactic_emission_vel_ranges.py script should be in the same directory.

Usage:
1. from get_gal_em_range import find_rms_range
2. channel_range = find_rms_range(<loc_dir>, <field_name>, <splinefits>)
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import constants as const
from astropy.wcs import WCS
from astropy import wcs
from astropy.coordinates import SkyCoord, Supergalactic, Galactic, match_coordinates_sky
from spectral_cube import SpectralCube
from get_exp_gal_em import get_lb, get_gal_vel
import os.path
import os, sys
import astropy.table as table
from astropy.table import Table, QTable
from astropy.io import ascii
from io import StringIO
from io import BytesIO
import glob
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['text.usetex'] = False
matplotlib.rcParams['font.size'] = 15

#fits header params, which I think are the same for all cube3 fields (check in CARTA as you read them in, 
#and note them in field_header_vals.txt
CRVAL3 = 1.414366893369E+09
CRPIX3 = 1
CDELT3 = 1.220703125000E+04

#helpful values
HI_restfreq = 1420405751.77 * u.Hz
range_hist=[1.4177e9, 1.422e9] #analysis range (rough [wide] range where Galactic emission will be)

#adapted from https://github.com/kmhess/SoFiA-image-pipeline/blob/master/src/modules/functions.py
def chan2freq(channels):
    #converting from channels to frequencies
    frequencies = (CDELT3 * channels + CRVAL3) * u.Hz 
    return frequencies

#adapted from https://github.com/kmhess/SoFiA-image-pipeline/blob/master/src/modules/functions.py
def freq2chan(frequencies):
    #converting from frequencies to channels
    channels = (((frequencies) - CRVAL3)/CDELT3)

    #transforming into right format
    if hasattr(channels, "__len__"):
        channels = np.array([int(channel) for channel in channels])
    else:
        channels  = int(channels)
    return channels

#adapted from https://github.com/kmhess/SoFiA-image-pipeline/blob/master/src/modules/functions.py
def freq2vel(frequencies):
    redshift = (HI_restfreq-(frequencies*u.Hz))/(frequencies*u.Hz)
    velocities = (redshift*const.c).to(u.km/u.s)
    return velocities.value

#adapted from https://github.com/kmhess/SoFiA-image-pipeline/blob/master/src/modules/functions.py
def vel2freq(velocities):
    redshift = (velocities*(u.km/u.s))/(const.c)
    frequencies = HI_restfreq/(1+redshift)
    return frequencies.value

#plots RMS noise as a function of channel and channels with high RMS in the Galactic emission region to exclude
def plot_rms_channel(loc_dir, field_name, splinefits, rms_table, chan_range, high_rms_mask):
    #plotting RMS vs. frequency/velocity
    fig, ax = plt.subplots()

    ax.scatter(rms_table['Frequency'][rms_table['RMS']<3], rms_table['RMS'][rms_table['RMS']<3], s=4)

    #plotting range of channels used in the high-rms analysis/calculation
    ax.vlines(range_hist[0], 0.975,1.1, color='black', 
               linestyles='--', label='Analysis Range')
    ax.vlines(range_hist[1], 0.975,1.1, color='black', linestyles='--')

    #warning user if channel range is large or 0
    ex_range_width = abs(freq2vel(chan2freq(chan_range[1]).value)-freq2vel(chan2freq(chan_range[0]).value))
    
    if ex_range_width>150:
        print('WARNING: The calculated Galactic velocity range is large (>150 km/s). Please check to make sure the exclusion channel range is correct.')
    elif ex_range_width==0:
        print('No channels were excluded.')

    #checking if there are any high-rms channels 
    if (ex_range_width!=0):
        #plotting calculated range of channels to exclude
        ax.scatter(rms_table['Frequency'][high_rms_mask], rms_table['RMS'][high_rms_mask], s=4, color='red')
        ax.vlines(chan2freq(chan_range[0]).value, 0.975,1.1, color='green', 
                   linestyles='--', label='Exclusion Range', alpha=0.3)
        ax.vlines(chan2freq(chan_range[1]).value, 0.975,1.1, color='green', linestyles='--', alpha=0.3)
    
    #plotting expected velocity ranges from Galactic emission model
    l, b = get_lb(splinefits)
    exp_vel = get_gal_vel(l, b)
    ax.fill_betweenx([0.975,1.1],
                     vel2freq(exp_vel[0]), vel2freq(exp_vel[1]), color='pink', 
               label='Expected Range', alpha=1.0, zorder=0)

    ax.legend(loc='upper right')
    ax.set_xlabel('Freq [Hz]')
    ax.set_ylabel('RMS [Jy/beam]')
    ax.set_ylim(0.975,1.1)
    secax = ax.secondary_xaxis('top', functions=(freq2vel, vel2freq))
    secax.set_xlabel('Velocity [km/s]')
    
    plt.title('Field {field}: Channels {min_chan} to {max_chan} excluded'.format(field=field_name, min_chan=chan_range[0], 
                                                                                max_chan=chan_range[1]))
    plt.savefig(loc_dir + '{field}_rms_per_channel.png'.format(field=field_name), bbox_inches='tight')


#finds and plots RMS noise as a function of channel in the cube, along with selected channels with high RMS in the Galactic emission region to exclude
def find_rms_range(loc_dir=None, field_name=None, splinefits=None):
    #checking if rms is on file. If not, calculates and saves it.
    rms_path = loc_dir + '{field_name}_RMS_mean.txt'.format(field_name = field_name)
    if os.path.exists(rms_path):
        print('RMS file found! Reading it in.')
        rms_table = ascii.read(rms_path)
        rms = rms_table['RMS']
    else: 
        print('RMS file not found. Calculating and saving RMS.')
        #reading in cube
        filtered_cube = SpectralCube.read(splinefits)
        filtered_cube.allow_huge_operations=True
        #calculating rms and mean
        rms = np.std(filtered_cube, axis=(1, 2))
        mean = np.mean(filtered_cube, axis=(1, 2))
        
        #saving rms to file 
        rms_table = Table([np.array(filtered_cube.spectral_axis), rms, mean], names=['Frequency', 'RMS', 'Mean'])
        rms_table.write(rms_path, format='ascii',overwrite=True)

    freq_and_reasonable_rms_mask = ((rms_table['Frequency']<range_hist[1])&
                                 (rms_table['Frequency']>range_hist[0])&
                                 (rms_table['RMS']<3))
    rms_mean1 = np.nanmean(rms_table['RMS'][freq_and_reasonable_rms_mask])
    rms_std1 = np.nanstd(rms_table['RMS'][freq_and_reasonable_rms_mask])
    
    #calculating standard deviation with 3-sigma outlier detection
    new_mask = (freq_and_reasonable_rms_mask & 
                (rms_table['RMS']>(rms_mean1-(3*rms_std1)))&
                (rms_table['RMS']<(rms_mean1+(3*rms_std1))))
    rms_mean = np.nanmean(rms_table['RMS'][new_mask])
    rms_std = np.nanstd(rms_table['RMS'][new_mask])

    #getting range of frequencies with high rms, choosing min/max channels as the +/- 5 channels around it
    high_rms = rms_mean+4*rms_std#1.0+1.5*rms_std
    high_rms_mask = ((rms_table['RMS']<3)&
                     (rms_table['RMS']>high_rms)&
                     #(rms_table['Mean']>0)&
                    (rms_table['Frequency']<range_hist[1])&
                    (rms_table['Frequency']>range_hist[0]))
    chan_offset= 5

    #checking if there are any high-rms channels 
    if len(rms_table[high_rms_mask])>0:
        #converting edges of high-rms frequency range into channels 
        max_chan_edges = freq2chan(np.array([min(rms_table['Frequency'][high_rms_mask]),
                                         max(rms_table['Frequency'][high_rms_mask])]))
        chan_range = (min(max_chan_edges)-chan_offset,max(max_chan_edges)+chan_offset)
    else:
        chan_range = (0,0)

    #calculating expected velocity ranges from Galactic emission model
    l, b = get_lb(splinefits)
    exp_vel = get_gal_vel(l, b)

    #updating channel range based on expected ranges
    exp_chan_range = freq2chan(vel2freq(exp_vel))

    final_chan_range = (min(np.append(chan_range,exp_chan_range)),
                        max(np.append(chan_range,exp_chan_range)))

    #plotting
    plot_rms_channel(loc_dir, field_name, splinefits, rms_table, final_chan_range, high_rms_mask)
        
    #returning endpoints of the range of high RMS channels
    return int(final_chan_range[0]), int(final_chan_range[1])
