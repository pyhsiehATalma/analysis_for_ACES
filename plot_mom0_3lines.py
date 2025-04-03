import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import AsinhStretch, PercentileInterval
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord

import matplotlib.colors as mcolors
import pylab as pl

# Set font, I am used to the IDL style plotting, so the fonts I choose are, Arial, Helvetica, or Nimbus Sans.

#plt.rcParams["font.family"] = "Nimbus Sans"
plt.rcParams["font.family"] = "Arial"

# Create the custom colormap by stacking two colormaps
colors1 = pl.cm.gray_r(np.linspace(0., 1, 128))  # Reversed gray colormap
colors2 = pl.cm.hot(np.linspace(0, 1, 128))  # Hot colormap

colors = np.vstack((colors1, colors2))  # Stack both colormaps
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)  # Create colormap

# Define FITS files and their corresponding PercentileIntervals
# fits_files=[File name, line (for title), min_PercentileIntervals, max_PercentileIntervals], you can change the mix/max.

fits_files = [
    ('CS21_CubeMosaic_masked_hlsig_dilated_mom0.fits',  'CS(2-1)', 1, 99.9),
    ('SO32_CubeMosaic_masked_hlsig_dilated_mom0.fits', 'SO 3(2)-2(1)', 1, 99.9),
    ('CH3CHO_CubeMosaic_masked_hlsig_dilated_mom0.fits','CH3CHO 5(3)-4(3)', 1, 99.0)
]

# Initialize figure with 3 rows and 1 column
fig = plt.figure(figsize=(12, 18))
fig.subplots_adjust(hspace=0.015)  # Reduce vertical spacing between rows

# Loop through FITS files and plot each with different PercentileInterval
#for i, (fits_file, interval, title) in enumerate(fits_files):
for i, (fits_file, title, p_min, p_max) in enumerate(fits_files):
    with fits.open(fits_file) as hdulist:
        data = hdulist[0].data
        wcs = WCS(hdulist[0].header)

    # Create subplot with WCS projection
    ax = fig.add_subplot(3, 1, i+1, projection=wcs)

    # Plot the FITS data
# Compute vmin and vmax using specified percentiles

    vmin = np.nanpercentile(data, p_min)   # Use the specified min percentile
    vmax = np.nanpercentile(data, p_max)   # Use the specified max percentile

    im = ax.imshow(data, origin='lower', cmap=mymap, vmin=vmin, vmax=vmax)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.01725, pad=0.025)
    cbar.set_label(r"Intensity (Jy km s$^{-1}$ beam$^{-1}$)", fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    
    # Set axis labels and title
    ax.set_xlabel('Galactic Longitude', fontsize=16)
    ax.set_ylabel('Galactic Latitude', fontsize=16)
    ax.set_title(f"{title}", fontsize=16)

    # Tick labels formatting (fixing WCS support)
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.coords[0].set_ticklabel(size=16)
    ax.coords[1].set_ticklabel(size=16)

    # Set tick spacing
    ax.coords[0].set_ticks(spacing=0.2 * u.deg)
    ax.coords[1].set_ticks(spacing=0.2 * u.deg)

    ax.coords[0].set_ticks(number=10)
    ax.coords[1].set_ticks(number=5)

    ax.coords[0].display_minor_ticks(True)
    ax.coords[1].display_minor_ticks(True)

    ax.coords[0].set_minor_frequency(4)
    ax.coords[1].set_minor_frequency(2)

    ax.coords[0].tick_params(which='minor', length=5)
    ax.coords[1].tick_params(which='minor', length=5)

    # Tick directions inward
    ax.coords[0].tick_params(which='both', direction='in')
    ax.coords[1].tick_params(which='both', direction='in')

    # Galactic coordinates for labels (Galactic Longitude, Galactic Latitude)
    label_coords = [
        (0.679, 0.13, "Sgr B"),  
        (-0.490, 0.13, "Sgr C"),  
        (359.925, 0.22, "Polar Arc"),
    ]
    for lon, lat, label in label_coords:
        sky_coord = SkyCoord(lon, lat, frame='galactic', unit='deg')
        x_pixel, y_pixel = wcs.world_to_pixel(sky_coord)

        # Add label at the corresponding WCS location
        ax.text(x_pixel, y_pixel, label, color='white', fontsize=17, ha='center', va='center', bbox=dict(facecolor='black', alpha=0.2, edgecolor='none'))

    # arrow_coords=[start_x, start_y (text position), end_x, end_y (source position)]
    arrow_coords = [
        (359.8, 0.1, 359.944, -0.046, "CND and Sgr A*"),  
        (359.8, -0.2, 359.872, -0.087, "20 MC"),
        (359.95, -0.2, 359.978, -0.087, "50 MC"),
        (0.47, -0.177, 0.16, -0.06, "The Quintuplet Cluster"),
        (0.25, 0.125, 0.12, 0.02, "The Arches Cluster"),
        (0.5, 0.18, 0.253, 0.016, "The Brick\n(G0.253+0.016)"),
        (0.3, -0.225, 0.145, -0.086, "Straw Cloud\n(G0.145-0.086)"),
        (0.2, -0.3, 0.106, -0.082, "Sticks Cloud\n(G0.106-0.082)"),
        (0.05, -0.25, 0.068, -0.075, "Stone Cloud\n(G0.068-0.075)"),
    ]

    for lon_start, lat_start, lon_end, lat_end, label in arrow_coords:
        # Convert Galactic coordinates to pixel coordinates
        sky_coord_start = SkyCoord(lon_start, lat_start, frame='galactic', unit='deg')
        sky_coord_end = SkyCoord(lon_end, lat_end, frame='galactic', unit='deg')
        x_start, y_start = wcs.world_to_pixel(sky_coord_start)
        x_end, y_end = wcs.world_to_pixel(sky_coord_end)

        # Draw the arrow with a label
        ax.annotate(label, xy=(x_end, y_end), xytext=(x_start, y_start),
                    arrowprops=dict(arrowstyle="-", linestyle="--",color="black", lw=1),
                    fontsize=10, color="black")


# Save and show the plot
plt.savefig("fits_3row_plot.png", dpi=300, bbox_inches='tight')
#plt.show()
