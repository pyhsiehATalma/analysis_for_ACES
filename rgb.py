# RGB images
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import (AsinhStretch, PercentileInterval, HistEqStretch, SinhStretch, SqrtStretch, PowerStretch)

# Load the FITS files
data1 = fits.getdata('CS21_CubeMosaic_max.fits')  # Red  
data2 = fits.getdata('SO32_CubeMosaic_max.fits')  # Green  
#data3 = fits.getdata('CH3CHO_CubeMosaic_max.fits')  # Blue  
data3 = fits.getdata('CH3CHO_CubeMosaic_edgelessmax.fits')  # "Blue" 

# Normalize using percentile interval

# Stretch & Normalize Data

#LinearStretch()	No transformation, just normal scaling	Data with uniform distribution
#PowerStretch(a)	Power-law transformation, enhances bright areas more	Highlighting bright features
#LogStretch()	Logarithmic scaling, expands dark regions more	Data with large dynamic range
#AsinhStretch()	Similar to log, but smoother in dark regions	Preserving details in both bright & dark areas
#SinhStretch()	Opposite of AsinhStretch, expands bright regions more	Data with weak background noise
#SqrtStretch()	Square-root scaling, enhances contrast	Moderate dynamic range
#HistEqStretch()	Equalizes histogram for uniform intensity distribution	Enhancing faint details

interval1 = PercentileInterval(97.5)
interval2 = PercentileInterval(99.5)
interval3 = PercentileInterval(99.5)

stretch1 = AsinhStretch()
stretch2 = AsinhStretch()
stretch3 = AsinhStretch()
#stretch3 = PowerStretch(0.3)

data1_stretched = stretch1(interval1(data1))
data2_stretched = stretch2(interval2(data2))
data3_stretched = stretch3(interval3(data3))

# Stack Data into a Custom 3 -Color Image
rgb_image = np.stack([data1_stretched, data2_stretched, data3_stretched], axis=-1)

# Create a Single Plot for the False-Color Composite
fig, ax = plt.subplots(figsize=(20, 10))
ax.imshow(rgb_image, origin='lower')
ax.set_title('False Color (R:CS21, G: SO32, B: CH3CHO)')
ax.axis('off')

# Save and Show
plt.savefig("FalseColor.png", dpi=300, bbox_inches='tight')
plt.show()
