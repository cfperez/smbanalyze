from smbanalyze import Image, FileIO
import matplotlib.pyplot as plt

############################################
## CHANGE filename to what you want
## then run this file
############################################
filename = 'test.img'

# Create a Stack from a file
image = Image.Stack(filename)

# Create a named ROI (default is 'donor' => maybe defaults are bad?)
donor = Image.ROI( (28,258), (39,267), name='donor' )

# Add the named ROI to the Stack for later processing
image.addROI(donor)

# OR can name ROI when adding to Stack
# roi=Image.ROI( (28,258), (39,267) )
# image.addROI(roi, 'donor')

# This uses absolute coordinates; default is 'relative' coordinates i.e. wrt subimage
acceptor = Image.ROI( 95,105, 85, 96, origin='absolute', name='acceptor')
image.addROI(acceptor)

# To avoid setting the same ROIs every time you load a new img file, set the default ROIs 
# to use for every new Stack. Add them in any order. ROI with same name are overwritten.
Image.setDefaultROI(donor, acceptor)

# Pull out the subimages for each ROI. Because they are named, can use multiple ROIs.
# this returns another Stack object, meaning you can look at the donor/acceptor subimages
# by using e.g. donor[0].show() (donor.show() with slider to be implemented)
donor = image['donor']

# Plot here. You can use .counts on ANY Stack object, but generally only makes sense
# on subregions returned using the ROI API
plt.clf()
plt.subplot(211)
plt.plot( donor.counts, 'b-', label='donor' )

# or if using default names for ROIs (which can be changed), then simply
plt.plot( image.acceptor, 'r-' , label='acceptor')
plt.legend()
# which is shorthand for image['donor'].counts

# calculate FRET using properly ROI'ed Image.Stack => this function uses .donor
# and .acceptor, so requires following proper conventions
fret = Image.calcFRET(image)

plt.subplot(212)
plt.plot( fret, 'g-', label='FRET' )

# show the plots!
plt.title('FRET')
plt.show()

# You can use FileIO.savedat() to save to text
FileIO.savedat('testing.txt', (donor.counts,image.acceptor), header='donor acceptor', fmt='%u')

# or use this shorthand to also calculate and save FRET data
Image.saveFRETdata('testing2.txt', image)
