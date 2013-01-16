# smbanalyze

## Usage

Load the package:

    from smbanalyze import *

Loading a single image:

```python
image = Image.fromFile('filename.img')
background = Image.fromFile('bg_filename.img', background=True)

image_bg = Image.fromFile('filename.img', background='bg_filename.img')

# is true!
image - background == image_bg
```

Get out the donor/acceptor counts:

```python
ROIs = Image.ROI.fromFile('roi.txt')
image_bg.addROI(*ROIs)

# or bottom,left and top,right and origin
donorROI = Image.ROI( (5,5), (20,25), origin='relative', name='donor' )
# left, right, bottom, top
acceptorROI = Image.ROI.fromCorners( 50, 60, 90, 100, origin='absolute', name='acceptor' )
image_bg.addROI(donorROI,acceptorROI)

donor,acceptor = image_bg.donor,image_bg.acceptor
```
