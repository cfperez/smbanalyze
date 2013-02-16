# smbanalyze

Load the package:

    from matplotlib.pyplot import *
    from smbanalyze import *

## Advanced Usage

Process images using give roi file and background in the current directory:

```python
img_files = FileIO.flist('SJF', 'img')
fret_data = FRET.processFiles(img_files, roi='roi.txt', background='SJF_background.img')
```
Remember you can change directories with `cd` and list contents with `ls`.

Load a pulling experiment:

```python
pull = Experiment.Pulling.fromFile('SJF4_s1m1_4')
pull.plot()

# or from multiple files
pulls = Experiment.fromMatch('SJF4', 's1m1')
figure()
pulls[0].plot(FEC=True)
```

Fit to the section before the rip and automatically plot if the figure is open:
```python
# x is the MINIMUM extension and f is the MAXIMUM force to fit to
# can also use x=(min,max) and f=(min,max)
fit = pulls[0].fitForceExtension(x=750, f=9)
# can also get the fit and parameters through
pulls[0].fit.parameters
pulls[0].fit['Lc']
```
Note that the default is the regular Marko Siggia curve. But Hallelujah! I figured out what
was going on with the MMS curve fit, and made my own (better) version
```python
fitMMS = pulls[0].fitForceExtension(x=750, f=9, fitfunc='MMS')
```
Note that this fit is more testy (e.g. K parameter is auto-fixed to 1200 by default).

## Basic Usage

Loading a single image:

```python
image = Image.fromFile('test_s1m1.img')
background = Image.fromFile('background.img', background=True)

image_bg = Image.fromFile('test_s1m1.img', background='background.img')

# is true!
image - background == image_bg

# Properties
image.frames
image.width
image.height
image.times

# display frames -- starts with 0
image_bg.show(0)
image_bg.show(10)
```

Get out the donor/acceptor counts:

```python
ROIs = Image.ROI.fromFile('roi1.txt')
image_bg.addROI(*ROIs)

# or bottom,left and top,right and origin
donorROI = Image.ROI( (5,5), (20,25), origin='relative', name='donor' )
# left, right, bottom, top
acceptorROI = Image.ROI.fromCorners( 50, 60, 90, 100, origin='absolute', name='acceptor' )
image_bg.addROI(donorROI,acceptorROI)

donor,acceptor = image_bg.donor,image_bg.acceptor
```

Origin reflects whether the pixels are numbered with respect to the CCD origin ('absolute') or the image origin ('relative'). Absolute ROIs are robust to changes in the subimage coordinates.

To calculate fret:

```python
fret = FRET.calculate(image_bg) # optional: beta= , gamma=

# Access as named fields
fret.time, fret.donor, fret.acceptor, fret.fret

# or as a tuple (order matters!)
T,donor,acceptor,f = fret

# and save to a file
FRET.toFile('saveTo.fret', fret)
```

`fret` is a special "named tuple" from the collections package in the python library with more flexible usage, as shown above. Don't be confused; it's just a tuple in which each position also has a name which you can see with `fret._fields`.

Plot the data:

```python
FRET.plot(fret, title='test1')

# Rerunning this command overwrites the current figure
# but if you want a new figure
figure()
FRET.plot(fret)

# and if you have pulling data
pull = FileIO.loadstr('test_s1m1.str')
FRET.plot(fret, pull=pull)

# and if you want to see an FEC
FRET.plot(pull, FEC=True)
# this works too, though they don't align
FRET.plot(fret, pull, FEC=True)
```
