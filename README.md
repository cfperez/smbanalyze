# smbanalyze

Load the package:

    from matplotlib.pyplot import *
    from smbanalyze import *

## New Features

Save or load an experiment:
```python
pull = experiment.Pulling.load('filename.exp')
pull.save('filename2.exp')
```
## Experiment Workflow

### Tips

1. Use generically named yet descriptive variables to allow for some automation with IPython macros.

2. The new experiment.List collection object has some useful features demonstrated below. We will continue to add features to make manipulating multiple experiments easier. Requests welcome!

### Example

```
pulls = experiment.fromMatch('SJ2UL')
mol = pulls.matching('s1m1')

# if you need to remove an experiment
del mol[0]
# or add an experiment
mol += [pulls[0]]

mol.adjustForceOffset()
# soon => mol.adjustExtOffset()
# or for both => mol.adjustOffset()

figure()
mol.plotall('pull')

# individual plots--there may be alot!
mol.plot()

mol.fitHandles(830, 10)
mol.fitRip(1000, 15)
```

```python
pulls = experiment.fromMatch('SJ2UL')
mol = pulls.matching('s1m1')

hstart = 850  # extension to start fit
hend = 9      # force to stop fit
p = mol.next() # Imagine this is line [3] in IPython

# For now, simply print this force offset
mean( p.fec.force[10:20] )
p.fitHandles(hstart, hend)
p.plot()

# make a macro of previous
macro fitp 4-6
macro nextp 3

rstart = 1000
rend = 15
p.fitRip(rstart, rend) # line [11]

# Optional. Now just adjust rstart as needed before running this macro
macro fitrip 11

nextp
fitp
rstart = 1010
fitrip
```

## Advanced Usage

Process images using give roi file and background in the current directory:

```python
fcalc.processMatch('SJF4', roi='roi.txt', background='SJF_background.img')
```
Remember, in IPython, you can change directories with `cd` and list contents with `ls`.

Load a pulling data:

```python
pull = experiment.fromFile('SJF4_s1m1_4')
pull.plot()

# or from multiple files
pulls = experiment.fromMatch('SJF4', 's1m1')
for a_pull in pulls:
  figure()
  a_pull.plot(FEC=True)
```

Fit to the section before the rip and automatically plot if the figure is open:
```python
# x is the MINIMUM extension and f is the MAXIMUM force to fit to
# can also use x=(min,max) and f=(min,max)
fit = pull.fitHandles(x=750, f=9)
# can also get the fit and parameters through
pull.handles.parameters
pull.handles['Lc']
```

Fit the upper portion to get the rip size:
```python
pull.fitRip(f=(10,20))
```

## Basic Usage

Loading a single image:

```python
img = image.fromFile('test_s1m1.img')
background = image.fromFile('background.img', background=True)

image_bg = image.fromFile('test_s1m1.img', background='background.img')

# is true!
image - background == image_bg

# Properties
img.frames
img.width
img.height
img.times

# display frames -- starts with 0
image_bg.show(0)
image_bg.show(10)
```

Get out the donor/acceptor counts:

```python
ROIs = image.ROI.fromFile('roi1.txt')
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
fret = fcalc.calculate(image_bg) # optional: beta= , gamma=

# Access as named fields
fret.time, fret.donor, fret.acceptor, fret.fret

# or as a tuple (order matters!)
T,donor,acceptor,f = fret

# and save to a file
fcalc.toFile('saveTo.fret', fret)
```

`fret` is a special "named tuple" from the collections package in the python library with more flexible usage, as shown above. Don't be confused; it's just a tuple in which each position also has a name which you can see with `fret._fields`.

Plot the data:

```python
fplot.plot(fret, title='test1')

# Rerunning this command overwrites the current figure
# but if you want a new figure
figure()
fret.plot(fret)

# and if you have pulling data
pull = TrapData.fromFile('test_s1m1.str')
fplot.plot(fret, pull=pull)

# and if you want to see an FEC
fplot.plot(pull, FEC=True)
# this works too, though they don't align
fplot.plot(fret, pull, FEC=True)
```
