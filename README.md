# smbanalyze

Load the package for interactive use:
```python
    from matplotlib.pyplot import *
    from smbanalyze.shell import *
```
Gives you immediate access to most packages plus some convenience functions. This is the recommended way of using the package for interactive use.

See the notebooks/ for (somewhat hacky) example code.

experiment.py and datatypes.py define the data types experiment.Pulling, TrapData, and FretData (these last two may by simplified in the future.)

fec.py fits WLC to Pulling experiments and extracts features.

refolding.py is used for calculating binding curves from refolding time experiments.

db/ is used for storing/retrieving experiments saved in databases. Currently, this only houses mongoDB functions, but in the future, it can be used to create a local SQLlite database for potentially more robust data crunching (if needed.)

fplot.py has the plotting functions for pretty graphs.

fcalc.py calculates the fret from the image data.

## New features

Currently under active development (i.e. expect lots of bugs and caveats if you don't know exactly what's going on) is an interface to a MongoDB instance using pymongo.
```python
datadb = db.connect()
db.find(datadb.tpp, construct='wt2at')
```
## Experiment Workflow

### Example

First, process the .img files into .fret files for loading. Then load them up, align them, and pick out the different bound states.
```python
construct = 'WT2at'
slide_id,mol_id = 2,5
mol_name = 's%dm%d' % (slide_id,mol_id)

beta,gamma = 0.13, 1.1
fcalc.processMatch(construct, mol_name, 'refold', 'up', 
background='background.img', 
beta=beta, gamma=gamma,
roi='roi.txt')

exp = experiment.fromMatch(construct, mol_name, 'refold', 'up')

# Add extra metadata from conditions
metadata = {'conditions': {
        'tpp': 10e-6,
        'mg': 4e-3,
        'tx': 3e-3,
        'tq': 0.19,
        'edta': 0.1e-3},
    'construct': construct.lower(),
    'mol': mol_id,
    'slide': slide_id
    }
for p in exp:
    p.metadata.update(metadata)

exp.plot()

# Align traces still a little buggy!!

exp.plot('-', show_fret=False)

align_to=pick_pt()
align_cutoff,f_cutoff = align_to
f_range = (f_cutoff,f_cutoff+1)

to_align = to_align.has_value(trap_ext_atleast=align_cutoff)
to_align.adjustOffset(force_x_range=(720,800), ext_f_range=f_range, ext_x_range=align_cutoff)

# Split off strongly bound from the rest
# Pick ONE POINT on the FEC where strongly bound go ABOVE
# and weakly bound go BELOW
strong_pt = pick_pts(1)
weak_unbound, strong = fec.split(exp, *strong_pt)

# split off weakly bound from unbound
# Here, I'm using TWO POINTS: pick points where unbound go
# BELOW BOTH POINTS, leaving weakly bound to go above atleast ONE POINT
weak_pt = pick_pts(2)
unbound, weak = fec.split(weak_unbound, *weak_pt)

# Plot them
fig('TPP binding')
fplot.fplotall(trap=unbound, style='k-')
fplot.fplotall(trap=weak, style='y-')
fplot.fplotall(trap=strong, style='r-')

# Add there binding state to their metadata
def set_all(list_, key, value):
    for x in list_:
        x[key] = value
set_all(weak, 'bind_tpp', 1)
set_all(strong, 'bind_tpp', 2)
set_all(unbound, 'bind_tpp', 0)

# More to come!!
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
