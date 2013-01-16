# smbanalyze

## Usage

Load the package:

    from smbanalyze import *

Loading a single image:

```python
image = Image.fromFile('filename.img')
background = Image.fromFile('bg_filename.img', background=True)

image_bg_subtracted = Image.fromFile('filename.img', background='bg_filename.img')

# is true!
image - background == image_bg_subtracted
'''
