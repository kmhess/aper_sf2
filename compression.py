from glob import glob
import os

weights_cubes = glob('/mnt/data/mos_*/*weights.fits')
noise_cubes = glob('/mnt/data/mos_*/*noise.fits')
mask_cubes = glob('/mnt/data/mos_*/*mask.fits')

for w in weights_cubes:
    if not os.path.isfile(w + '.fz') & os.path.isfile(w):
        os.system('fpack {}'.format(w))
        print('fpack {}'.format(w))
    if os.path.isfile(w + '.fz') & os.path.isfile(w):
        os.system('rm -r {}'.format(w))
        print('rm -r {}'.format(w))

for n in noise_cubes:
    if not os.path.isfile(n + '.fz') & os.path.isfile(n):
        os.system('fpack {}'.format(n))
        print('fpack {}'.format(n))
    if os.path.isfile(n + '.fz') & os.path.isfile(n):
        os.system('rm -r {}'.format(n))
        print('rm -r {}'.format(n))

for m in mask_cubes:
    if not os.path.isfile(m + '.fz') & os.path.isfile(m):
        os.system('fpack {}'.format(m))
        print('fpack {}'.format(m))
    if os.path.isfile(m + '.fz') & os.path.isfile(m):
        os.system('rm -r {}'.format(m))
        print('rm -r {}'.format(m))
