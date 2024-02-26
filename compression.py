from glob import glob
import os

weights_cubes = glob('/mnt/data/mos_*/*weights.fits')
mask_cubes = glob('/mnt/data/mos_*/*mask.fits')

print(weights_cubes[:20])
print(mask_cubes[:20])

for w in weights_cubes[:5]:
    if not os.path.isfile(w + '.fz') & os.path.isfile(w):
        os.system('fpack {}'.format(w))
        print('fpack {}'.format(w))
    if os.path.isfile(w + '.fz') & os.path.isfile(w):
        os.system('rm -r {}'.format(w))
        print('rm -r {}'.format(w))

for m in mask_cubes[:5]:
    if not os.path.isfile(w + '.fz') & os.path.isfile(m):
        os.system('fpack {}'.format(m))
        print('fpack {}'.format(m))
    if os.path.isfile(m + '.fz') & os.path.isfile(m):
        os.system('rm -r {}'.format(m))
        print('rm -r {}'.format(m))        