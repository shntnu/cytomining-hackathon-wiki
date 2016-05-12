import numpy as np
from keras.utils import np_utils

def random_crop(img, size=192):
    _, ny, nx = img.shape
    i = np.random.randint(0, ny - size)
    j = np.random.randint(0, nx - size)
    return img[:,
               i:(i + size),
               j:(j + size)]


def load_cropped_site(grp, size=192):
    return random_crop(np.array([grp[dset].value for dset in grp.keys()]),
                       size)


def load_image(grp, size=192):
    sites = grp.keys()
    site = np.random.choice(sites, 1)[0]
    return load_cropped_site(grp[site], size)


def load_minibatch(datastore, frame, target_mapper, size=192):
    y = np_utils.to_categorical(map(lambda x: target_mapper[x], frame['compound']), len(target_mapper))
    x = np.array([load_image(datastore[drug][conc][rep], size=size)
                  for drug, conc, rep in zip(frame['compound'], frame.concentration, frame.replicate)])
    return x.astype('float32') / 65535., y


def minibatch_generator(datastore, frame, target_mapper, size=192):
    while 1:
        yield load_minibatch(datastore, frame, target_mapper, size=192)
