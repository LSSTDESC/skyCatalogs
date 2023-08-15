import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sncosmo


from skycatalogs.skyCatalogs import SkyCatalog, open_catalog
from skycatalogs.objects.base_object import BaseObject, load_lsst_bandpasses
from skycatalogs.utils.sn_tools import SncosmoModel

PIXEL = 9556

def explore_lc(obj):
    if obj.object_type != 'sn':
        print('No light curve for object of type ', obj.object_type)
        return

    for native in ['id', 'ra', 'dec', 'salt2_params']:
        print(native,'=',obj.get_native_attribute(native))

    # get fluxes for some days around t0
    params = obj.get_native_attribute('salt2_params')
    t0 = params['t0']


    fluxes_mjd = dict()
    dt_rng = '-5_30_3'
    #for dt in np.arange(0, 40, 2):
    fluxes_mjd[-50] = obj.get_LSST_fluxes(cache=False, mjd=(t0 - 50))
    for dt in np.arange(-5, 30, 3):
        fluxes_mjd[dt] = obj.get_LSST_fluxes(cache=False, mjd=(t0 + dt))

    # Also add far-out points which could be outside time interval sn
    # is active

    fluxes_mjd[50] = obj.get_LSST_fluxes(cache=False, mjd=(t0 + 50))

    plt.figure()
    by_band =  dict()
    for b in 'ugrizy':
        print('fluxes for band ', b)
        band_fluxes = []
        for dt in fluxes_mjd.keys():
            print(f'dt: {dt}  flux: {fluxes_mjd[dt][b]}')
            band_fluxes.append(fluxes_mjd[dt][b])
        # Also plot flux vs dt
        plt.plot(list(fluxes_mjd.keys()), band_fluxes, label=f'{b}')
    plt.legend(fontsize='x-small', ncol=2)
    plt.title(f'{obj.id}')
    plt.xlabel('dt (days)')
    plt.ylabel('flux')
    # _nm suffix probably means "no magnorm"
    plt.savefig(f'{obj.id}_dt_fluxes_nm.png')

    plt.close()

    # also plot SEDs
    sn_obj = SncosmoModel(params=params)


    plt.figure()
    for dt in np.arange(-5, 30, 3):
        sed = sn_obj.get_sed(t0 + dt)
        plt.plot(sed.wave_list, sed(sed.wave_list), label=f'{dt}')
    plt.yscale('log')
    ##plt.ylim(3e-8, 2e-3)
    plt.legend(fontsize='x-small', ncol=2)
    plt.xlabel('wavelength (nm)')
    plt.ylabel('photons/nm/cm^2/s')
    plt.savefig(f'{obj.id}_seds.png')
    plt.close()


def make_sncosmo_lc(obj):
    params = obj.get_native_attribute('salt2_params')

    # Set up the sncosmo source
    src = sncosmo.Model(source='salt2-extended')
    src.set(**params)

    bandpasses = load_lsst_bandpasses()
    sncosmo_bandpasses = []
    for nm,val in bandpasses.items():
        snc_bp = sncosmo.Bandpass(val.wave_list, [val(wv) for wv in val.wave_list], name=nm, wave_unit='nm')
        sncosmo_bandpasses.append(snc_bp)

    dt_start = -5
    dt_end  = 30
    dt_incr = 3
    t0 = params['t0']

    dt_rng = '-5_30_3'
    plt.figure()
    #for b_name, b in bandpasses.items():
    for b in sncosmo_bandpasses:
        fluxes = []
        rel_times = []
        for dt in np.arange(dt_start, dt_end, dt_incr):
            fluxes.append(src.bandflux(b, t0 + dt))
            rel_times.append(dt)
        plt.plot(rel_times, fluxes, label=f'{b.name}')
    plt.legend(fontsize='x-small', ncol=2)
    plt.xlabel('dt (days)')
    plt.ylabel('flux')
    plt.savefig(f'sn_cosmo_{obj.id}_dt{dt_rng}_fluxes.png')


def explore(cat, obj_type, ix_list=[0]):
    obj_list = cat.get_object_type_by_hp(PIXEL, obj_type)

    collects = obj_list.get_collections()
    print('Number of collections: ', len(collects))

    icoll = 0
    for c in collects:
        print(f'Object type: {c._object_type_unique}')
        print('Native columns:')
        print(c.native_columns)

        len_coll = len(c)
        obj0 = c[0]
        print(f"\nFor hpid {c.get_partition_id()} collection {icoll} found {len_coll} objects")
        print("First object: ")
        print('id=', obj0.id, ' ra=', obj0.ra, ' dec=', obj0.dec,
              ' belongs_index=', obj0._belongs_index)

        all_cols = obj0.native_columns
        if obj0.object_type == 'sn':
            extras = {'lsst_flux_u', 'lsst_flux_g', 'lsst_flux_r',
                      'lsst_flux_i', 'lsst_flux_z', 'lsst_flux_y', 'mjd'}
            all_cols.difference_update(extras)

        for native in all_cols:           # obj0.native_columns:
            print(native,'=',obj0.get_native_attribute(native))
        icoll += 1

        if obj0.object_type == 'sn':
            for ix in ix_list:
                print('sn object index is ', ix)
                explore_lc(c[ix])
                ##make_sncosmo_lc(c[ix])
            # # get fluxes for some days around t0
            # params = obj0.get_native_attribute('salt2_params')
            # t0 = params['t0']
            # fluxes = dict()
            # for dt in np.arange(-5, 30, 3):
            #     fluxes[dt] = obj0.get_LSST_fluxes(cache=False, mjd=(t0 + dt))

            # for b in 'ugrizy':
            #     print('fluxes for band ', b)
            #     for dt in fluxes.keys():
            #         print(f'dt: {dt}  flux: {fluxes[dt][b]}')

skycatalog_root = os.path.join(Path(__file__).resolve().parents[1],
                               'skycatalogs', 'data')
config_path = os.path.join(skycatalog_root, 'ci_sample', 'skyCatalog.yaml')
#skycatalog_root = os.getenv('CFS_SKY_ROOT')
#config_path = os.path.join(skycatalog_root, 'point_test', 'skyCatalog.yaml')

cat = open_catalog(config_path, skycatalog_root=skycatalog_root)

#print('Explore star collection')
#explore(cat, {'star'})

print('explore sn collection')
###explore(cat, {'sn'}, ix_list = [0, 3, 100])
## have tried ix 3, 7, 100. 105, 10
## For ix 100 there are no visible seds. For 105 almost none
##explore(cat, {'sn'}, ix_list = [3,7,10,100, 105])
##explore(cat, {'sn'}, ix_list = [202])
explore(cat, 'star', ix_list = [10])
#print('explore both sn and star')
#explore(cat, {'sn', 'star'})
