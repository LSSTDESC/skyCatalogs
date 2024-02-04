import os
import time
import healpy
import logging # expected by make_galaxy_schema
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
# Following imports were used in old GCR_siminterface rotate class
#  from lsst.sims.utils import angularSeparation
#  from lsst.sims.utils import rotationMatrixFromVectors
#  from lsst.sims.utils import cartesianFromSpherical, sphericalFromCartesian

#  Possible alternatives are in rubin_schedule.utils.coordinate_transformations
#
from rubin_scheduler.utils import cartesian_from_spherical, spherical_from_cartesian
from rubin_scheduler.utils import rotation_matrix_from_vectors, angular_separation
from skycatalogs.skyCatalogs import open_catalog
from skycatalogs.utils.parquet_schema_utils import make_galaxy_schema
from skycatalogs.utils.creator_utils import make_MW_extinction_av
from skycatalogs.utils.shapes import Disk

DC2_RA = [68, 64, 58,
          68.5, 62, 55,
          67, 62.5, 57]
DC2_DEC = [-30.0, -32.0, -31.0,
           -34.5, -36, -36.5,
           -40, -42, -41]
FIELD_NAMES = ['COSMOS', 'DEEP_AO', 'DESI_SV3_R1',
               'Rubin_SV_095_-25', 'Rubin_SV_125_-15', 'Rubin_SV_225__-40',
               'Rubin_SV_250_2', 'Rubin_SV_280_-48', 'Rubin_SV_300_-41']
FIELD_RA = [150.1, 216.0, 179.6, 95.0, 125.0, 225.0, 250.0, 280.0, 300.0]
FIELD_DEC = [2.1819444444444445, -12.5, 0.0, -25.0, -15.0, -40.0, 2.0,
             -48.0, -41.0]
DISK_RADIUS_DG = 1.4


class FieldRotator(object):

    def __init__(self, ra0, dec0, ra1, dec1):
        """
        Parameters
        ----------
        ra0, dec0 are the coordinates of the original field
        center in degrees

        ra1, dec1 are the coordinates of the new field center
        in degrees

        The transform() method of this class operates by first
        applying a rotation that carries the original field center
        into the new field center.  Points are then transformed into
        a basis in which the unit vector defining the new field center
        is the x-axis.  A rotation about the x-axis is applied so that
        a point that was due north of the original field center is still
        due north of the field center at the new location.  Finally,
        points are transformed back into the original x,y,z bases.
        """

        # do we actually need to do the rotation, or is the simulation
        # already in the right spot?
        self._needs_to_be_rotated = True
        rot_dist = angular_separation(ra0, dec0, ra1, dec1)
        if rot_dist < 1.0/3600.0:
            self._needs_to_be_rotated = False
            return

        # find the rotation that carries the original field center
        # to the new field center
        xyz = cartesian_from_spherical(np.radians(ra0), np.radians(dec0))
        xyz1 = cartesian_from_spherical(np.radians(ra1), np.radians(dec1))
        if np.abs(1.0-np.dot(xyz, xyz1)) < 1.0e-10:
            self._transformation = np.identity(3, dtype=float)
            return

        first_rotation = rotation_matrix_from_vectors(xyz, xyz1)

        # create a basis set in which the unit vector
        # defining the new field center is the x axis
        xx = np.dot(first_rotation, xyz)
        rng = np.random.RandomState(99)
        mag = np.NaN
        while np.abs(mag) < 1.0e-20 or np.isnan(mag):
            random_vec = rng.random_sample(3)
            comp = np.dot(random_vec, xx)
            yy = random_vec - comp*xx
            mag = np.sqrt((yy**2).sum())
            yy /= mag

        zz = np.cross(xx, yy)

        to_self_bases = np.array([xx,
                                  yy,
                                  zz])

        out_of_self_bases = to_self_bases.transpose()

        # Take a point due north of the original field
        # center.  Apply first_rotation to carry it to
        # the new field.  Transform it to the [xx, yy, zz]
        # bases and find the rotation about xx that will
        # make it due north of the new field center.
        # Finally, transform back to the original bases.
        d_dec = 0.1
        north = cartesian_from_spherical(np.radians(ra0),
                                         np.radians(dec0+d_dec))

        north = np.dot(first_rotation, north)

        # print(np.degrees(sphericalFromCartesian(north)))

        north_true = cartesian_from_spherical(np.radians(ra1),
                                              np.radians(dec1+d_dec))

        north = np.dot(to_self_bases, north)
        north_true = np.dot(to_self_bases, north_true)
        north = np.array([north[1], north[2]])
        north /= np.sqrt((north**2).sum())
        north_true = np.array([north_true[1], north_true[2]])
        north_true /= np.sqrt((north_true**2).sum())

        c = north_true[0]*north[0]+north_true[1]*north[1]
        s = north[0]*north_true[1]-north[1]*north_true[0]
        norm = np.sqrt(c*c+s*s)
        c = c/norm
        s = s/norm

        # nprime = np.array([c*north[0]-s*north[1],
        #                    s*north[0]+c*north[1]])

        yz_rotation = np.array([[1.0, 0.0, 0.0],
                                [0.0, c, -s],
                                [0.0, s, c]])

        second_rotation = np.dot(out_of_self_bases,
                                 np.dot(yz_rotation,
                                        to_self_bases))

        self._transformation = np.dot(second_rotation,
                                      first_rotation)

    def transform(self, ra, dec):
        """
        ra, dec are in degrees; return the RA, Dec coordinates
        of the point about the new field center
        """
        xyz = cartesian_from_spherical(np.radians(ra),
                                       np.radians(dec)).transpose()
        xyz = np.dot(self._transformation, xyz).transpose()
        ra_out, dec_out = spherical_from_cartesian(xyz)
        return np.degrees(ra_out), np.degrees(dec_out)


class GalaxyRotator:
    def __init__(self, field_rotator, field_name, obj_list):
        self._field_rotator = field_rotator
        self._field_name = field_name
        self._obj_list = obj_list
        self._stride = 1000000

    def output_field_pixels(self, output_dir, arrow_schema):
        # First make necessary new columns and keep track of output
        # healpixels
        field_distinct = set()
        colls = self._obj_list.get_collections()
        new_ra = []
        new_dec = []
        new_av_ext = []
        hps = []
        hps_distinct = []
        # Make a dict indexed by hp number.  Value is df to be written
        # to file for that hp
        hp_out_dict = dict()
        for c in colls:
            ra, dec = self._field_rotator.transform(c._ra, c._dec)
            new_ra.append(ra)
            new_dec.append(dec)
            new_av_ext.append(make_MW_extinction_av(ra, dec))
            rot_hp = np.array(healpy.pixelfunc.ang2pix(NSIDE, ra, dec,
                                                       lonlat=True))
            hps.append(rot_hp)
            distinct = {h for h in rot_hp}
            hps_distinct.append(sorted(distinct))
            field_distinct = set.union(field_distinct, distinct)

        # hp_sorted = sorted(field_distinct)
        native = colls[0].native_columns

        for i, c in enumerate(colls):
            to_get = [n for n in native if (n != 'ra') and (n != 'dec') and (n != 'MW_av') and n != 'sed_val_bulge' and n != 'sed_val_disk']
            dat = c.get_native_attributes(to_get)
            dat['ra'] = new_ra[i]
            dat['dec'] = new_dec[i]
            dat['MW_av'] = new_av_ext[i]
            dat['hp'] = hps[i]
            # special handling for for tophat SEDs
            for k in ['sed_val_bulge', 'sed_val_disk', 'sed_val_knots']:
                dat[k] = c.get_native_attribute(k, no_np=True)

            df = pd.DataFrame(dat)
            df.sort_values('hp', inplace=True)

            # Now for hps represented in the collection, copy relevant
            # rows to df for that hp
            for hp in hps_distinct[i]:
                sub_df = (df.loc[df['hp'] == hp]).copy()
                if hp in hp_out_dict:
                    hp_out_dict[hp] = pd.concat([hp_out_dict[hp], sub_df])
                else:
                    hp_out_dict[hp] = sub_df
            del df

        # Now write out healpixels for this collection (sets of healpixels
        # belonging to different fields are disjoint)
        # for hp in hps_distinct:
        for hp in field_distinct:
            outpath = os.path.join(output_dir, f'galaxy_{hp}.parquet')
            writer = pq.ParquetWriter(outpath, arrow_schema)
            dlen = len(hp_out_dict[hp]['galaxy_id'])
            if dlen == 0:
                continue
            u_bnd = min(self._stride, dlen)
            l_bnd = 0
            while u_bnd > l_bnd:
                out_df = hp_out_dict[hp].iloc[l_bnd:u_bnd]
                out_table = pa.Table.from_pandas(out_df, schema=arrow_schema)
                writer.write_table(out_table)
                l_bnd = u_bnd
                u_bnd = min(l_bnd + self._stride, dlen)

            writer.close()
            print(f'{time.asctime()} Completed pixel {hp}')


input_config = '/pscratch/sd/j/jrbogart/desc/skycatalogs/rehearsal/skyCatalog.yaml'
rotated_dir = '/pscratch/sd/j/jrbogart/desc/skycatalogs/rotated'
cat = open_catalog(input_config)
NSIDE = 32    # should really read from config

for dc2_ra, dc2_dec, f_ra, f_dec, name in zip(DC2_RA, DC2_DEC, FIELD_RA, FIELD_DEC, FIELD_NAMES):
    disk = Disk(dc2_ra, dc2_dec, DISK_RADIUS_DG * 3600)
    obj_list = cat.get_object_type_by_region(disk, 'galaxy')
    print('Field ', name)
    print(f'# objects: {len(obj_list)}')
    print(f'# collections: {obj_list.collection_count}')

    rotator = FieldRotator(dc2_ra, dc2_dec, f_ra, f_dec)
    galaxy_schema = make_galaxy_schema('rotate_logger')
    galaxy_rotator = GalaxyRotator(rotator, name, obj_list)
    galaxy_rotator.output_field_pixels(rotated_dir, galaxy_schema)

    break    ### for testing stop after 1 field
