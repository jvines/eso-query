import argparse
from pyvo.dal import tap
from astropy import units as u

ESO_TAP_OBS = 'http://archive.eso.org/tap_obs'
tap_obs = tap.TAPService(ESO_TAP_OBS)


def arg_parse():
    """Parse command line arguments."""
    p = argparse.ArgumentParser()
    p.add_argument('ra', help='Target RA in degrees.')
    p.add_argument('dec', help='Target DEC in degrees.')
    p.add_argument('radius', help='Box search radius.')
    return p.parse_args()


if __name__ == '__main__':
    args = arg_parse()
    ra = float(args.ra) * u.deg
    dec = float(args.dec) * u.deg
    radius = (float(args.radius) * u.arcmin).to(u.deg)
    ra_min = (ra - radius).value
    ra_max = (ra + radius).value
    dec_min = (dec - radius).value
    dec_max = (dec + radius).value
    n = 100

    query = """
    SELECT TOP {} object, ra, dec, prog_id, instrument, telescope, exp_start, exposure, mjd_obs, dp_cat
    from dbo.raw
    where ra between {ra_min} and {ra_max}
      and dec between {dec_min} and {dec_max}
      and dp_cat='SCIENCE'
      and (instrument='ESPRESSO' or instrument='HARPS' or instrument='FEROS')
    """.format(n, ra_min, ra_max, dec_min, dec_max)
    print(query)
    res = tap_obs.search(query=query)
    print(res)