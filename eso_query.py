import argparse
from pyvo.dal import tap
from astropy import units as u
from astroquery.mast import Catalogs
import numpy as np
from astropy.time import Time

ESO_TAP_OBS = 'http://archive.eso.org/tap_obs'
tap_obs = tap.TAPService(ESO_TAP_OBS)


def arg_parse():
    """Parse command line arguments."""
    desc = '''
    Code to search the ESO archive with RA/DEC (or a TIC ID if available).
    '''
    p = argparse.ArgumentParser('ESO QUERY', description=desc)
    p.add_argument('-ra', help='Target RA in degrees.', type=float)
    p.add_argument('-dec', help='Target DEC in degrees.', type=float)
    p.add_argument('-radius', help='Box search radius in arcminutes.',
                   default=2, type=float)
    p.add_argument('-out', help='Output directory.', default='./query')
    p.add_argument('-t', '--tic-id', required=False, help='TIC ID.', type=int)
    return p.parse_args()


def summarise_obs(t):
    if len(t) == 0:
        print('No results to show.')
        return
    # First, lets check to see if more than one object exists
    if np.unique(t['object']).shape[0] > 1:
        text = f'''
        Summary includes multiple targets: {",".join(np.unique(t['object']))}
        '''
    else:
        text = f'Summary for {t["object"][0]}'
    for instrument in np.unique(t['instrument']):
        text +='\n'
        t_ = t[t['instrument']==instrument]
        start = str(Time(np.min(t_['mjd_obs']), format='mjd').datetime)[:10]
        end = str(Time(np.max(t_['mjd_obs']), format='mjd').datetime)[:10]
        text += f'{instrument} {start} -> {end} [{len(t_):,} points]'
        for prog_id in np.unique(t_['prog_id']):
            t__ = t_[t_['prog_id']==prog_id]
            start = str(
                Time(np.min(t__['mjd_obs']), format='mjd').datetime
                )[:10]
            end = str(
                Time(np.max(t__['mjd_obs']), format='mjd').datetime
                )[:10]
            tmp = f'{prog_id} ({t__["pi_coi"][0].split("/")[0]})'
            text += f'\n\t{tmp:<40} {start} -> {end} [{len(t__):,} points]'
    return text


def do_query(ra, dec, radius):
    query = f"""
    SELECT *
    FROM
    (
    SELECT
        target, object, ra, dec, pi_coi, prog_id, instrument, telescope,
        exp_start, exposure, mjd_obs, dp_cat, datalink_url
        FROM dbo.raw
        WHERE dp_cat='SCIENCE'
        AND (instrument='ESPRESSO' OR instrument='HARPS' OR instrument='FEROS')
        AND dec BETWEEN -90 AND 90
    ) AS sub
    WHERE 1=CONTAINS(
                point('', sub.ra, sub.dec),
                circle('', {ra}, {dec}, {radius}))
    """
    print(query)
    res = tap_obs.search(query=query)
    return res.to_table()


if __name__ == '__main__':
    verbose = False
    args = arg_parse()
    # import pdb; pdb.set_trace()
    if args.tic_id is None:
        ra = float(args.ra) * u.deg
        dec = float(args.dec) * u.deg
    else:
        cat = Catalogs.query_object(f'TIC {args.tic_id}',
                                    radius=1e-3, catalog="TIC")
        ra = cat['ra'][0] * u.deg
        dec = cat['dec'][0] * u.deg
    radius = (float(args.radius) * u.arcmin).to(u.deg)
    out = args.out

    res = do_query(ra.value, dec.value, radius.value)
    res.write(f'{out}.dat', format='ascii.ecsv', overwrite=True)
    text = summarise_obs(res)
    print(text)
