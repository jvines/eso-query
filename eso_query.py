#!/Users/jayvains/miniforge3/bin/python
import argparse
from pyvo.dal import tap
from astropy import units as u
from astroquery.mast import Catalogs
import numpy as np
from astropy.time import Time
from difflib import SequenceMatcher
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import re

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
                   default=3, type=float)
    p.add_argument('-out', help='Output directory.', default='./query')
    p.add_argument('-t', '--tic-id', required=False, help='TIC ID.', type=int,
                   nargs='+')
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


def text_similarity(s1, s2):
    return SequenceMatcher(None, s1, s2).ratio()


def get_tic_id(target):
    catalog = Catalogs.query_object(target,
                                    radius=1e-3, catalog='TIC')
    return f'TIC-{catalog["ID"][0]}'


def summarize_multiple_observations(table):
    table['object'] = list(map(lambda x: fix_tic(x), table['object']))

    text = ''
    for target in np.unique(table['object']):
        mask = table['object'] == target
        n_points = len(table[mask])
        instruments = ';'.join(np.unique(table[mask]['instrument']))
        start_date = (Time(np.min(table[mask]['mjd_obs']), format='mjd')
                      .datetime
                      .strftime('%Y-%m-%d'))
        end_date = (Time(np.max(table[mask]['mjd_obs']), format='mjd')
                      .datetime
                      .strftime('%Y-%m-%d'))
        pis = []
        pi_cois = np.sort(table[mask]['pi_coi'])
        for pi in pi_cois:
            pi = pi.split('/')[0]
            if pi not in pis:
                if len(list(filter(lambda x: text_similarity(pi, x) > .75, pis))):
                    continue
                pis.append(pi)
        pis = f"{';'.join(pis)}"
        target_name = table[mask]['target'][0]
        text += f'{target_name}|({instruments})|{start_date} -> {end_date} [{n_points} points]|({pis})\n'
    return text


def do_query(ra, dec, radius):
    query = f"""
    SELECT *
    FROM
    (
        SELECT
            target
            , object
            , ra
            , dec
            , pi_coi
            , prog_id
            , instrument
            , telescope
            , exp_start
            , exposure
            , mjd_obs
        --    , dp_cat
            , datalink_url
        FROM dbo.raw
        WHERE dp_cat='SCIENCE'
            AND (instrument='ESPRESSO' OR instrument='HARPS' OR instrument='FEROS')
            AND dec BETWEEN -90 AND 90
    ) AS sub
    WHERE 1=CONTAINS(
                point('', sub.ra, sub.dec),
                circle('', {ra}, {dec}, {radius}))
    """
    res = tap_obs.search(query=query)
    return res.to_table()


def build_circle_condition(ras, decs, radius):
    template = "1=CONTAINS(point('', ra, dec), circle('', {}, {}, {}))"
    condition = template.format(ras[0].value, decs[0].value, radius)
    for ra, dec in zip(ras[1:], decs[1:]):
        condition += ' OR ' + template.format(ra.value, dec.value, radius)
    return condition


def do_multiple_query(ra, dec, radius):
    condition = build_circle_condition(ra, dec, radius)
    query = f"""
    SELECT
        target
        , object
        , ra
        , dec
        , pi_coi
        , prog_id
        , instrument
        , telescope
        , exp_start
        , exposure
        , mjd_obs
    --    , dp_cat
        , datalink_url
    FROM dbo.raw
    WHERE dp_cat='SCIENCE'
        AND (instrument='ESPRESSO' OR instrument='HARPS' OR instrument='FEROS')
        AND dec BETWEEN -90 AND 90
        AND ({condition})
    """
    res = tap_obs.search(query=query)
    return res.to_table()


def fix_tic(tic_id):
    if 'TIC' in tic_id:
        pattern = r'(TIC)\s?(\d+)'
        replacement = r"\1-\2"
        return re.sub(pattern, replacement, tic_id)
    elif 'TOI' in tic_id:
        return get_tic_id_from_toi(tic_id.split('-')[-1])
    else:
        try:
            return get_tic_id(tic_id)
        except:
            return tic_id


def get_tic_id_from_toi(toi):
    catalog = NasaExoplanetArchive.query_criteria(table='toi',
                                                  where=f'toipfx = {toi}')
    return f'TIC-{catalog["tid"][0]}'


def get_ra_dec_from_tic_id(tic_id):
    catalog = Catalogs.query_object(f'TIC {tic_id}',
                                    radius=1e-3, catalog='TIC')
    return catalog['ra'][0] * u.deg, catalog['dec'][0] * u.deg


def get_multiple_tic_ids_ra_dec(tic_ids):
    catalog = Catalogs.query_criteria(catalog='TIC', ID=tic_ids)
    return catalog['ra'] * u.deg, catalog['dec'] * u.deg


if __name__ == '__main__':
    args = arg_parse()
    if args.tic_id is None:
        ra = float(args.ra) * u.deg
        dec = float(args.dec) * u.deg
    elif type(args.tic_id) == list:
        ras, decs = get_multiple_tic_ids_ra_dec(args.tic_id)
    else:
        ra, dec = get_ra_dec_from_tic_id(args.tic_id)
    radius = (float(args.radius) * u.arcmin).to(u.deg)
    out = args.out
    if len(ras) and len(decs):
        results = do_multiple_query(ras, decs, radius.value)
        print(summarize_multiple_observations(results))
    else:
        resuls = do_query(ra.value, dec.value, radius.value)
        print(summarise_obs(results))
    # res.write(f'{out}.dat', format='ascii.csv', overwrite=True)
    # text = summarise_obs(res)
    # print(text)
