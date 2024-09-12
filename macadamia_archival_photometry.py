import sys,os
import argparse
from macadamia_functions import *


def parseArguments():
    parser = argparse.ArgumentParser(description = 'Command line arguments for MACADAMIA archival photometry script')
    parser.add_argument('base_path', help='Base path', type=str)
    parser.add_argument('config_file', help='Configuration file', type=str)
    parser.add_argument('target_name', help='Target name', type=str)
    parser.add_argument('telinst', help='Telescope/Instrument', type=str)
    parser.add_argument('-fr', '--field_source_phot_radius_arcsec', help = 'Field source photometry radius (arcsec)', required = False, default = 0)
    parser.add_argument('-ffwhm', '--daofind_fwhm', help = 'FWHM for DAOStarFinder field source detection', required = False, default = 0)
    parser.add_argument('-fsig', '--daofind_sigma_threshold', help = 'Sigma threshold for DAOStarFinder field source detection', required = False, default = 0)
    parser.add_argument('-fskyin', '--sky_inner_r_arcsec', help = 'Inner sky annulus radius for field source photometry (arcsec)', required = False, default = 0)
    parser.add_argument('-fskyout', '--sky_outer_r_arcsec', help = 'Outer sky annulus radius for field source photometry (arcsec)',required = False, default = 0)
    parser.add_argument('-tr', '--target_phot_radii_arcsec', nargs='*', type=float, help = 'Target photometry radii (arcsec)', required = False, default = [])
    parser.add_argument('-tfwhm', '--daofind_fwhm_target', help = 'FWHM for DAOStarFinder target candidate detection', required = False, default = 0)
    parser.add_argument('-tsig', '--daofind_sigma_threshold_target', help = 'Sigma threshold for DAOStarFinder target candidate detection', required = False, default = 0)
    parser.add_argument('-tskyin', '--target_sky_inner_r_arcsec', help = 'Inner sky annulus radius for target photometry (arcsec)', required = False, default = 0)
    parser.add_argument('-tskyout', '--target_sky_outer_r_arcsec', help = 'Outer sky annulus radius for field source photometry (arcsec)', required = False, default = 0)
    args = parser.parse_args()
    return args

def main():

    args = parseArguments()

    param_dict = {
        'base_path':args.base_path,
        'config_file':args.config_file,
        'target_name':args.target_name,
        'telinst':args.telinst
    }

    print(param_dict)
    if args.daofind_fwhm != 0:
        param_dict['daofind_fwhm'] = args.daofind_fwhm
    if args.daofind_sigma_threshold != 0:
        param_dict['daofind_sigma_threshold'] = args.daofind_sigma_threshold
    if args.sky_inner_r_arcsec != 0:
        param_dict['sky_inner_r_arcsec'] = args.sky_inner_r_arcsec
    if args.sky_outer_r_arcsec != 0:
        param_dict['sky_outer_r_arcsec'] = args.sky_outer_r_arcsec
    if args.field_source_phot_radius_arcsec != 0:
        param_dict['field_source_phot_radius_arcsec'] = args.field_source_phot_radius_arcsec

    if args.daofind_fwhm_target != 0:
        param_dict['daofind_fwhm_target'] = args.daofind_fwhm_target
    if args.daofind_sigma_threshold_target != 0:
        param_dict['daofind_sigma_threshold_target'] = args.daofind_sigma_threshold_target
    if args.target_sky_inner_r_arcsec != 0:
        param_dict['target_sky_inner_r_arcsec'] = args.target_sky_inner_r_arcsec
    if args.target_sky_outer_r_arcsec != 0:
        param_dict['target_sky_outer_r_arcsec'] = args.target_sky_outer_r_arcsec
    if args.target_phot_radii_arcsec != []:
        param_dict['target_phot_radii_arcsec'] = args.target_phot_radii_arcsec
        
    search_extract_archival_photometry(param_dict)

    return None


if __name__ == '__main__':
    main()

