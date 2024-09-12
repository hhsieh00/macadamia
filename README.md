# MACADAMIA Archival Photometry

- MACADAMIA: Multi-Archive Catalog of Asteroid Detections And Measurements for Interactive Access
- written by Henry H. Hsieh (`hhsieh@psi.edu` / `hhsieh@gmail.com`), 2024-09-12


--------
Overview
--------

This python code uses the Canadian Astronomy Data Centre's Solar System Object Image Search (SSOIS) tool to search
for detections of known asteroids (and comets) in archival data, downloads the relevant data, and returns
photometry of those detections.

At the present time, only SDSS, CFHT/MegaCam, and Blanco/DECam archival data are supported, but adding additional
archives is relatively trivial (just need to add code blocks to specify basic data parameters and specific FITS
header keywords for certain metadata parameters).  The main caveat is that this code does not perform any image
processing (e.g., bias subtraction and flatfield correction) on downloaded data (SDSS, MegaCam and DECam archival
data are available already in processed form), and so photometry of archival data from other facilities may not
be reliable.


--------
Workflow
--------

The code carries out the following steps:

1. Places a query for a specific target from a specific facility to the SSOIS search tool.

2. Parses the query results to identify available data (along with metadata), rejecting data deemed to be unsuitable
   (at the moment, with exposure times of <30 seconds, or using u-band or Y-band filters which are unlikely to
   contain usable detections of small solar system objects).

3. Downloads data using links from the SSOIS query results.

4. Retrieves the target object's predicted position in each image using a query to the Horizons ephemeris tool
   using the midpoint of each image's exposure length as the ephemeris time.

5. Extracts data from individual FITS file extensions, computes astrometric solutions for each extension, and checks
   to see if the object's predicted position lies within that extension's field of view (FOV).  If an object's
   predicted position is found to be within a particular extension's FOV, that extracted extension data is saved
   and the code proceeds to the next step.

6. Identifies field sources in the data using DAOStarFinder, matches these against the Refcat2 photometric catalog,
   and computes a zeropoint for the image extension in question, provided that the image used a standard broadband
   filter that can be calibrated with the Refcat2 catalog (i.e., g', r', i', z', B, V, R, or I).  Also computes the
   point-source and surface brightness detection limits of each image using the computed zeropoint and background
   sky statistics.

7. Identifies sources close (<5 arcsec) to the expected position of the target object in each image and performs
   multi-aperture photometry on identified sources.


------
Output
------

The code produces a number of output files that are automatically sorted into various folders.  The organizational
structure of these folders and specific output files are as follows:

1. A folder containing all searches for a given object named "`archival_search_OBJECTNAME`" in the specified base
   directory.  Searches for different instruments and/or at different times will all be placed in this folder.

2. Individual folders for specific telescopes, search start date, and search end date named
   "`TELESCOPE_YYYYMMDD_YYYYMMDD`" in each object folder.

3. A file named `log_archival_photometry_YYYYMMDD_HHMMSS.txt` logging various processing data for each run of this
   code.

4. A file named `ssois_results_TARGETNAME_YYYYMMDD_YYYYMMDD.txt`: a table of SSOIS search results

5. A file named `ssois_results_TARGETNAME_YYYYMMDD_YYYYMMDD_ephems.txt`: a table of SSOIS search results with
   additional observational geometry data, including expected magnitudes according to Horizons, to aid in manual
   identification of usable data.

6. Original downloaded FITS files (compressed using `fpack`), located in a sub-folder named "`data_files_orig`".

7. Extracted extension data containing the predicted positions of the target object (as FITS files compressed using
   `fpack`), located in a sub-folder named "`data_files_ext_wcs`", where the extension number NN of each extension is
   recorded in each FITS file's filename (i.e., `*_NN_wcs.fits.fz`).

8. Background-subtracted extension data containing the predicted positions of the target object (as FITS files
   compressed using `fpack`) that were used for aperture photometry, located in a sub-folder named
   "`data_files_ext_wcs_bgsub`"

9. Various files pertaining to field source photometry, located in a sub-folder named "`field_source_photometry`".
   Included files include files following the formats of (where prefixes for all of these files indicate the
   corresponding image file):
   
   `*.field_source_photometry.txt`: x,y and RA,Dec coordinates, fluxes, signal-to-noise ratios, and computed
   				  magnitudes of field sources)

   `*.zeropoint_histogram?.pdf`: histogram plots of zeropoints corresponding to each field source used for absolute
                               photometric calibration, where histogram1 plots show the initial distribution of
			       calculated zeropoints and histogram2 plots show the distribution of calculated
			       zeropoints after outlier rejection (which is the distribution used to compute the
			       final zeropoint of each image)

   `*.field_source_apertures.pdf`: plots of each image with circles around detected field sources

   `*.field_sources.txt`: x,y and RA,Dec coordinates of detected field sources

   `*.refcat_calib_stars.txt`: table of field sources matched to Refcat2 catalog reference stars with catalogued
                             magnitudes in Pan-STARRS1 filters, and computed equivalent magnitudes in SDSS and
			     Johnson-Cousins broadband filters)

10. Various files pertaining to target photometry, located in a sub-folder named "`target_photometry`". Included files
    include files following the formats of:

    `*.target_candidate_sources.txt`: x,y and RA/Dec coordinates of detected target candidate sources

    `*.target_ap*arcsec.pdf`: plots of each image with circles around target detections showing the photometry
                            aperture and annuli used to measure the sky background

    `target_phot_results.ap*arcsec.txt`: tables of detection metadata and output photometry for all usable images
                                       found by SSOIS and have had photometry successfully measured


------------
Requirements
------------

1. A user-customized configuration file (named `macadamia.cfg` in the example file provided here) specifying file
   paths to certain required local installations of external software tools.  If this configuration file is
   located in the same location as the data processing, just the filename needs to be specified, but otherwise,
   the entire filepath should be specified.

2. A local installation of fpack and funpack (https://heasarc.gsfc.nasa.gov/fitsio/fpack/).  Paths to these
   commands should be specified in the user configuration file.

3. A local installation of the refcat search tool and catalog data (https://archive.stsci.edu/hlsp/atlas-refcat2).
   Paths to the search program and catalogs should be specified in the user configuration file.


------
Usage
------

The code can be run from the command line as follows:

`python3 macadamia_archival_photometry.py BASE_PATH macadamia.cfg TARGET_NAME TELESCOPE`

where the currently supported values for TELESCOPE are `SDSS`, `CFHT/MegaCam`, and `CTIO-4m/DECam`
(case-sensitive).  Additional optional parameters can be specified as follows:

    -start_date             (start date for SSOIS search)
    -end_date               (end date for SSOIS search)
    -fr [float]             (field source photometry radius in arcsec)
    -ffwhm [float]          (FWHM for DAOStarFinder field source detection)
    -fsig [float]           (sigma threshold for DAOStarFinder field source detection)
    -fskyin [float]         (inner sky annulus radius for field source photometry in arcsec)
    -fskyout [float]        (outer sky annulus radius for field source photometry in arcsec
    -tr [float] [float] ... (list of target photometry radii in arcsec)
    -tfwhm [float]          (FWHM for DAOStarFinder target candidate detection)
    -tsig [float]           (sigma threshold for DAOStarFinder target candidate detection)
    -tskyin [float]         (inner sky annulus radius for target photometry in arcsec)
    -tskyout [float]        (outer sky annulus radius for field source photometry in arcsec)

All downloaded data and photometry output will be placed in sub-folders under `BASE_PATH` so users should make
sure that sufficient disk space for downloaded data is available at that path.  If no start and end dates are
specified for SSOIS searches, these dates will be set to 1990-01-01 and the current date, respectively.

Efforts have been made to make processing is relatively resilient to interruptions to avoid unnecessary
re-processing.  If a SSOIS query results file corresponding to the same search terms is detected, SSOIS will not
be queried again.  Similarly, the code also checks to see if each individual data file has already been downloaded
and/or decompressed before downloading it.  We note that in the event of an interruption, it is the user's
responsibility to ensure that incompletely downloaded files are re-downloaded in their entirety (e.g., deleted
before running the code so that the code does not automatically skip them).

We also note that if processing over more than one day is interrupted, the code may not recognize that data was
previously downloaded because a new search for the current day will be created upon re-running the code.  To
continue the earlier processing, the user can simply manually change the end date of the search in the folder name
and ssois_results* files to reflect the current date.

There is currently no mechanism for resuming processing of data files interrupted mid-way through a given
multi-extension file (i.e., re-processing will always begin with the first extension).  If entire files have
already been processed, however, or the user wishes to skip certain files due to the `ssois_results_*_ephems.txt`
file indicating that the file is unlikely to contain usable data (e.g., the object is excessively faint at the
time), these can be deleted from the `ssois_results*.txt` file, which will lead the code to omit them from
processing.  In these cases, we recommend that the user saves a copy of the original `ssois_results_*.txt` file in
case they wish to revisit those omitted files.

Photometry will only be conducted for images using SDSS or Johnson-Cousins filters.  While images using certain
other filters (u and Y, for now) will be excluded from processing entirely, images using other wideband or
narrowband filters (e.g., VR, gri, CaHK) will still be searched for target objects, with positions recorded in the
processing log if found to enable the user to conduct their own manual searches for and analyses of detections of
the object at the predicted positions if desired.


-------
Caveats
-------

1. This code relies on accurate time data being stored in image FITS headers.  If times are incorrect, objects may
   not be found, or may be mis-identified (i.e., associated with an incorrect field source).

2. This code also relies on accurate ephemeris predictions.  If an object's ephemeris position has a large
   uncertainty at the time of a particular archival observation or is incorrect (e.g., due to un-accounted-for
   non-gravitational effects or spurious astrometry data being used to compute its orbit), it may not be found, or
   may be mis-identified (i.e., associated with an incorrect field source).
