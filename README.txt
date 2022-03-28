CALCIUM (STATIONARY CELLS)
DOI: 10.5281/zenodo.6390919

OVERVIEW:

This code analyses intracellular calcium release. It’s designed for fluorescence microscopy experiments where cells are labelled with a calcium indicator. With this version it’s important that cells are stationary throughout and in place from the start. If some cells are moving use calcium.m instead.

INSTRUCTIONS:

Change the input parameters at the top of the main file (calciumStationaryCells.m). Make sure the input directory contains a single TIFF file.

DESCRIPTION:

The code reads in a single TIFF file representing a 3D time series. After background subtraction, flat-field correction and Gaussian smoothing, local maxima in each frame are identified. These are recorded (for a subset of frames) in Peaks.tif.

Maxima from each frame are combined into tracks using a nearest-neighbour approach. If the standard deviation of the position is greater than half the cell radius, the track is discarded. Otherwise the mean position is assumed to be the location of a cell, and recorded in Cells.tif and Labelled_cells.tif.

The mean intensity of a cell-sized disk around each of these fixed positions is calculated, forming intensity traces for each cell. Noise is reduced by Gaussian smoothing. Sharp increases in intensity are identified from the first derivative and a user-defined threshold determines which of these are classified as calcium spikes. If a minimum spike intensity and/or minimum spike duration have been specified at the start, spikes not meeting these criteria are discarded.

A spike starts at a peak in the first derivative and ends when the intensity drops below its level at the start of the spike. A second spike can only begin once the first has ended. Graphs of intensity and intensity gradient with respect to time are saved in the Calcium_release and No_calcium_release folders.

The duration, maximum intensity and integrated intensity of each spike are recorded. The areas corresponding to integrated intensity are shaded on the intensity/time graphs in grey. Data relating to individual spikes are recorded in Spikes.xlsx and data relating to individual cells in Results.xlsx. The % of "triggered" cells (cells which have released calcium) is plotted with respect to time in Triggered.png.

AUTHOR:

Jane Humphrey (janehumphrey@outlook.com)

LICENSE:

MIT (see LICENSE.txt file for details)

DEPENDENCIES:

MATLAB v9.5
Signal Processing Toolbox v8.1
Image Processing Toolbox v10.3
Statistics and Machine Learning Toolbox v11.4

The code has been tested on Windows 8, Windows 10 and MacOS operating systems.

INSTALLATION:

MATLAB can be obtained from https://www.mathworks.com/products/matlab.html. No further installation is required.

ACKNOWLEDGEMENTS:

The functions pkfnd.m and track.m were written by Daniel Blair and Eric Dufresne. More information can be found at http://site.physics.georgetown.edu/matlab.

The export_fig toolbox was written by Yair Altman (https://github.com/altmany/export_fig).

Thanks to Dr Aleks Ponjavic for inspiration and Prof. Sir David Klenerman for support.
