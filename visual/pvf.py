#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt
import h5py
import os
import sys
import plot_compose as pc

cmap = 'viridis'
pcolor = 'default'
plotdir = 'frames'
sctype = 'linear'
cu, cx, cy, cz = False, 0.0, 0.0, 0.0
zmin, zmax = 0.0, 0.0
draw_part = False
draw_data = False
dnames = ''
uaxes = ''
nbins = 1
psize = 0
exten = '.png'

print('PIERNIK VISUALIZATION FACILITY')


def print_usage():
    print('Usage: ./pvf.py <file> [options]')
    print('')
    print('Options:')
    print(' -h, \t\t--help \t\t\tprint this help')
    print(' -a UNIT, \t--axes UNIT \t\tscale plot axes with UNIT [default: dataset units]')
    print(' -b BINS, \t--bins BINS \t\tmake a 2D histogram plot using BINS number instead of scattering particles [default: 1, which leads to scattering]')
    print(' -c CX,CY,CZ, \t--center CX,CY,CZ \tplot cuts across given point coordinates CX, CY, CZ [default: computed domain center]')
    print(' -d VAR[,VAR2], --dataset VAR[,VAR2] \tspecify one or more datafield(s) to plot [default: print available datafields; all or _all_ to plot all available datafields]')
    print(' -D COLORMAP, \t--colormap COLORMAP \tuse COLORMAP palette [default: viridis]')
    print(' -e EXTENSION, \t--extension EXTENSION \tsave plot in file using filename extension EXTENSION [default: png]')
    print(' -l SCALETYPE, \t--scale SCALETYPE \tdump use SCALETYPE scale type for displaying data (possible values: 0 | linear, 1 | symlin, 2 | log, 3 | symlog) [default: linear]')
    print(' -o OUTPUT, \t--output OUTPUT \tdump plot files into OUTPUT directory [default: frames]')
    print(' -p,\t\t--particles\t\tscatter particles onto slices [default: switched-off]')
    print(' -P,\t\t--particle-color\tuse color for particles scattering or colormap for particles histogram plot [default: #1f77b4 (blue) or viridis]')
    print(' -s,\t\t--particle-sizes\tmarker sizes for scattering particles onto slices [default: switched-off]')
    print(' -z ZMIN,ZMAX, \t--zlim ZMIN,ZMAX \tlimit colorscale to ZMIN and ZMAX [default: computed data maxima symmetrized]')


def cli_params(argv):
    try:
        opts, args = getopt.getopt(argv, "a:b:c:d:D:e:hl:o:pP:s:z:", ["help", "axes=", "bins=", "center=", "colormap=", "dataset=", "extension=", "output=", "particles", "particle-color=", "particle-sizes=", "scale=", "zlim="])
    except getopt.GetoptError:
        print("Unrecognized options: %s \n" % argv)
        print_usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print_usage()
            sys.exit()

        elif opt in ("-a", "--axes"):
            global uaxes
            uaxes = str(arg)

        elif opt in ("-b", "--bins"):
            global nbins
            nbins = int(arg)

        elif opt in ("-c", "--center"):
            global cx, cy, cz, cu
            cx, cy, cz = arg.split(',')
            cu, cx, cy, cz = True, float(cx), float(cy), float(cz)

        elif opt in ("-d", "--dataset"):
            global dnames
            global draw_data
            dnames = str(arg)
            draw_data = True

        elif opt in ("-D", "--colormap"):
            global cmap
            cmap = str(arg)

        elif opt in ("-e", "--extension"):
            global exten
            exten = '.' + str(arg)
            print(exten)

        elif opt in ("-l", "--scale"):
            global sctype
            sctype = str(arg)

        elif opt in ("-o", "--output"):
            global plotdir
            plotdir = str(arg)
            print('PLOTDIR: ', plotdir)

        elif opt in ("-p", "--particles"):
            global draw_part
            draw_part = True

        elif opt in ("-P", "--particle-color"):
            global pcolor
            pcolor = str(arg)

        elif opt in ("-s", "--particle-sizes"):
            global psize
            psize = float(arg)

        elif opt in ("-z", "--zlim"):
            global zmin, zmax
            zmin, zmax = arg.split(',')
            zmin = float(zmin)
            zmax = float(zmax)
            print("zmin, zmax = ", zmin, zmax)


if (len(sys.argv) < 2):
    print_usage()
    exit()

iw = 1
for word in sys.argv[1:]:
    if word[0] == '-':
        break
    iw += 1

cli_params(sys.argv[iw:])
if pcolor == 'default':
    if nbins > 1:
        pcolor = 'viridis'
    else:
        pcolor = '#1f77b4'
options = zmin, zmax, cmap, pcolor, psize, sctype, cu, cx, cy, cz, draw_data, draw_part, nbins, uaxes
if not os.path.exists(plotdir):
    os.makedirs(plotdir)


files_list = sys.argv[1:iw]
for pthfilen in files_list:
    print('')
    file_exists = os.path.exists(pthfilen)
    if not file_exists:
        print('The file %s does not exist!' % pthfilen)
        continue
    if pthfilen.split('.')[-1] == 'h5':
        h5f = h5py.File(pthfilen, 'r')
    else:
        continue
    if not (draw_data or draw_part) or (draw_data and dnames == ''):
        partincl = ''
        if 'particle_types' in list(h5f):
            partincl = 'and particles'
        print('Available datafields in the file %s: \n' % pthfilen, list(h5f['field_types'].keys()), partincl)
        h5f.close()
        continue
    filen = pthfilen.split('/')[-1]

    print("Reading file: %s" % pthfilen)
    prd, prp = '', '',
    if draw_data:
        if dnames == "_all_" or dnames == "all":
            varlist = h5f['field_types'].keys()
        else:
            varlist = dnames.split(',')
        prd = 'datasets: %s' % varlist
        if draw_part:
            prp = 'particles and '
    else:
        varlist = ['part']
        prp = 'particles only'
    print('Going to read ' + prp + prd)

    for var in varlist:
        if (not draw_data or var in list(h5f['field_types'].keys())):
            # output = plotdir+'/'+filen.split('/')[-1].replace('.h5',"_%s.png" % var)
            fnl = filen.split('/')[-1]
            output = plotdir + '/' + '_'.join(fnl.split('_')[:-1]) + '_' + var + '_' + fnl.split('_')[-1].replace('.h5', exten)
            pc.plotcompose(pthfilen, var, output, options)
        else:
            print(var, ' is not available in the file ', pthfilen)

    h5f.close()
