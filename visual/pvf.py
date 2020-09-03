#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt, h5py, os, sys
import plot_compose as pc

if (len(sys.argv) < 3):
    print('PIERNIK VISUALIZATION FACILITY')
    print('Usage: ./pvf.py <file> <varname,[varname,...]> [options]')
    if len(sys.argv) < 2:
        exit()

cmap    = 'viridis'
plotdir = 'frames'
zmin = 0.0
zmax = 0.0

def cli_params(argv):
    try:
        opts,args=getopt.getopt(argv,"ho:r:z:",["help","colormap=","output=","zlim="])
    except getopt.GetoptError:
        print("Unidentified error.")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(" -h, \t\t--help \t\t\tprint this help \n  \
               -o OUTPUT, \t--output OUTPUT \tdump plot files into OUTPUT directory \n \
               -r COLORMAP, \t--colormap COLORMAP \tuse COLORMAP palette \n \
               -z ZMIN,ZMAX, \t--zlim ZMIN,ZMAX \tlimit colorscale to ZMIN and ZMAX")
            sys.exit()

        elif opt in ("-o", "--output"):
            global plotdir
            plotdir = str(arg)
            print('PLOTDIR: ', plotdir)

        elif opt in ("-r", "--colormap"):
            global cmap
            cmap = str(arg)

        elif opt in ("-z", "--zlim"):
            global zmin, zmax
            zmin, zmax = arg.split(',')
            zmin = float(zmin)
            zmax = float(zmax)
            print("zmin, zmax = ", zmin, zmax)

cli_params(sys.argv[3:])

pthfilen = sys.argv[1]
filen  = pthfilen.split('/')[-1]
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

h5f = h5py.File(pthfilen,'r')
if len(sys.argv) < 3:
    print("Available datafields: ", list(h5f['field_types'].keys()))
    exit(1)
if sys.argv[2] == "_all_":
    varlist = h5f['field_types'].keys()
else:
    varlist  = sys.argv[2].split(',')

print(varlist)

print("Reading file: %s" % pthfilen)
for var in varlist:
    #output = plotdir+'/'+filen.split('/')[-1].replace('.h5',"_%s.png" % var)
    fnl = filen.split('/')[-1]
    output = plotdir+'/'+'_'.join(fnl.split('_')[:-1])+'_'+var+'_'+fnl.split('_')[-1].replace('.h5',".png")
    options = zmin, zmax, cmap
    pc.plotcompose(pthfilen, var, output, options)

h5f.close()
