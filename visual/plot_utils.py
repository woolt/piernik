#!/usr/bin/env python
import numpy as np


def fsym(vmin, vmax):
    vmx = np.max([np.abs(vmin), np.abs(vmax)])
    vmn = -1.0 * vmx
    if vmn == vmx:
        vmn = vmn - 0.00001
        vmx = vmx + 0.00001
    return vmn, vmx


def scale_manage(sctype, refis, umin, umax, d2min, d2max):

    symmin = 1.0
    if (umin == 0.0 and umax == 0.0):
        vmin, vmax = d2min, d2max
    else:
        vmin, vmax = umin, umax

    if (sctype == '1' or sctype == 'symlin'):
        vmin, vmax = fsym(vmin, vmax)

    elif (sctype == '2' or sctype == 'log'):
        if (vmin > 0.0):
            vmin = np.log10(vmin)
        else:
            vmin = np.log10(check_minimum_data(refis, False)[0])
        if (vmax > 0.0):
            vmax = np.log10(vmax)
        else:
            vmax = -1.
    elif (sctype == '3' or sctype == 'symlog'):
        if (umin > 0.0 and umax > 0.0):
            symmin = umin
            vmax = np.log10(umax / umin)
            vmin = np.log10(vmax)
        else:
            if (d2min * d2max > 0.0):
                smin, smax = d2min, d2max
            else:
                smin, smax = check_minimum_data(refis, True)
            symmin = min(np.abs(smin), np.abs(smax))
            vmax = np.log10(max(np.abs(smin), np.abs(smax)) / symmin)
        vmin = -vmax
        print('SYMMIN value for SYMLOG scaletype: %s' % symmin)

    return vmin, vmax, symmin


def check_minimum_data(refis, extended):
    cmdmin, cmdmax = np.inf, -np.inf
    for blks in refis:
        for bl in blks:
            binb, bxyz = bl[0:2]
            for ncut in range(3):
                if binb[ncut]:
                    cmdmin = min(cmdmin, np.min(bxyz[ncut], initial=np.inf, where=(bxyz[ncut] > 0.0)))
                    if extended:
                        cmdmax = max(cmdmax, np.max(bxyz[ncut], initial=-np.inf, where=(bxyz[ncut] < 0.0)))
    return cmdmin, cmdmax


def scale_plotarray(pa, sctype, symmin):
    if (sctype == '2' or sctype == 'log'):
        pa = np.log10(pa)
    elif (sctype == '3' or sctype == 'symlog'):
        pa = np.sign(pa) * np.log10(np.maximum(np.abs(pa) / symmin, 1.0))
    return pa


def list3_division(l3, divisor):
    return l3[0] / divisor, l3[1] / divisor, l3[2] / divisor


def labelx():
    return lambda var: '$' + str(var)[2:-1].replace('**', '^') + '$'


def take_nonempty(lst):
    for it in lst:
        if it != []:
            return it
    return []


def colorbar_mode(drawd, drawh):
    if drawd and drawh:
        cbar_mode = 'none'
    elif drawd or drawh:
        cbar_mode = 'single'
    else:
        cbar_mode = 'none'
    return cbar_mode


def color_axes(wax, color):
    wax.spines['top'].set_color(color)
    wax.spines['bottom'].set_color(color)
    wax.spines['left'].set_color(color)
    wax.spines['right'].set_color(color)
    wax.tick_params(axis='x', colors=color)
    wax.tick_params(axis='y', colors=color)
    return


def detindex(nd, cxyz, smin, smax):
    return int(np.floor(nd * (cxyz - smin) / (smax - smin)))


def ind_limits(nd, cxyz, smin, smax):
    return int(min(nd - 1, max(0, detindex(nd, cxyz, smin, smax))))


def isinbox(cxyz, smin, smax, warn, cc):
    isin = (cxyz >= smin and cxyz <= smax)
    if not isin and warn:
        print('Domain edges %s %s used to plot as the given plot center %s coordinate (%s) is outside the domain.' % (smin, smax, cc, cxyz))
    return isin


def find_indices(nd, cxyz, smin, smax, warn):
    inb = isinbox(cxyz[0], smin[0], smax[0], warn, 'CX'), isinbox(cxyz[1], smin[1], smax[1], warn, 'CY'), isinbox(cxyz[2], smin[2], smax[2], warn, 'CZ')
    icc = ind_limits(nd[0], cxyz[0], smin[0], smax[0]), ind_limits(nd[1], cxyz[1], smin[1], smax[1]), ind_limits(nd[2], cxyz[2], smin[2], smax[2])
    return inb, icc


def convert_units(infile, toplot):
    au_cm = 1.49597870700e13
    pc_au = 206264.806248712
    pc_cm = pc_au * au_cm
    if infile == 'pc':
        cm = 1.0 / pc_cm
        pc = 1.0
    elif infile == 'au':
        cm = 1.0 / au_cm
        pc = pc_cm * cm
    elif infile == 'kpc':
        cm = 1.0 / (1.0e3 * pc_cm)
        pc = 0.001
    elif infile == 'm':
        cm = 1.0 / 1.0e2
        pc = pc_cm * cm
    else:
        return 1., False

    if toplot == 'cm':
        return cm, True
    elif toplot == 'metr':
        return 1.0e2 * cm, True
    elif toplot == 'km':
        return 1.0e5 * cm, True
    elif toplot == 'au':
        return au_cm * cm, True
    elif toplot == 'pc':
        return pc, True
    elif toplot == 'kpc':
        return 1.0e3 * pc, True
    elif toplot == 'Mpc':
        return 1.0e6 * pc, True
    elif toplot == 'lyr':
        return 9.4605e17 * cm, True
    else:
        return 1., False


def change_units(fromfile, toplot):
    infile = fromfile.decode('utf-8')
    if infile == toplot or toplot == '':
        return 1., fromfile, False
    if toplot == 'k':
        return 1.e3, b"".join([b'k', fromfile]), True
    if toplot == 'M':
        return 1.e6, b"".join([b'M', fromfile]), True
    if toplot == 'G':
        return 1.e9, b"".join([b'G', fromfile]), True
    if toplot == 'm':
        return 1.e-3, b"".join([b'm', fromfile]), True
    if toplot == 'mu':
        return 1.e-6, b"".join([b'mu', fromfile]), True
    conv, chan = convert_units(infile, toplot)
    if chan:
        return conv, bytes(toplot, 'utf-8'), chan
    else:
        return conv, fromfile, chan
