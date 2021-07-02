#!/usr/bin/env python
import h5py as h5
import numpy as np
import plot_utils as pu


def reconstruct_uniform(h5f, var, level, gridlist, cu, center, smin, smax):
    dset, nd, levelmet = collect_dataset(h5f, var, level, gridlist)
    if not levelmet:
        return [], []

    if cu:
        inb, ind = pu.find_indices(nd, center, smin, smax, True)
        print('Ordered plot center', center[0], center[1], center[2], ' gives following uniform grid indices:', ind[0], ind[1], ind[2])
    else:
        ind = int(nd[0] / 2), int(nd[1] / 2), int(nd[2] / 2)

    xy = dset[:, :, ind[2]].swapaxes(0, 1)
    xz = dset[:, ind[1], :].swapaxes(0, 1)
    yz = dset[ind[0], :, :].swapaxes(0, 1)

    d3min, d3max = np.min(dset), np.max(dset)
    d2max = max(np.max(xz), np.max(xy), np.max(yz))
    d2min = min(np.min(xz), np.min(xy), np.min(yz))
    block = [True, True, True], [yz, xz, xy], smin, smax, level

    return [[block, ], ], [[d2min], [d2max], [d3min], [d3max]]


def collect_dataset(h5f, dset_name, level, gridlist):
    print('Reading', dset_name)
    attrs = h5f['domains']['base'].attrs
    nd = [i * 2**level for i in attrs['n_d']]
    nxd, nyd, nzd = attrs['n_d'][0:3]
    dset = np.zeros((nd[0], nd[1], nd[2]))

    print('Reconstructing domain from cg parts')
    levelmet = False
    for ig in gridlist:
        h5g = h5f['data']['grid_' + str(ig).zfill(10)]
        if h5g.attrs['level'] == level:
            levelmet = True
            off = h5g.attrs['off']
            ngb = h5g.attrs['n_b']
            n_b = [int(ngb[0]), int(ngb[1]), int(ngb[2])]
            ce = n_b + off
            dset[off[0]:ce[0], off[1]:ce[1], off[2]:ce[2]] = h5g[dset_name][:, :, :].swapaxes(0, 2)

    return dset, nd, levelmet


def collect_gridlevels(h5f, var, refis, extr, maxglev, plotlevels, gridlist, cgcount, center, usc, getmap):
    l2, h2, l3, h3 = extr
    for iref in range(maxglev + 1):
        if iref in plotlevels:
            print('REFINEMENT ', iref)
            blks = []
            for ib in gridlist:
                levok, block, extr = read_block(h5f, var, ib, iref, center, usc, getmap)
                if levok:
                    blks.append(block)
                    if getmap:
                        l2.append(extr[0])
                        h2.append(extr[1])
                        l3.append(extr[2])
                        h3.append(extr[3])
            if blks != []:
                refis.append(blks)
    return refis, [l2, h2, l3, h3]


def read_block(h5f, dset_name, ig, olev, oc, usc, getmap):
    h5g = h5f['data']['grid_' + str(ig).zfill(10)]
    level = h5g.attrs['level']
    levok = (level == olev)
    if not levok:
        return levok, [], []

    ledge = h5g.attrs['left_edge']
    redge = h5g.attrs['right_edge']
    ngb = h5g.attrs['n_b']
    inb, ind = pu.find_indices(ngb, oc, ledge, redge, False)
    if not any(inb):
        return False, [], []
    if not getmap:
        return levok, [inb, [], ledge / usc, redge / usc, olev], []
    clen = h5g.attrs['dl']
    off = h5g.attrs['off']
    n_b = [int(ngb[0]), int(ngb[1]), int(ngb[2])]
    ce = n_b + off
    dset = h5g[dset_name][:, :, :].swapaxes(0, 2)
    xy = dset[:, :, ind[2]].swapaxes(0, 1)
    xz = dset[:, ind[1], :].swapaxes(0, 1)
    yz = dset[ind[0], :, :].swapaxes(0, 1)
    d3min, d3max = np.min(dset), np.max(dset)
    d2max = max(np.max(xz), np.max(xy), np.max(yz))
    d2min = min(np.min(xz), np.min(xy), np.min(yz))

    return levok, [inb, [yz, xz, xy], ledge / usc, redge / usc, olev], [d2min, d2max, d3min, d3max]


def collect_particles(h5f, drawh, center, player, uupd, usc, plotlevels, gridlist):
    if 'particle_types' not in list(h5f):
        return False, [], []
    print('Reading particles')
    px, py, pz, pm = np.array([]), np.array([]), np.array([]), np.array([])
    for gn in h5f['data']:
        if h5f['data'][gn].attrs['level'] in plotlevels and int(gn.split('grid_')[-1]) in gridlist:
            if str(player[1]) == '0' and str(player[2]) == '0' and str(player[3]) == '0':
                px = np.concatenate((px, h5f['data'][gn]['particles']['stars']['position_x'][:]))
                py = np.concatenate((py, h5f['data'][gn]['particles']['stars']['position_y'][:]))
                pz = np.concatenate((pz, h5f['data'][gn]['particles']['stars']['position_z'][:]))
                if drawh:
                    pm = np.concatenate((pm, h5f['data'][gn]['particles']['stars']['mass'][:]))
            else:
                apx = h5f['data'][gn]['particles']['stars']['position_x'][:]
                apy = h5f['data'][gn]['particles']['stars']['position_y'][:]
                apz = h5f['data'][gn]['particles']['stars']['position_z'][:]
                maskx = np.abs(apx - center[0]) <= float(player[1])
                masky = np.abs(apy - center[1]) <= float(player[2])
                maskz = np.abs(apz - center[2]) <= float(player[3])
                if player[0]:
                    auxm = str(player[1]) == '0' or str(player[2]) == '0' or str(player[3]) == '0'
                    for i in range(np.size(maskx)):
                        maskx[i] = maskx[i] or masky[i] or maskz[i] or auxm
                else:
                    for i in range(np.size(maskx)):
                        maskx[i] = maskx[i] and masky[i] and maskz[i]
                px = np.concatenate((px, apx[maskx]))
                py = np.concatenate((py, apy[maskx]))
                pz = np.concatenate((pz, apz[maskx]))
                if drawh:
                    pm = np.concatenate((pm, h5f['data'][gn]['particles']['stars']['mass'][maskx]))

    if px.size == 0:
        return False, [], []
    if uupd:
        return True, pu.list3_division([px, py, pz], usc), pm

    return True, [px, py, pz], pm
