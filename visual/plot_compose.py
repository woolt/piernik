#!/usr/bin/python
# -*- coding: utf-8 -*-

import h5py
import matplotlib
import numpy as np
import pylab as P
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import plot_utils as pu
import pvf_settings as ps
import read_dataset as rd
import time as timer
# matplotlib.use('cairo')      # choose output format


def plot_axes(ax, ulen, l1, min1, max1, l2, min2, max2):
    ax.set_xlim(min1, max1)
    ax.set_ylim(min2, max2)
    ax.set_xlabel("%s [%s]" % (l1, pu.labelx()(ulen)))
    ax.set_ylabel("%s [%s]" % (l2, pu.labelx()(ulen)))
    return ax


def plot1d(refis, field, parts, equip1d, ncut, n1, n2):
    smin, smax, zoom, ulen, umin, umax, linstyl, output, timep = equip1d
    fig1d = P.figure(ncut + 2, figsize=(10, 8))
    ax = fig1d.add_subplot(111)
    P.xlim(zoom[1][ncut], zoom[2][ncut])
    label = []
    axis = "xyz"[ncut]

    if field[0]:
        vmin, vmax, sctype, symmin, cmap, autsc, labf = field[1:]
        label.append(labf)
        if not autsc:
            P.ylim(vmin, vmax)

        for blks in refis:
            for bl in blks:
                binb, ble, bre, level, b1d = bl[1:]
                if binb[n1] and binb[n2]:
                    if b1d != []:
                        bplot = pu.scale_plotarray(b1d[ncut], sctype, symmin)
                        dxh = (bre[ncut] - ble[ncut]) / float(len(b1d[ncut])) / 2.0
                        vax = np.linspace(ble[ncut] + dxh, bre[ncut] - dxh, len(b1d[ncut]))
                        ax.plot(vax, bplot, linstyl[level], color=ps.plot1d_linecolor)

    if parts[0]:
        pxyz, pm, nbins, pcolor, psize, player, labh = parts[1:]
        pn1, pmm = pxyz[ncut], pm
        ax = plot1d_particles(ax, pn1, pmm, nbins, [smin[ncut], smax[ncut]], pcolor, psize)
        label.append(labh)

    if not field[0]:
        if umin != umax:
            P.ylim(umin, umax)

    P.ylabel('  |  '.join(label))
    P.xlabel("%s [%s]" % (axis, pu.labelx()(ulen)))
    P.title(timep)
    P.tight_layout()
    P.draw()
    out1d = output[0] + axis + '_' + output[1]
    fig1d.savefig(out1d, facecolor=ps.f_facecolor)
    print(out1d, "written to disk")
    P.clf()


def draw_plotcomponent(ax, refis, field, parts, equip2d, ncut, n1, n2):
    smin, smax, zoom, ulen, drawg, gcolor, center = equip2d
    ag, ah = [], []
    if field[0] or drawg:
        if field[0]:
            vmin, vmax, sctype, symmin, cmap = field[1:-2]
        for blks in refis:
            for bl in blks:
                bxyz, binb, ble, bre, level = bl[:-1]
                if binb[ncut]:
                    if bxyz != []:
                        bplot = pu.scale_plotarray(bxyz[ncut], sctype, symmin)
                        ag = ax.imshow(bplot, origin="lower", extent=[ble[n1], bre[n1], ble[n2], bre[n2]], vmin=vmin, vmax=vmax, interpolation='nearest', cmap=cmap)
                    if drawg:
                        ax.plot([ble[n1], ble[n1], bre[n1], bre[n1], ble[n1]], [ble[n2], bre[n2], bre[n2], ble[n2], ble[n2]], '-', linewidth=ps.grid_linewidth, alpha=0.1 * float(level + 1), color=gcolor[level], zorder=4)
    if parts[0]:
        pxyz, pm, nbins, pcolor, psize, player = parts[1:-1]
        if player[0] and player[ncut + 1] != '0':
            pmask = np.abs(pxyz[ncut] - center[ncut]) <= float(player[ncut + 1])
            pn1, pn2 = pxyz[n1][pmask], pxyz[n2][pmask]
            if nbins > 1:
                pmm = pm[pmask]
            else:
                pmm = pm
        else:
            pn1, pn2, pmm = pxyz[n1], pxyz[n2], pm
        ax, ah = draw_particles(ax, pn1, pn2, pmm, nbins, [smin[n1], smax[n1], smin[n2], smax[n2]], field[0], pcolor, psize)
    ax = plot_axes(ax, ulen, "xyz"[n1], zoom[1][n1], zoom[2][n1], "xyz"[n2], zoom[1][n2], zoom[2][n2])
    ax.set_xticks([center[n1]], minor=True)
    ax.set_yticks([center[n2]], minor=True)
    ax.tick_params(axis='x', which='minor', color=ps.centertick_color, bottom='on', top='on', width=ps.centertick_width, length=ps.centertick_length)
    ax.tick_params(axis='y', which='minor', color=ps.centertick_color, left='on', right='on', width=ps.centertick_width, length=ps.centertick_length)
    return ax, ag, ah


def draw_particles(ax, p1, p2, pm, nbins, ranges, drawd, pcolor, psize):
    if nbins > 1:
        ah = ax.hist2d(p1, p2, nbins, weights=pm, range=[ranges[0:2], ranges[2:4]], norm=matplotlib.colors.LogNorm(), cmap=pcolor, alpha=ps.hist2d_alpha)
        if not drawd:
            ax.set_facecolor(ps.hist2d_facecolor)
    else:
        ah = []
        if psize <= 0:
            psize = matplotlib.rcParams['lines.markersize']**2
        ax.scatter(p1, p2, c=pcolor, marker=ps.particles_marker, s=psize)
    return ax, ah


def plot1d_particles(ax, p1, pm, nbins, ranges, pcolor, psize):
    if nbins > 1:
        ax.hist(p1, nbins, weights=pm, range=ranges, color=ps.particles_color, alpha=ps.hist2d_alpha)
    else:
        if psize <= 0:
            psize = matplotlib.rcParams['lines.markersize']**2
        ax.scatter(p1, p1 * 0., c=pcolor, marker=ps.particles_marker, s=psize)
    return ax


def add_cbar(figmode, cbar_mode, grid, ab, ic, clab):
    icb = pu.cbsplace[figmode][ic]
    clf = [ps.cbar_hist2d_label_format, ps.cbar_plot2d_label_format][ic]
    if cbar_mode == 'none':
        axg = grid[icb]
        if figmode == 4 or figmode == 5:
            axg.set_xlim((grid[icb - 1].get_xlim()[0] / 2., grid[icb - 1].get_xlim()[1] / 2.))
        pu.color_axes(axg, 'white')
        fr = [ps.cbar_hist2d_shift, ps.cbar_plot2d_shift]
        bar = inset_axes(axg, width='100%', height='100%', bbox_to_anchor=(fr[ic], 0.0, 0.06, 1.0), bbox_transform=axg.transAxes, loc=2, borderpad=0)
        cbarh = P.colorbar(ab, cax=bar, format=clf, drawedges=False)
    else:
        if cbar_mode == 'each':
            bar = grid.cbar_axes[icb]
        elif cbar_mode == 'single':
            bar = grid.cbar_axes[0]
        bar.axis["right"].toggle(all=True)
        cbarh = P.colorbar(ab, cax=bar, format=clf, drawedges=False)
    cbarh.ax.set_ylabel(clab)
    if cbar_mode == 'none':
        cbarh.ax.yaxis.set_label_coords(ps.cbar_label_coords[0], ps.cbar_label_coords[1])


def plotcompose(pthfilen, var, output, options):
    axc, umin, umax, cmap, pcolor, player, psize, sctype, cu, center, drawg, drawd, drawu, drawa, drawp, nbins, uaxes, zoom, plotlevels, gridlist, gcolor, linstyl = options
    labh = ps.particles_label
    drawh = drawp and nbins > 1
    h5f = h5py.File(pthfilen, 'r')
    time = h5f.attrs['time'][0]
    utim = h5f['dataset_units']['time_unit'].attrs['unit']
    ulenf = h5f['dataset_units']['length_unit'].attrs['unit']
    usc, ulen, uupd = pu.change_units(ulenf, uaxes)
    if drawd:
        uvar = h5f['dataset_units'][var].attrs['unit']
    if drawh:
        umass = h5f['dataset_units']['mass_unit'].attrs['unit']
        labh = labh[:-1] + ' mass histogram' + " [%s]" % pu.labelx()(umass)
    smin = h5f['simulation_parameters'].attrs['domain_left_edge']
    smax = h5f['simulation_parameters'].attrs['domain_right_edge']
    n_d = h5f['simulation_parameters'].attrs['domain_dimensions']
    if uupd:
        smin = pu.list3_division(smin, usc)
        smax = pu.list3_division(smax, usc)
    cgcount = int(h5f['data'].attrs['cg_count'])
    glevels = h5f['grid_level'][:]
    maxglev = max(glevels)

    timep = "time = %5.2f %s" % (time, pu.labelx()(utim))
    print(timep)

    parts, field = [drawp, ], [drawd, ]

    if not cu:
        center = (smax[0] + smin[0]) / 2.0, (smax[1] + smin[1]) / 2.0, (smax[2] + smin[2]) / 2.0

    drawa, drawu = pu.choose_amr_or_uniform(drawa, drawu, drawd, drawg, drawp, maxglev, gridlist)
    plotlevels = pu.check_plotlevels(plotlevels, maxglev, drawa)
    linstyl = pu.linestyles(linstyl, maxglev, plotlevels)
    if drawg:
        gcolor = pu.reorder_gridcolorlist(gcolor, maxglev, plotlevels)
    gridlist = pu.sanitize_gridlist(gridlist, cgcount)

    if drawp:
        pinfile, pxyz, pm = rd.collect_particles(h5f, drawh, center, player, uupd, usc, plotlevels, gridlist)
        parts = pinfile, pxyz, pm, nbins, pcolor, psize, player, labh
        drawh = drawh and pinfile

    draw1D, draw2D, figmode = pu.check_1D2Ddefaults(axc, n_d, drawd and drawh)

    refis = []
    if drawd or drawg:
        extr = [], [], [], [], [], []
        if drawu:
            if len(plotlevels) > 1:
                print('For uniform grid plotting only the firs given level!')
            print('Plotting base level %s' % plotlevels[0])
            refis, extr = rd.reconstruct_uniform(h5f, var, plotlevels[0], gridlist, cu, center, smin, smax, draw1D, draw2D)

        if drawa or drawg:
            refis, extr = rd.collect_gridlevels(h5f, var, refis, extr, maxglev, plotlevels, gridlist, cgcount, center, usc, drawd, draw1D, draw2D)

        if refis == []:
            drawd = False
        else:
            if drawd:
                d1min, d1max, d2min, d2max, d3min, d3max = min(extr[0]), max(extr[1]), min(extr[2]), max(extr[3]), min(extr[4]), max(extr[5])
                vmin, vmax, symmin, autsc = pu.scale_manage(sctype, refis, umin, umax, any(draw2D), [d1min, d1max, d2min, d2max])

                print('3D data value range: ', d3min, d3max)
                if any(draw2D):
                    print('Slices  value range: ', d2min, d2max)
                if any(draw1D):
                    print('1D data value range: ', d1min, d1max)
                print('Plotted value range: ', vmin, vmax)
                vlab = var + " [%s]" % pu.labelx()(uvar)
                field = drawd, vmin, vmax, sctype, symmin, cmap, autsc, vlab

    h5f.close()

    if not (parts[0] or drawd or drawg):
        print('No particles or levels to plot. Skipping.')
        return

    cbar_mode = pu.colorbar_mode(drawd, drawh, figmode)

    if not zoom[0]:
        zoom = False, smin, smax

    equip1d = smin, smax, zoom, ulen, umin, umax, linstyl, output, timep
    equip2d = smin, smax, zoom, ulen, drawg, gcolor, center

    if draw1D[0]:
        plot1d(refis, field, parts, equip1d, 0, 1, 2)
    if draw1D[1]:
        plot1d(refis, field, parts, equip1d, 1, 0, 2)
    if draw1D[2]:
        plot1d(refis, field, parts, equip1d, 2, 0, 1)

    if any(draw2D):
        fig = P.figure(1, figsize=pu.figsizes[figmode])

        grid = AxesGrid(fig, 111, nrows_ncols=pu.figrwcls[figmode], axes_pad=ps.plot2d_axes_pad, aspect=True, cbar_mode=cbar_mode, label_mode="L",)
        ag0, ag2, ag3 = [], [], []

        if draw2D[0]:
            ax = grid[pu.figplace[figmode][0]]
            ax, ag3, ah = draw_plotcomponent(ax, refis, field, parts, equip2d, 0, 1, 2)

        if draw2D[1]:
            ax = grid[pu.figplace[figmode][1]]
            ax, ag2, ah = draw_plotcomponent(ax, refis, field, parts, equip2d, 1, 0, 2)

        if draw2D[2]:
            ax = grid[pu.figplace[figmode][2]]
            ax, ag0, ah = draw_plotcomponent(ax, refis, field, parts, equip2d, 2, 0, 1)
        ax.set_title(timep)

        if drawh:
            add_cbar(figmode, cbar_mode, grid, ah[3], 0, labh)

        if drawd:
            add_cbar(figmode, cbar_mode, grid, pu.take_nonempty([ag0, ag2, ag3]), 1, vlab)

        P.draw()
        out2d = output[0] + pu.plane_in_outputname(figmode, draw2D) + output[1]
        P.savefig(out2d, facecolor=ps.f_facecolor)
        print(out2d, "written to disk")
        P.clf()
