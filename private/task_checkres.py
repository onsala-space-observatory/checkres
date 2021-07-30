#############################################################################
#
#    CHECK_RES: Quick check of residuals in Fourier space
#
#    Copyright (C) 2016  Ivan Marti-Vidal (Nordic ARC Node, OSO, Sweden)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#############################################################################

import matplotlib.pyplot as pl
import numpy as np
import os
from matplotlib.widgets import Slider
import casatools

ia = casatools.image()
ms = casatools.ms()
tb = casatools.table()

def checkres(vis='', residual=''):
    auxout = checkresAux(vis, residual)
    return auxout


def checkresAux(vis='', residual=''):
    # if True:

    RESIM = residual
    VIS = vis

# FOR TESTING:
#  RESIM = '/home/marti/NO_BACKUP/MERLIN/toto.psf'
#  VIS = '/home/marti/NO_BACKUP/MERLIN/MERLIN.ms'
#  RESIM = '/media/marti/LaCie_1/WORKAREA/ALMA_REDUCED/M100_B3/M100_Band3_CalibratedData/M100line.residual'
#  VIS = '/media/marti/LaCie_1/WORKAREA/ALMA_REDUCED/M100_B3/M100_Band3_CalibratedData/M100all_lores.ms.contsub'

    Const = {'rad': 180./np.pi*3600., 'deg': 3600.}

    print(VIS)

    ms.open(VIS)
    ms.selectinit(datadescid=0)
    uvdat = ms.getdata(['u', 'v', 'time', 'axis_info'], ifraxis=True)
    ms.close()

    tb.open(os.path.join(VIS, 'ANTENNA'))
    antnames = list(tb.getcol('NAME'))
    antcoords = tb.getcol('POSITION')
    tb.close()

    # r = np.sqrt(antcoords[0]**2. + antcoords[1]**2.)
    R = np.sqrt(antcoords[0]**2. + antcoords[1]**2. + antcoords[2]**2.)
    lat = np.arcsin(antcoords[2]/R)*180./np.pi
    avlat = lat - np.average(lat)
    lon = np.arctan2(antcoords[1], antcoords[0])*180./np.pi
    avlon = lon - np.average(lon)

    offmax = max(np.max(np.abs(avlon)), np.max(np.abs(avlat)))

    if offmax < 0.01:
        avlat *= 3600.
        avlon *= 3600.
        aunit = 'as'
        offmax *= 3600.
    else:
        aunit = 'deg'

    ia.open(RESIM)
    imagedat = ia.getchunk()
    peakim = np.max(imagedat)
    arcx = float(ia.summary()['shape'][0])*ia.summary()['incr'][0] / \
        2.*Const[ia.summary()['axisunits'][0]]
    arcy = float(ia.summary()['shape'][1])*ia.summary()['incr'][1] / \
        2.*Const[ia.summary()['axisunits'][1]]
    os.system('rm -rf CHECKRES_TEMP.resfft')
    ia.fft(amp='CHECKRES_TEMP.resfft')
    ia.close()

    ia.open('CHECKRES_TEMP.resfft')
    resids = ia.getchunk()
    peakres = np.max(resids)
    freqs = ia.summary()['refval'][-1]+np.arange(ia.summary()['shape'][-1])*ia.summary()['incr'][-1]
    # maxfreq = np.max(freqs)
    # minfreq = np.min(freqs)
    lambdasx = np.arange(
        int(ia.summary()['shape'][0]/2.))*ia.summary()['incr'][0]
    lambdasy = np.arange(
        int(ia.summary()['shape'][1]/2.))*ia.summary()['incr'][1]

    nchan = len(freqs)
    ia.close()

    basname = uvdat['axis_info']['ifr_axis']['ifr_name']
    basants = [filter(lambda x: antnam in x, basname) for antnam in antnames]

    hours = 24.*(uvdat['time']/86400. - np.floor(uvdat['time']/86400.))

    lambdasx *= 3.e8/np.average(freqs)
    lambdasy *= 3.e8/np.average(freqs)

    basmax = max(np.max(np.abs(uvdat['u'])), np.max(
        np.abs(uvdat['v'])))  # *3.e8/np.average(freqs)
    imbasmax = max(np.max(lambdasx), np.max(np.abs(lambdasy)))

    if basmax > 100.:
        basmax /= 1000.
        imbasmax /= 1000.
        lambdasx /= 1000.
        lambdasy /= 1000.
        uvdat['u'] /= 1000.
        uvdat['v'] /= 1000.
        lunit = 'km'
    else:
        lunit = 'm'

    fig = pl.gcf()  # pl.figure(figsize=(8,8))
    fig.clf()
    # sub1 = fig.add_axes([0.07,0.2,0.35,0.7],aspect='equal')
    # sub2 = fig.add_axes([0.65,0.2,0.35,0.7],aspect='equal',sharex=sub1,sharey=sub1)
    # sub3 = fig.add_axes([0.47,0.65,0.15,0.3],aspect='equal')
    # sub4 = fig.add_axes([0.47,0.25,0.15,0.3],aspect='equal')
    sub1 = fig.add_subplot(221, aspect='equal')
    sub2 = fig.add_subplot(222, aspect='equal', sharex=sub1, sharey=sub1)
    sub3 = fig.add_subplot(223, aspect='equal')
    sub4 = fig.add_subplot(224, aspect='equal')
    fig.subplots_adjust(bottom=0.2)

    slax = fig.add_axes([0.2, 0.05, 0.5, 0.05])
    slid = Slider(slax, 'Channel', 0, nchan, valinit=int(nchan/2))
    # but1 = fig.add_axes([0.9,0.05,0.075,0.05])
    # plus = Button(but1,'+ chan')
    # but2 = fig.add_axes([0.8,0.05,0.075,0.05])
    # minus = Button(but2,'- chan')

    canvas = fig.canvas

    basplot = []
    basplot2 = []
    for bi, bas in enumerate(basname):
        mask = np.logical_or(uvdat['u'][bi, :] != 0.0,
                             uvdat['v'][bi, :] != 0.0)
        basplot.append(sub1.plot(uvdat['u'][bi, mask], uvdat['v'][bi, mask],
                                 '.r', picker=3, markersize=2, markeredgewidth=0.0)[0])
        basplot2.append(sub1.plot(-uvdat['u'][bi, mask], -uvdat['v'][bi, mask],
                                  '.b', picker=3, markersize=2, markeredgewidth=0.0)[0])

    basplotaux = sub2.plot([], [], 'ow', markersize=2, markeredgewidth=0.0)[0]
    antplotaux = sub4.plot([0], [0], '-b')[0]
    basplotaux2 = sub1.plot([], [], 'xk', markersize=6)[0]

    text = sub1.text(0.05, 1.02, "Click on any UV point",
                     transform=sub1.transAxes)
    text2 = sub2.text(0.05, 1.02, " ", transform=sub2.transAxes)

    sub1.set_xlabel('U (%s)' % lunit)
    sub1.set_ylabel('V (%s)' % lunit)

    if len(np.shape(resids)) < 4:
        resids = resids[:, np.newaxis]
        imagedat = imagedat[:, np.newaxis]
        nchan = 1

    resim = sub2.imshow(np.transpose(resids[:, :, 0, int(nchan/2)]),
                                     extent=(-imbasmax, imbasmax, -imbasmax, imbasmax),
                                             vmax=peakres, vmin=0.0,
                                             interpolation='gaussian')
    image = sub3.imshow(np.transpose(imagedat[:, :, 0, int(nchan/2)]), extent=(
        arcx, -arcx, -arcy, arcy), vmax=peakim, vmin=0.0, interpolation='gaussian', origin='lower')

    # sub2.set_xlabel('U (%s)'%lunit)
    # sub2.set_ylabel('V (%s)'%lunit)
    pl.setp(sub2.get_xticklabels(), visible=False)
    pl.setp(sub2.get_yticklabels(), visible=False)

    sub3.set_xlabel('RA offset (as)')
    sub3.set_ylabel('Dec offset (as)')

    sub4.plot(avlon, avlat, 'or', picker=3)

    sub4.set_xlabel('Longitude offset (%s)' % aunit)
    sub4.set_ylabel('Latitude offset (%s)' % aunit)
    for i, nam in enumerate(antnames):
        sub4.annotate(nam, xy=(
            avlon[i], avlat[i]), xytext=(-5, 5), textcoords='offset points', fontsize=10)

    sub4.set_xlim((-offmax*1.1, offmax*1.1))
    sub4.set_ylim((-offmax*1.1, offmax*1.1))

    sub2.set_xlim((basmax*1.1, -basmax*1.1))
    sub2.set_ylim((-basmax*1.1, basmax*1.1))

    def _onPick(event):

        if event.mouseevent.inaxes == sub4:

            iant = event.ind[0]

        elif event.mouseevent.inaxes == sub1:
            if event.artist in basplot:
                selbas = basplot.index(event.artist)
            else:
                selbas = basplot2.index(event.artist)
            if len(event.ind) > 1:
                idx = event.ind[0]
            else:
                idx = int(event.ind)
            text.set_text('%s at %.3f hours' % (basname[selbas], hours[idx]))

            selants = filter(lambda x: x in basname[selbas], antnames)
            if event.mouseevent.button == 1:
                iant = antnames.index(selants[0])
                iant2 = antnames.index(selants[1])
            else:
                iant = antnames.index(selants[1])
                iant2 = antnames.index(selants[0])

            antplotaux.set_data(
                ([avlon[iant], avlon[iant2]], [avlat[iant], avlat[iant2]]))
            basplotaux2.set_data(([event.artist.get_data()[0][idx]], [
                                 event.artist.get_data()[1][idx]]))

        else:
            return

        alldataX = []
        alldataY = []
        text2.set_text('Highlight baselines to %s' % antnames[iant])
        # print iant, basants[iant]
        for bas in basants[iant]:
            ibas = np.where(basname == bas)[0][0]
            # print ibas
            X, Y = basplot[ibas].get_data()
            alldataX.append(X)
            alldataY.append(Y)
            X, Y = basplot2[ibas].get_data()
            alldataX.append(X)
            alldataY.append(Y)
        basplotaux.set_data(
            (np.concatenate(alldataX), np.concatenate(alldataY)))

        pl.draw()

    def _onSlide(channel):
        if channel != int(channel):
            slid.set_val(int(channel))
            return
        resim.set_data(np.transpose(resids[:, :, 0, int(channel)]))
        image.set_data(np.transpose(imagedat[:, :, 0, int(channel)]))

        pl.draw()

    def _addchan(event):
        currchan = slid.val
        if currchan + 1 < nchan:
            slid.set_val(currchan+1)

    def _subchan(event):
        currchan = slid.val
        if currchan > 0:
            slid.set_val(currchan-1)

    slid.on_changed(_onSlide)
    canvas.mpl_connect('pick_event', _onPick)
    # plus.on_clicked(_addchan)
    # minus.on_clicked(_subchan)

    pl.show()

    return slid  # , plus, minus, but1, but2
