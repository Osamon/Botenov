# -*- coding: utf-8 -*-

"""
Зная БГК, расчитываем те матрицы, которые можем получить сверткой
"""

import numpy as np

from mc2py.zgiw import getv, GetHgrp3, GetG, IcamRD, GetTin, GetTout, GetDp

from akpm1200.data.koefs_code.bgk_koeffs import bgk_koeffs


def get_koeffs(bgk, erheat, wnom, plant='kaln', Hnums=[-3,-2,-1], Icam0=None):
    """Получение коэффициентов из текущего состояния"""
    v = getv()
    koeff = bgk_koeffs(bgk)
    koeff["t"] = float(v.YMTIME)
    koeff["h"] = GetHgrp3(Hnums)
    koeff["G"] = GetG()
    koeff["Tin"] = GetTin()
    koeff["Tout"] = GetTout()
    koeff["dP"] = GetDp()
    koeff['offset'] = float(v.YMOFFSET)
    koeff["wt"] = (1 - np.sum(v.YMRP_POW[...]) / wnom) * (np.sum(v.YMFISF[...]) / wnom) * 100
    if Icam0 is not None:
        koeff["I"] = IcamRD(plant) / Icam0
    else:
        koeff["I"] = IcamRD(plant)
    return koeff


def make_matrix(dims, bgk, inipoint, defaults, norm_koef):
    """
    Расчет свернутых операторов  для перевода одних переменных в другие
    (учитываются заявленные размерности базисов)
    """
    res_dict = {}
    erheat = defaults['erheat'][-1]

    wnom = defaults['Wnom'][-1]
    w_norm = norm_koef['w_norm'][-1]
    xe_norm = norm_koef['xe_norm'][-1]
    id_norm = norm_koef['id_norm'][-1]

    pairs = [["Id", "Xe"], ["w", "Id"]]
    for nmv1, nmv2 in pairs:
        n1 = dims["N" + nmv1]
        n2 = dims["N" + nmv2]
        v1 = bgk[nmv1][:n1]
        v2 = bgk[nmv2][:n2]
        res = np.dot(v2, np.transpose(v1))
        outname = nmv1 + "2" + nmv2
        res_dict[outname] = [outname, outname, res]

    # тензор для смеси мощности и ксенона
    Xeb = np.array(bgk["Xe"][:dims["NXe"]])
    wb = np.array(bgk["w"][:dims["Nw"]])

    #sec = inipoint['p3d'] / w_norm / wnom * erheat * inipoint['xe_init'] / xe_norm / (-inipoint['xe_init'] / xe_norm
    #                                                                                  + inipoint['id_init'] / id_norm)
    sec = np.array(inipoint['p3d']) / w_norm / wnom * erheat * np.array(inipoint['xe_init']) / xe_norm / (-np.array(inipoint['xe_init'] ) / xe_norm
                                                                                      + np.array(inipoint['id_init']) / id_norm)

    nXeXew = (dims["NXe"], dims["NXe"], dims["Nw"])
    wXe2Xe = np.zeros(nXeXew, dtype='f')
    for (k0, k1, k2, i) in np.ndindex(nXeXew + (len(sec),)):
        wXe2Xe[k0, k1, k2] += Xeb[k0, i] * Xeb[k1, i] * wb[k2, i] / sec[i]
    res_dict['wXe2Xe'] = ['wXe2Xe', 'wXe2Xe', wXe2Xe]

    bgkWn = wb.reshape((-1, dims['NTVS'], dims['NZ']))

    w2Kz = np.transpose(np.sum(bgkWn, 1))

    if dims['NZ'] == 20:
        buf = []
        for i in range(10):
            buf.append(w2Kz[i*2]+w2Kz[i*2+1])
        w2Kz = np.array(buf)

    res_dict['w2Kz'] = ['w2Kz', 'w2Kz', w2Kz]

    return res_dict


if __name__ == '__main__':
    from mc2py.scutil import ld
    dims = ld(r'F:\Documents\Projects\akpm1200\data\data\kaln\b03\k01\out\basedims.json')
    bgk = ld(r'F:\Documents\Projects\akpm1200\data\data\kaln\b03\k01\out\bgk.pkl')
    inipoint = ld(r'F:\Documents\Projects\akpm1200\data\data\kaln\b03\k01\out\initialpoint.pkl')
    defaults = ld(r'F:\Documents\Projects\akpm1200\data\data\kaln\b03\k01\out\defaults.pkl')
    norm_koef = ld(r'F:\Documents\Projects\akpm1200\data\data\kaln\b03\k01\out\norm_koeffs.pkl')
    res_dict = make_matrix(dims, bgk, inipoint, defaults, norm_koef)
