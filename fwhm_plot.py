# -*- coding: utf-8 -*-
"""Plot FWHM for the provided intensity file.

Author: Maksim Rakitin, BNL
Date: 2016-08-24
"""
import os

import matplotlib.pyplot as plt
import numpy as np

PICS_DIR = 'pics'


def get_width(dat_file, y_idx='middle', plot=True, variate_steps=0, show=False, suffix='', axis='y'):
    list_2d, x, y = prepare_data(dat_file)

    if axis == 'x':
        idx = 0
        axis_values = np.copy(x)
        other_values = np.copy(y)
        position = 'Vertical'
    elif axis == 'y':
        idx = 1
        axis_values = np.copy(y)
        other_values = np.copy(x)
        position = 'Horizontal'
    else:
        raise ValueError('"{}" is not a valid option. Use either "x" or "y".'.format(axis))

    if y_idx == 'middle':
        y_idx = np.where(list_2d == list_2d.max())[idx][0]

    factor = 0.5
    m_to_um = 1e6

    widths = []
    for i in range(y_idx - variate_steps, y_idx + variate_steps + 1):
        try:
            y_coor = axis_values[i] * m_to_um
        except IndexError:
            raise ValueError('Index #{} is beyond the size of the array {}'.format(i, len(y)))

        if axis == 'x':
            x_at_zero_y = list_2d[i, :]
        else:
            x_at_zero_y = list_2d[:, i]
        max_value = x_at_zero_y.max()
        frac_of_max = max_value * factor
        idx_ge_frac_of_max = np.where(x_at_zero_y >= frac_of_max)[0]
        x_ge_frac_of_max = other_values[idx_ge_frac_of_max]

        # Rescaling for plotting:
        x_um = np.copy(other_values) * m_to_um
        x_ge_frac_of_max_um = np.copy(x_ge_frac_of_max) * m_to_um

        vals_ge_frac_of_max = x_at_zero_y[idx_ge_frac_of_max]
        width = abs(x_ge_frac_of_max_um[-1] - x_ge_frac_of_max_um[0])  # um
        widths.append(width)

        if plot:
            if not os.path.isdir(PICS_DIR):
                os.mkdir(PICS_DIR)
            units_um = u'\u00B5m'
            frac_of_max_text = '{} * max intensity ({})'.format(factor, max_value)
            plt.figure(figsize=(10, 6))
            plt.title(
                u'Width at {}: {} {}'.format(frac_of_max_text, round(width, 3), units_um))
            plt.xlabel(u'{} Position [{}]'.format(position, units_um))
            plt.ylabel(u'Intensity [ph/s/.1%bw/mm²]')

            int_dat = os.path.join(PICS_DIR, 'intensity_{}_{}.dat'.format(suffix, axis))
            with open(int_dat, 'w') as f:
                for i in range(len(x_um)):
                    f.write('{}\t{}\n'.format(x_um[i], x_at_zero_y[i]))

            plt.plot(x_um, x_at_zero_y, label=u'All X values at Y={} {} (index #{})'.format(y_coor, units_um, i))
            plt.plot(x_ge_frac_of_max_um, vals_ge_frac_of_max, c='red', label='Values >= {}'.format(frac_of_max_text))
            half_width = abs(x_ge_frac_of_max_um[-1] - x_ge_frac_of_max_um[0]) * 0.5
            plt.xlim([x_ge_frac_of_max_um[0] - half_width, x_ge_frac_of_max_um[-1] + half_width])
            plt.grid()
            plt.axvline(0.0, c='black')
            plt.legend()
            if show:
                plt.show()
            else:
                suf = '{}_'.format(suffix) if suffix else ''
                pic_name = os.path.join(PICS_DIR, 'x_{}{:05d}_{}.png'.format(suf, i, axis))
                plt.savefig(pic_name, dpi=100)
            plt.cla()
            plt.clf()
            plt.close()

    return widths


def prepare_data(dat_file):
    list_1d, x_range, y_range = _read_data(dat_file)
    list_2d = _convert_1d_to_2d(list_1d, x_range, y_range)

    x = np.linspace(x_range[0], x_range[1], x_range[2])
    y = np.linspace(y_range[0], y_range[1], y_range[2])

    return list_2d, x, y


def _convert_1d_to_2d(list_1d, x_range, y_range):
    tot_len = int(x_range[2] * y_range[2])
    len_1d = len(list_1d)
    if len_1d > tot_len:
        list_1d = np.array(list_1d[0:tot_len])
    elif len_1d < tot_len:
        aux_list = np.zeros(len_1d)
        for i in range(len_1d):
            aux_list[i] = list_1d[i]
        list_1d = np.array(aux_list)
    list_1d = np.array(list_1d)
    list_2d = list_1d.reshape(x_range[2], y_range[2], order='F')
    return list_2d


def _parse_header(row, data_type):
    return data_type(row.split('#')[1].strip())


def _read_data(dat_file, skip_lines=11):
    list_1d = np.loadtxt(dat_file)
    with open(dat_file, 'r') as f:
        content = f.readlines()[:skip_lines]
        x_range = [
            _parse_header(content[4], float),
            _parse_header(content[5], float),
            _parse_header(content[6], int),
        ]
        y_range = [
            _parse_header(content[7], float),
            _parse_header(content[8], float),
            _parse_header(content[9], int),
        ]
    return list_1d, x_range, y_range


if __name__ == '__main__':
    for suffix in ['es1', 'es2']:
        dat_file = 'SMI/me_fwhm/res_int_pr_me_{}.dat'.format(suffix)
        for axis in ['x', 'y']:
            widths = get_width(dat_file, y_idx='middle', variate_steps=0, show=False, suffix=suffix, axis=axis)
            print('Widths for {} ({}): {}'.format(suffix, axis, widths))
