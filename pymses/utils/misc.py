# -*- coding: utf-8 -*-
#   This file is part of PyMSES.
#
#   PyMSES is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PyMSES is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.
"""
miscellaneous utilities module
"""

import numpy


def deleteObject(obj):
    del obj


def concatenate_reorder(array_list, final_indices_list=None, axis=0):
    r"""
    Concatenates a list of numpy arrays over a given `axis` and reorders the final array according to
    `final_indices_list`.

    Parameters
    ----------

    array_list : ``list`` of ``array``
        arrays to be concatenated
    final_indices_list : ``list`` of ``array`` of ``int``
        For each array `arr` from `array_list` and its corresponding `idx` from	`final_indices_list`, one has::

        >>> cat_arrays[idx[i]] = arr[i]

    Returns
    -------

    cat_arrays : ``array``
        concatenation of the elements of `array_list`. reordered if `final_indices_list` is not None,

    """

    # Concatenate the arrays
    cat_arrays = numpy.concatenate(array_list, axis=axis)
    out_shape = cat_arrays.shape

    if final_indices_list is None:
        return cat_arrays

    slc = [slice(n) for n in out_shape]

    for arr, idx in zip(array_list, final_indices_list):
        slc[axis] = idx
        cat_arrays[tuple(slc)] = arr

    return cat_arrays


def fname_iterator(path, file_regexp):
    r"""
    Returns a filename iterator function yielding every filename matching
    a given pattern of files found in a given directory.

    Parameters
    ----------
    path        : ``string``
        path of the directory to search into
    file_regexp : ``string``
        regular expression of the filename to be matched

    Returns
    -------
    fname_iter : iterator ``function`` yielding ``string`` filepath
        iterates over files matching the search pattern in the search directory.

    """
    import os
    import re

    assert (os.path.isdir(path)), "'%s' directory doesn't exist" % path
    regexp = re.compile(file_regexp)
    file_list = []
    ls = os.listdir(path)
    for file in ls:
        res = regexp.findall(file)
        if len(res) > 0:
            if res[0] == file:
                file_list.append(res[0])

    assert (len(file_list) != 0), \
        "No file found matching '%s' in '%s' directory." % (file_regexp, path)

    def fname_iterator():
        for fname in file_list:
            yield os.path.join(path, fname)

    return fname_iterator


def chunks(cpu_list, node_number):
    chunked_cpu_list = []
    j = 0
    init = True
    for cpu in cpu_list:
        if init:
            chunked_cpu_list.append([])
        chunked_cpu_list[j].append(cpu)
        j += 1
        if j == node_number:
            j = 0
            init = False
    return chunked_cpu_list


def balanced_cpu_list(cpu_cost_list, node_number):
    """
    Build a cost balanced cpu list

    Parameters
    ----------
    cpu_cost_list : list of cpu with associated cost structure [cost, cpu_number]
    node_number : length of the balanced cpu list returned

    Example
    --------
        In  [1]: balanced_cpu_list( [[0, 3],[46, 1],[13, 2], [63, 5],[96, 4]], 3)
        Out [1]: [[4], [5], [1, 2, 3]]

    """
    balanced_cpu_cost_list = []
    cpu_cost_list.sort()
    init = True
    while len(cpu_cost_list) > 0:
        if init:
            for i in range(min(node_number, len(cpu_cost_list))):
                balanced_cpu_cost_list.append([0, []])
                cost, cpu = cpu_cost_list.pop()
                balanced_cpu_cost_list[i][0] += cost
                balanced_cpu_cost_list[i][1].append(cpu)
        else:
            cost, cpu = cpu_cost_list.pop()
            balanced_cpu_cost_list[0][0] += cost
            balanced_cpu_cost_list[0][1].append(cpu)
        balanced_cpu_cost_list.sort()
        init = False
    # print "balanced_cpu_cost_list",balanced_cpu_cost_list
    balanced_cpu_list = []
    for i in range(len(balanced_cpu_cost_list)):
        balanced_cpu_list.append(balanced_cpu_cost_list.pop()[1])
    return balanced_cpu_list


__all__ = ["concatenate_reorder", "fname_iterator", "chunks", "balanced_cpu_list"]
