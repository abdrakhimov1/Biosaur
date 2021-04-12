from copy import copy
from collections import defaultdict
import numpy as np
from scipy.signal import medfilt
import math
import itertools


def meanfilt(data, window_width):
    cumsum_vec = np.cumsum(np.insert(data, 0, 0))
    ma_vec = (cumsum_vec[window_width:] -
              cumsum_vec[:-window_width]) / window_width
    ma_vec = data[:1] + list(ma_vec) + data[-1:]
    return ma_vec


class ready_hill:

    def __init__(self, intensity, scan_id, mass, ion_mobility):
        # self.mz = np.median(mass)
        # self.mz = np.mean(mass)
        # self.mz_std = np.std(mass)
        self.intensity = intensity
        self.scan_id = scan_id
        self.scan_set = set(scan_id)
        self.mass = mass
        self.diff_for_output = 0
        tmp = max(range(len(self.intensity)), key=self.intensity.__getitem__)
        self.scan_of_max_intensity = self.scan_id[tmp]
        self.max_intensity = self.intensity[tmp]
        # self.mz = np.average(self.mass, weights=self.intensity)
        # self.mz = sum(weight * value for weight, value in zip(self.intensity, self.mass)) / sum(self.intensity)
        # self.mz = self.mass[tmp]
        # self.max_intensity = sum(self.intensity)
        if not (ion_mobility is None):
            self.ion_mobility = ion_mobility
            self.opt_ion_mobility = self.ion_mobility[tmp]
        else:
            self.ion_mobility = None
            self.opt_ion_mobility = None
        self.scan_len = len(self.scan_id)

        self.idict = dict()
        for i, j in zip(self.scan_id, self.intensity):
            self.idict[i] = j

        # self.sqrt_of_i_sum_squares = math.sqrt(
        #     sum(v**2 for v in self.idict.values()))
        intensity_np = np.array(intensity)
        self.sqrt_of_i_sum_squares = np.sqrt(np.sum(np.power(intensity_np, 2)))

class next_peak:

    def __init__(
            self,
            next_mz_array,
            next_intensity_array,
            next_scan_id,
            next_ion_mobility_array):

        self.next_mz_array = next_mz_array
        self.next_intensity_array = next_intensity_array
        self.next_ion_mobility_array = next_ion_mobility_array
        self.next_scan_id = next_scan_id


class peak_ion_mobility:

    def __init__(self, mz, intensity, ion_mobility):
        self.mz_array = [mz, ]
        self.mass_array = [[mz, ], ]
        self.intensity_array = [[intensity, ]]
        self.ion_mobility_array = [[ion_mobility, ]]
        self.intensity_max = [intensity, ]
        self.ion_mobility_opt = [ion_mobility, ]
        self.ion_mobility_max = [ion_mobility, ]
        self.ion_mobility_min = [ion_mobility, ]
        self.total = 1

    def get_nearest_values(self, value):
        return np.argsort(np.abs(self.mz_array) - value)

    def extend(self, mz, intensity, ion_mobility):
        self.mz_array.append(mz)
        self.mass_array.append([mz, ])
        self.intensity_array.append([intensity, ])
        self.intensity_max.append(intensity)
        self.ion_mobility_opt.append(ion_mobility)
        self.ion_mobility_array.append([ion_mobility, ])
        self.ion_mobility_max.append(ion_mobility)
        self.ion_mobility_min.append(ion_mobility)
        self.total += 1

    def append_and_recalc(self, mz, intensity, ion_mobility, index):
        self.mass_array[index].append(mz)
        self.intensity_array[index].append(intensity)
        self.ion_mobility_array[index].append(ion_mobility)
        self.recalc(index)

    def recalc(self, index):
        self.mz_array[index] = np.mean(self.mass_array[index])
        self.ion_mobility_max[index] = max(self.ion_mobility_array[index])
        self.ion_mobility_min[index] = min(self.ion_mobility_array[index])
        # if self.intensity_array[index][-1] > self.intensity_max[index]:
        #     # self.mz_array[index] = 
        #     self.intensity_max[index] = self.intensity_array[index][-1]
        #     self.ion_mobility_opt[index] = self.ion_mobility_array[index][-1]

    def push_me_to_the_peak_ion_mob(self, mz, intensity, ion_mobility, diff):
        # nearest_ids = self.get_nearest_values(mz)
        flag = 0

        nearest_id = self.total - 1
        mass_accuracy = diff * 1e-6 * mz
        while nearest_id >= 0:
            tmp_diff = abs(self.mz_array[nearest_id] - mz)
            # tmp_diff = abs(self.mz_array[nearest_id] - mz) / mz
            # if tmp_diff <= diff * 1e-6:
            if tmp_diff <= mass_accuracy:
                if abs(
                        self.ion_mobility_max[nearest_id] -
                        ion_mobility) <= 0.05 or abs(
                        self.ion_mobility_min[nearest_id] -
                        ion_mobility) <= 0.05:
                    flag = 1
                    self.append_and_recalc(
                        mz, intensity, ion_mobility, nearest_id)
                    break
            else:
                break
            nearest_id -= 1

        if not flag:
            self.extend(mz, intensity, ion_mobility)


class peak:
    def __init__(
            self,
            mz_array,
            intensity,
            scan_id,
            start_id,
            ion_mobility_array):

        self.mz_array = copy(mz_array)

        self.scan_id = [[scan_id, ] for _ in range(len(mz_array))]
        # self.scan_id = []

        # for _ in range(len(mz_array)):
        #     self.scan_id.append([scan_id, ])

        self.intensity = [[i, ] for i in intensity]
        if not (ion_mobility_array is None):
            self.ion_mobility = [[i, ] for i in ion_mobility_array]
        else:
            self.ion_mobility = None
        # self.intensity = []
        # for i in intensity:
        #     self.intensity.append([i, ])

        self.mass_array = [[i, ] for i in mz_array]
        # self.mass_array = []
        # for i in mz_array:
        #     self.mass_array.append([i, ])

        self.finished_hills = []
        self.crosslinked_hills = []

        self.intervals = [start_id, ]
        self.actual_degree = 0

        self.medar = [1.0, ]


    def get_potential_isotope_id(self, i_fast, i_idx):
        tmp = self.finished_hills_fast_dict.get(i_fast, [])
        # tmp.remove(i_idx)
        return tmp

    def recalc_fast_array_for_finished_hills(self, mz_step):
        m_koef = mz_step
        im_koef = 0.02
        self.finished_hills_fast_array = [int(fh.mz/m_koef) for fh in self.finished_hills]
        self.finished_hills_fast_dict = defaultdict(set)
        for idx, fm in enumerate(self.finished_hills_fast_array):
            self.finished_hills_fast_dict[fm-1].add(idx)
            self.finished_hills_fast_dict[fm+1].add(idx)
            self.finished_hills_fast_dict[fm].add(idx)

    def recalc_fast_array(self, mz_step):
        m_koef = mz_step
        im_koef = 0.02
        # self.fast_array = [int(tm/m_koef) for tm in self.mz_array]
        self.fast_array = (self.mz_array/m_koef).astype(int)
        self.fast_dict = defaultdict(set)
        for idx, fm in enumerate(self.fast_array):
            self.fast_dict[fm-1].add(idx)
            self.fast_dict[fm+1].add(idx)
            self.fast_dict[fm].add(idx)

    def concat_peak_with(self, second_peak):

        self.mz_array = self.mz_array + second_peak.mz_array
        self.intensity = self.intensity + second_peak.intensity
        if not (self.ion_mobility is None):
            self.ion_mobility = self.ion_mobility + second_peak.ion_mobility
        self.mass_array = self.mass_array + second_peak.mass_array
        self.finished_hills = self.finished_hills + second_peak.finished_hills
        self.crosslinked_hills = self.crosslinked_hills + \
            second_peak.crosslinked_hills
        self.intervals = self.intervals + second_peak.intervals

    def crosslink_simple(self, mass_accuracy):

        mz_step = mass_accuracy * 1e-6 * 2500

        crosslink_counter = 0
        self.finished_hills = sorted(
            self.finished_hills,
            key=lambda x: x.scan_id[0])

        allowed_ids = set()
        for i in self.intervals:
            allowed_ids.add(i - 1)
            allowed_ids.add(i - 2)

            
        allowed_ids2 = set()
        for i in self.intervals:
            allowed_ids2.add(i)
            allowed_ids2.add(i+1)

        map_ids_1 = defaultdict(list)
        map_ids_2 = defaultdict(set)

        self.finished_hills_fast_dict = defaultdict(set)
        m_koef = mz_step

        for i, hill in enumerate(self.finished_hills):

            end_scan = hill.scan_id[-1]
            if end_scan in allowed_ids:
                map_ids_1[end_scan].append(i)
                fm = int(hill.mz / m_koef)
                self.finished_hills_fast_dict[fm-1].add(i)
                self.finished_hills_fast_dict[fm+1].add(i)
                self.finished_hills_fast_dict[fm].add(i)

                
            start_scan = hill.scan_id[0]
            if start_scan in allowed_ids2:
                map_ids_2[start_scan].add(i)
                fm = int(hill.mz / m_koef)
                self.finished_hills_fast_dict[fm-1].add(i)
                self.finished_hills_fast_dict[fm+1].add(i)
                self.finished_hills_fast_dict[fm].add(i)

        banned_ids = set()
        way_to_combine = []

        for al_id in sorted(allowed_ids):

            for i in map_ids_1[al_id]:

                if i not in banned_ids:

                    hill = self.finished_hills[i]
                    fm = int(hill.mz / m_koef)
                    for j in self.finished_hills_fast_dict[fm]:

                        if (j in map_ids_2[al_id+1] or j in map_ids_2[al_id+2]) and j not in banned_ids:

                            hill2 = self.finished_hills[j]
                            if abs(hill.mz - hill2.mz) / \
                                    hill.mz <= mass_accuracy * 1e-6:

                                banned_ids.add(i)
                                banned_ids.add(j)
                                way_to_combine.append((i, j))

        for i, j in way_to_combine[::-1]:
            self.finished_hills[i] = ready_hill(
                intensity=hill.intensity +
                hill2.intensity,
                scan_id=hill.scan_id +
                hill2.scan_id,
                mass=hill.mass +
                hill2.mass,
                ion_mobility=(
                    hill.ion_mobility +
                    hill2.ion_mobility
                    if not (hill.ion_mobility is None)
                    else None))
            del self.finished_hills[j]

        for i in list(range(len(self.finished_hills)))[::-1]:
            if len(self.finished_hills[i].scan_id) < 3:
                del self.finished_hills[i]

    def crosslink(self, mass_accuracy):

        crosslink_counter = 0
        # crosslink_counter2 = 0
        self.finished_hills = sorted(
            self.finished_hills,
            key=lambda x: x.scan_id[0])

        i = 0
        ini_len = len(self.finished_hills)

        while i < ini_len:

            hill = self.finished_hills[i]
            j = i + 1

            while j < ini_len:

                hill2 = self.finished_hills[j]

                # if hill.scan_id[-1] == hill2.scan_id[0]:
                if abs(hill.scan_id[-1] - hill2.scan_id[0]) <= 1:
                    # crosslink_counter2 += 1
                    if abs(hill.mz - hill2.mz) / \
                            hill.mz <= mass_accuracy * 1e-6:

                        self.finished_hills[i] = ready_hill(
                            intensity=hill.intensity + hill2.intensity,
                            scan_id=hill.scan_id + hill2.scan_id,
                            mass=hill.mass + hill2.mass,
                            ion_mobility=hill.ion_mobility +
                            hill2.ion_mobility)
                        del self.finished_hills[j]
                        ini_len -= 1
                        crosslink_counter += 1
                elif hill2.scan_id[0] > hill.scan_id[-1] + 1:
                    break

                j += 1

            i += 1

        # print(crosslink_counter)
        # print(crosslink_counter2)

    def sort_finished_hills(self):
        self.finished_hills = sorted(self.finished_hills, key=lambda x: x.mz)

    def check_its_ready(self, id_real, check_degree, min_length):


        # ar_for_median = []
        # for m_ar, scan_ar in zip(self.mass_array, self.scan_id):
        #     if scan_ar[-1] == id_real - 1:
        #         if len(m_ar) >= 2:
        #             ar_for_median.append(m_ar[-1]/m_ar[-2])
        # # print(np.median(ar_for_median), 'median!')
        # if len(ar_for_median) >= 20:
        #     self.medar.append(np.median(ar_for_median))
        # else:
        #     self.medar.append(1.0)

        mask_to_del = [True] * self.mz_array.size
        set_to_del = set()
        for i in range(self.mz_array.size)[::-1]:

            # degree_actual = id_real - self.scan_id[i][0] - len(self.scan_id[i]) + 1
            degree_actual = id_real - self.scan_id[i][-1]
            # or (degree_actual == 2 and len(self.scan_id[i]) == 1):
            if degree_actual > check_degree:

                # degree_actual = id_real - self.scan_id[i][-1]
                # if degree_actual > check_degree or (degree_actual == 2 and
                # len(self.scan_id[i]) <= 3):

                # list_intensity = self.intensity.pop(i)
                list_intensity = self.intensity[i]
                if not (self.ion_mobility is None):
                    # list_ion_mobility = self.ion_mobility.pop(i)
                    list_ion_mobility = self.ion_mobility[i]
                else:
                    list_ion_mobility = None
                # list_scan_id = self.scan_id.pop(i)
                list_scan_id = self.scan_id[i]
                # list_mass = self.mass_array.pop(i)
                list_mass = self.mass_array[i]
                lsi = len(list_scan_id)
                if lsi >= min_length:
                    tmp_ready_hill = ready_hill(intensity=list_intensity,
                                                scan_id=list_scan_id,
                                                mass=list_mass,
                                                ion_mobility=list_ion_mobility,
                                                )
                    self.finished_hills.append(tmp_ready_hill)

                mask_to_del[i] = False
                set_to_del.add(i)

                # if len(tmp_ready_hill.scan_id) >= min_length:
                #     self.finished_hills.append(tmp_ready_hill)

        self.intensity = [i for j, i in enumerate(self.intensity) if j not in set_to_del]
        
        self.scan_id = [i for j, i in enumerate(self.scan_id) if j not in set_to_del]
        self.mass_array = [i for j, i in enumerate(self.mass_array) if j not in set_to_del]
        if not (self.ion_mobility is None):
            self.ion_mobility = [i for j, i in enumerate(self.ion_mobility) if j not in set_to_del]

        self.mz_array = self.mz_array[mask_to_del]

    def push_left(self, min_length):
        mask_to_del = [True] * self.mz_array.size
        for i in range(self.mz_array.size)[::-1]:

            tmp_ready_hill = ready_hill(
                intensity=self.intensity.pop(i),
                scan_id=self.scan_id.pop(i),
                mass=self.mass_array.pop(i),
                ion_mobility=(
                    self.ion_mobility.pop(i) if not (
                        self.ion_mobility is None) else None),
            )
            mask_to_del[i] = False

            if len(tmp_ready_hill.scan_id) >= min_length:
                self.finished_hills.append(tmp_ready_hill)

        self.mz_array = self.mz_array[mask_to_del]
        
        # self.medar.append(1.0)

    def get_nearest_value(self, value, mask):
        return np.argmin(np.abs(self.mz_array[mask] - value))

    def newid(self, nearest, mask):
        return np.nonzero(mask)[0][nearest]

    def get_potential_nearest(self, i_fast):
        return self.fast_dict.get(i_fast, None)

    def get_nearest_id(self, i, prev_nearest, diff, mz_array_l, ion_mobility, mask, mz_step):
        mass_diff = diff * 1e-6 * i
        best_diff = 2 * mass_diff
        best_id = False
        cur_md_abs = 0
        best_prev_nearest_id = False

        i_fast = int(i / mz_step)

        set_idx = self.get_potential_nearest(i_fast)

        if set_idx:
            for nearest_id in set_idx:
                if mask[nearest_id]:
            # nearest_id = prev_nearest
            # while nearest_id < mz_array_l:
                    cur_md = self.mz_array[nearest_id] - i
                    cur_md_abs = abs(cur_md)
                    if cur_md_abs <= mass_diff:
                        if not best_prev_nearest_id:
                            best_prev_nearest_id = int(nearest_id)
                        if (ion_mobility is None) or \
                                abs(ion_mobility -
                                    self.ion_mobility[nearest_id][-1]) <= 0.1:
                            if cur_md_abs <= best_diff:
                                best_diff = float(cur_md_abs)
                                best_id = int(nearest_id)
                        # prev_nearest = int(nearest_id)
                # elif cur_md > mass_diff:
                #     break

                # nearest_id += 1
        if not best_prev_nearest_id:
            best_prev_nearest_id = prev_nearest
        return best_id, best_diff / i, best_prev_nearest_id

    def get_arrays(self, tmp1):
        tmp1_nearest_id_arr = np.array([x[0] for x in tmp1])
        tmp1_idx_arr = np.array([x[1] for x in tmp1])
        tmp1_diff_arr = np.array([x[2] for x in tmp1])
        return tmp1_nearest_id_arr, tmp1_idx_arr, tmp1_diff_arr

    def push_me_to_the_peak(self, next_peak, diff, min_length, mz_step):

        next_mz_array = next_peak.next_mz_array
        next_intensity_array = next_peak.next_intensity_array
        next_ion_mobility_array = next_peak.next_ion_mobility_array
        next_scan_id = next_peak.next_scan_id

        self.check_its_ready(
            id_real=next_scan_id,
            check_degree=2,
            min_length=min_length)

        mask = [True] * (len(self.mz_array))
        tmp1 = []
        tmp2 = []

        prev_nearest = 0

        self.recalc_fast_array(mz_step)
        
        mask = [True] * (len(self.mz_array))

        mz_array_l = len(self.mz_array)
        for idx, i in enumerate(next_mz_array):
            best_id, \
                md_res, \
                prev_nearest = self.get_nearest_id(
                    i,
                    prev_nearest,
                    diff,
                    mz_array_l,
                    (next_ion_mobility_array[idx]
                        if not (
                        next_ion_mobility_array is None)
                        else None), mask, mz_step)
            if best_id:
                tmp1.append([best_id, idx, md_res])

        tmp1_nearest_id_arr, tmp1_idx_arr, tmp1_diff_arr = self.get_arrays(
            tmp1)

        sort_list = np.argsort(tmp1_diff_arr)  # try different kinds
        tmp1_nearest_id_arr = tmp1_nearest_id_arr[sort_list]
        tmp1_idx_arr = tmp1_idx_arr[sort_list]
        tmp1_diff_arr = tmp1_diff_arr[sort_list]

        saved_index = set()

        while tmp1:

            # tmp_id = tmp1_idx_arr[0]

            if tmp1_diff_arr.size == 0:
                break

            if tmp1_diff_arr[0] > diff * 1e-6:
                break

            tmp2.append((tmp1_nearest_id_arr[0], tmp1_idx_arr[0]))

            saved_index.add(tmp1_idx_arr[0])

            mask[tmp2[-1][0]] = False
            if any(mask):
                tmp1_nearest_id_arr = tmp1_nearest_id_arr[1:]

                tmp1_idx_arr = tmp1_idx_arr[1:]

                tmp1_diff_arr = tmp1_diff_arr[1:]

                if tmp1_diff_arr.size == 0:
                    break

                if tmp1_nearest_id_arr[0] in saved_index:

                    for idx, element in enumerate(tmp1_idx_arr):

                        if tmp1_nearest_id_arr[idx] in saved_index:

                            element_mz = next_mz_array[element]
                            element_im = (next_ion_mobility_array[element]
                                    if not (
                                    next_ion_mobility_array is None)
                                    else None)


                            # nearest = self.get_nearest_value(element_mz, mask)
                            # nearest_id_old = self.newid(nearest, mask)

                            nearest_id, \
                                md_res, \
                                prev_nearest = self.get_nearest_id(
                                    element_mz,
                                    0,
                                    diff,
                                    0,
                                    element_im, mask, mz_step)

                            if not nearest_id:
                                nearest_id = 0
                                md_res = 1e6

                            tmp1_nearest_id_arr[idx] = nearest_id

                            tmp1_diff_arr[idx] = md_res
                        else:
                            break
                    sort_list = np.argsort(
                        tmp1_diff_arr, kind='quicksort')  # try different kinds
                    tmp1_nearest_id_arr = tmp1_nearest_id_arr[sort_list]
                    tmp1_idx_arr = tmp1_idx_arr[sort_list]
                    tmp1_diff_arr = tmp1_diff_arr[sort_list]

            else:
                break

        for i, idx in tmp2:
            # FIXME
            # self.mz_array[i] = (self.mz_array[i] + next_mz_array[idx])/2
            self.scan_id[i].append(next_scan_id)
            self.intensity[i].append(next_intensity_array[idx])
            if not (self.ion_mobility is None):
                self.ion_mobility[i].append(next_ion_mobility_array[idx])
            self.mass_array[i].append(next_mz_array[idx])
            tmp_mass_array = self.mass_array[i][-3:]
            self.mz_array[i] = sum(tmp_mass_array)/len(tmp_mass_array)
            # self.mz_array[i] = np.average(self.mass_array[i][-3:], weights=self.intensity[i][-3:])

        added = set(x[1] for x in tmp2)
        mask2 = [(False if i in added else True)
                 for i in range(len(next_mz_array))]
        next_mz_array_size = next_mz_array[mask2].size
        self.mz_array = np.append(self.mz_array, next_mz_array[mask2])

        n_i_a_m = next_intensity_array[mask2]
        if not (self.ion_mobility is None):
            n_im_a_m = next_ion_mobility_array[mask2]
        n_m_a_m = next_mz_array[mask2]
        for i in range(next_mz_array_size):
            self.scan_id.append([next_scan_id, ])
            self.intensity.append([n_i_a_m[i], ])
            if not (self.ion_mobility is None):
                self.ion_mobility.append([n_im_a_m[i], ])
            self.mass_array.append([n_m_a_m[i], ])

        self.selfsort()

    def selfsort(self):
        idx = np.argsort(self.mz_array)
        self.mz_array = self.mz_array[idx]
        self.scan_id = [self.scan_id[i] for i in idx]
        self.intensity = [self.intensity[i] for i in idx]
        if not (self.ion_mobility is None):
            self.ion_mobility = [self.ion_mobility[i] for i in idx]
        self.mass_array = [self.mass_array[i] for i in idx]

    def cutting_down(self, intensity_propotion):

        for idx, peak in enumerate(self.finished_hills):

            max_intensity_propotion = peak.max_intensity * intensity_propotion
            # FIXME try "and"

            if (
                peak.intensity[0] >= max_intensity_propotion and
                    peak.intensity[-1] >= max_intensity_propotion):

                del self.finished_hills[idx]


    def split_peaks2(self, hillValleyFactor):

        set_to_del = set()
        new_hills = []
        for hill_idx, hill in enumerate(self.finished_hills):

            if len(hill.mass) >= 6:

                mz_diff = np.array([z - hill.mz for z in hill.mass])
                std_5 = np.std(np.diff(mz_diff))
                smothed_intensity = list(np.abs(np.diff(mz_diff))/std_5)

                c_len = len(smothed_intensity) - 3
                idx = 3
                min_idx_list = []
                min_val = 1.0
                while idx <= c_len:
                    mult_val = smothed_intensity[idx]
                    if mult_val >= hillValleyFactor:
                        # if not len(min_idx_list) or idx >= min_idx_list[-1] + 3:
                        #     min_idx_list.append(idx)
                        #     min_val = mult_val
                        # elif mult_val < min_val:
                        #     min_idx_list[-1] = idx
                        #     min_val = mult_val
                        if (not len(min_idx_list) or idx >= min_idx_list[-1] + 3) and max(hill.intensity[0:idx-1]) >= 1.5 * max(hill.intensity[0], hill.intensity[idx-1]) and max(hill.intensity[idx:]) >= 1.5 * max(hill.intensity[idx], hill.intensity[-1]):
                            min_idx_list.append(idx)
                            min_val = mult_val
                        elif (mult_val < min_val) and max(hill.intensity[0:idx-1]) >= 1.5 * max(hill.intensity[0], hill.intensity[idx-1]) and max(hill.intensity[idx:]) >= 1.5 * max(hill.intensity[idx], hill.intensity[-1]):
                            min_idx_list[-1] = idx
                            min_val = mult_val
                    idx += 1

                if len(min_idx_list):
                    set_to_del.add(hill_idx)
                    prev_idx = 1
                    for min_idx in min_idx_list:
                        new_hills.append(ready_hill(
                                            intensity=hill.intensity[prev_idx-1:min_idx],
                                            scan_id=hill.scan_id[prev_idx-1:min_idx],
                                            mass=hill.mass[prev_idx-1:min_idx],
                                            ion_mobility=(
                                                hill.ion_mobility[prev_idx-1:min_idx] if not
                                                (hill.ion_mobility is None) else
                                                None)))
                        prev_idx = min_idx

                    new_hills.append(ready_hill(
                                        intensity=hill.intensity[min_idx-1:],
                                        scan_id=hill.scan_id[min_idx-1:],
                                        mass=hill.mass[min_idx-1:],
                                        ion_mobility=(
                                            hill.ion_mobility[min_idx-1:] if not
                                            (hill.ion_mobility is None) else
                                            None)))

        print(len(self.finished_hills))

        for idx in sorted(list(set_to_del))[::-1]:
            del self.finished_hills[idx]
        
        print(len(self.finished_hills))
        self.finished_hills.extend(new_hills)

        print(len(self.finished_hills))

    def calc_accurate_mz(self):
        for hill in self.finished_hills:
            hill.mz = sum(weight * value for weight, value in zip(hill.intensity, hill.mass)) / sum(hill.intensity)

    def split_peaks(self, hillValleyFactor, min_length_hill):
        set_to_del = set()
        new_hills = []
        for hill_idx, hill in enumerate(self.finished_hills):

            hill_length = len(hill.intensity)

            if hill_length >= min_length_hill * 2:

            # smothed_intensity = hill.intensity

                smothed_intensity = meanfilt(hill.intensity, 2)

                # smothed_intensity = medfilt(smothed_intensity, 3)


                # smothed_intensity = medfilt(hill.intensity, 3)
                # smothed_intensity = meanfilt(smothed_intensity, 3)

                c_len = hill_length - min_length_hill
                idx = int(min_length_hill)
                # min_idx = False
                min_idx_list = []
                min_val = 0
                l_idx = 0

                while idx <= c_len:

                    if len(min_idx_list) and idx >= min_idx_list[-1] + min_length_hill:
                        l_idx = min_idx_list[-1]

                    l_r = max(smothed_intensity[l_idx:idx]) / float(smothed_intensity[idx])
                    if l_r >= hillValleyFactor:
                        r_r = max(smothed_intensity[idx:]) / float(smothed_intensity[idx])
                        if r_r >= hillValleyFactor:
                #     print(l_r, r_r)
                    # if l_r >= hillValleyFactor and r_r >= hillValleyFactor:
                            mult_val = l_r * r_r
                            # if mult_val < min_val:
                            #     min_val = mult_val
                            if not len(min_idx_list) or idx >= min_idx_list[-1] + min_length_hill:
                                min_idx_list.append(idx)
                                min_val = mult_val
                            elif mult_val > min_val:
                                min_idx_list[-1] = idx
                                min_val = mult_val
                                # min_idx = idx
                    idx += 1
                if len(min_idx_list):
                    set_to_del.add(hill_idx)
                    prev_idx = 0
                    for min_idx in min_idx_list:
                        new_hills.append(ready_hill(
                                            intensity=hill.intensity[prev_idx:min_idx+1],
                                            scan_id=hill.scan_id[prev_idx:min_idx+1],
                                            mass=hill.mass[prev_idx:min_idx+1],
                                            ion_mobility=(
                                                hill.ion_mobility[prev_idx:min_idx+1] if not
                                                (hill.ion_mobility is None) else
                                                None)))
                        prev_idx = min_idx

                    new_hills.append(ready_hill(
                                        intensity=hill.intensity[min_idx:],
                                        scan_id=hill.scan_id[min_idx:],
                                        mass=hill.mass[min_idx:],
                                        ion_mobility=(
                                            hill.ion_mobility[min_idx:] if not
                                            (hill.ion_mobility is None) else
                                            None)))
        # print(len(new_hills))
        # print(len(set_to_del))


        for idx in sorted(list(set_to_del))[::-1]:
            del self.finished_hills[idx]
        
        self.finished_hills.extend(new_hills)

    # self.finished_hills = result


class feature:

    def __init__(self, finished_hills, each, each_id, negative_mode, isotopes_mass_error_map, mass_accuracy):

        self.charge = each[1][0][1]
        self.shift = each[3]
        # self.mz = finished_hills[each[0]].mz

        # a_cus = 0.0033946045716987906 / 1000
        # b_cus = -1.8123641799696435

        mass_for_average2 = [np.average(finished_hills[each[0]].mass, weights=finished_hills[each[0]].intensity)]
        intensity_for_average2 = [finished_hills[each[0]].max_intensity, ]

        # for i_numb, ech in enumerate(each[1]):
        #     mass_for_average2.append(np.average(finished_hills[ech[0]].mass, weights=finished_hills[ech[0]].intensity) - (i_numb+1)*1.00335/ech[1])
        #     intensity_for_average2.append(finished_hills[ech[0]].max_intensity)

        # mass_for_average2 = [zm * (1 - 1e-6 * (a_cus * zi + b_cus)) for zm, zi in zip(mass_for_average2, intensity_for_average2)]
        self.mz = np.average(mass_for_average2, weights=intensity_for_average2)
        
        
        mass_acc = mass_accuracy
        self.mz_tol = mass_acc*1e-6*self.mz

        # mass_for_average = finished_hills[each[0]].mass + list(itertools.chain.from_iterable([(z * (1 - 1e-6 * isotopes_mass_error_map[i_numb+1][0]) - (i_numb+1)*1.00335/ech[1]) for z in finished_hills[ech[0]].mass] for i_numb, ech in enumerate(each[1])))
        # # mass_for_average = finished_hills[each[0]].mass + list(itertools.chain.from_iterable([(z - (i_numb+1)*1.00335/ech[1]) for z in finished_hills[ech[0]].mass] for i_numb, ech in enumerate(each[1])))
        intensity_for_average = finished_hills[each[0]].intensity + list(itertools.chain.from_iterable(finished_hills[ech[0]].intensity for ech in each[1]))
        # # mass_for_average = [zm * (1 - 1e-6 * (a_cus * zi + b_cus)) for zm, zi in zip(mass_for_average, intensity_for_average)]
        # scans_for_average = finished_hills[each[0]].scan_id + list(itertools.chain.from_iterable(finished_hills[ech[0]].scan_id for ech in each[1]))
        # # print(mass_for_average, intensity_for_average)
        # self.mz = np.average(mass_for_average, weights=intensity_for_average)
        # # self.mz = np.median(mass_for_average)

        scans_for_average = finished_hills[each[0]].scan_id + list(itertools.chain.from_iterable(finished_hills[ech[0]].scan_id for ech in each[1]))


        # self.mz = np.median(finished_hills[each[0]].mass)
        # self.mz = np.mean(finished_hills[each[0]].mass)
        self.negative_mode = negative_mode

        if negative_mode == True:
            self.neutral_mass = self.mz * self.charge + \
                1.0072765 * self.charge - self.shift * 1.00335
        else:
            self.neutral_mass = self.mz * self.charge - \
                1.0072765 * self.charge - self.shift * 1.00335

        self.isotopes_numb = len(each[1])

        self.scan_numb = len(finished_hills[each[0]].scan_id)
        self.scans = finished_hills[each[0]].scan_id
        self.id_for_scan = finished_hills[each[0]].intensity.index(
            max(finished_hills[each[0]].intensity))
        self.intensity = finished_hills[each[0]].max_intensity

        # self.mz = self.mz * (1 - 1e-6 * (a_cus * max(intensity_for_average2) + b_cus))

        # self.id_for_scan = intensity_for_average.index(
        #     max(intensity_for_average))
        # self.intensity = max(intensity_for_average)
        self.idict = finished_hills[each[0]].idict
        self.sqrt_of_i_sum_squares = math.sqrt(
            sum(v**2 for v in self.idict.values()))
        self.scan_set = finished_hills[each[0]].scan_set
        if not (finished_hills[each[0]].ion_mobility is None):
            self.ion_mobility = finished_hills[each[0]].opt_ion_mobility
        else:
            self.ion_mobility = None

        # self.scan_id = scans_for_average[self.id_for_scan]
        # self.scan_id = finished_hills[each[0]].scan_id[self.id_for_scan]
        # self.RT = self.scan_numb
        self.scan_id = int(np.average(scans_for_average, weights=intensity_for_average))
        self.RT = int(np.average(scans_for_average, weights=intensity_for_average))
        # self.sulfur = (1 if each[2] else 0)
        self.sulfur = (each[1][1][4] if len(each[1]) > 1 else -1)
        self.cos_corr = each[4][0]
        self.cos_corr_2 = each[4][1]
        self.corr_fill_zero = each[4][2]
        self.diff_for_output = each[4][3]
        self.intensity_1 = each[4][4]
        self.scan_id_1 = each[4][5]
        self.mz_std_1 = np.std(each[4][6])
        self.intensity_2 = each[4][7]
        self.scan_id_2 = each[4][8]
        self.mz_std_2 = np.std(each[4][9])
        self.id = each_id
        self.ms2_scan = []
        
    def targeted(self, scan):
        self.ms2_scan.append(scan)
