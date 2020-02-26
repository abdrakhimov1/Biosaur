from copy import copy
import numpy as np
from scipy.signal import medfilt
import math


def meanfilt(data, window_width):
    cumsum_vec = np.cumsum(np.insert(data, 0, 0))
    ma_vec = (cumsum_vec[window_width:] -
              cumsum_vec[:-window_width]) / window_width
    return ma_vec


class ready_hill:

    def __init__(self, intensity, scan_id, mass, ion_mobility):

        self.mz = np.median(mass)
        self.mz_std = np.std(mass)
        self.intensity = intensity
        self.scan_id = scan_id
        self.scan_set = set(scan_id)
        self.mass = mass
        self.diff_for_output = 0
        tmp = max(range(len(self.intensity)), key=self.intensity.__getitem__)
        self.scan_of_max_intensity = self.scan_id[tmp]
        self.max_intensity = self.intensity[tmp]
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

        self.sqrt_of_i_sum_squares = math.sqrt(
            sum(v**2 for v in self.idict.values()))


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
        if self.intensity_array[index][-1] > self.intensity_array[index][-2]:
            self.intensity_max[index] = self.intensity_array[index][-1]
            self.ion_mobility_opt[index] = self.ion_mobility_array[index][-1]

    def push_me_to_the_peak(self, mz, intensity, ion_mobility, diff):
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
                        ion_mobility) <= 0.1 or abs(
                        self.ion_mobility_min[nearest_id] -
                        ion_mobility) <= 0.1:
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

        crosslink_counter = 0
        # crosslink_counter2 = 0
        self.finished_hills = sorted(
            self.finished_hills,
            key=lambda x: x.scan_id[0])

        allowed_ids = set()
        for i in self.intervals:
            allowed_ids.add(i - 1)
            allowed_ids.add(i)
            allowed_ids.add(i + 1)

        i = 0
        ini_len = len(self.finished_hills)

        while i < ini_len:

            hill = self.finished_hills[i]

            if hill.scan_id[-1] in allowed_ids:
                j = i + 1

                while j < ini_len:

                    hill2 = self.finished_hills[j]

                    if hill2.scan_id[0] in allowed_ids:

                        # if hill.scan_id[-1] == hill2.scan_id[0]:
                        if abs(hill.scan_id[-1] - hill2.scan_id[0]) <= 1:
                            # crosslink_counter2 += 1
                            if abs(hill.mz - hill2.mz) / \
                                    hill.mz <= mass_accuracy * 1e-6:

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
                                ini_len -= 1
                                crosslink_counter += 1
                                j -= 1
                        elif hill2.scan_id[0] > hill.scan_id[-1] + 1:
                            break

                    elif hill2.scan_id[0] > hill.scan_id[-1] + 1:
                        break

                    j += 1

            i += 1

        # print(crosslink_counter)
        # print(crosslink_counter2)

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

        mask_to_del = [True] * self.mz_array.size
        for i in range(self.mz_array.size)[::-1]:

            degree_actual = id_real - self.scan_id[i][-1]
            # or (degree_actual == 2 and len(self.scan_id[i]) == 1):
            if degree_actual > check_degree:

                # degree_actual = id_real - self.scan_id[i][-1]
                # if degree_actual > check_degree or (degree_actual == 2 and
                # len(self.scan_id[i]) <= 3):

                list_intensity = self.intensity.pop(i)
                if not (self.ion_mobility is None):
                    list_ion_mobility = self.ion_mobility.pop(i)
                else:
                    list_ion_mobility = None
                list_scan_id = self.scan_id.pop(i)
                list_mass = self.mass_array.pop(i)
                if len(list_scan_id) >= min_length:
                    tmp_ready_hill = ready_hill(intensity=list_intensity,
                                                scan_id=list_scan_id,
                                                mass=list_mass,
                                                ion_mobility=list_ion_mobility,
                                                )
                    self.finished_hills.append(tmp_ready_hill)

                mask_to_del[i] = False

                # if len(tmp_ready_hill.scan_id) >= min_length:
                #     self.finished_hills.append(tmp_ready_hill)

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

    def get_nearest_value(self, value, mask):
        return np.argmin(np.abs(self.mz_array[mask] - value))

    def newid(self, nearest, mask):
        return np.nonzero(mask)[0][nearest]

    def get_nearest_id(self, i, prev_nearest, diff, mz_array_l, ion_mobility):
        mass_diff = diff * 1e-6 * i
        best_diff = 2 * mass_diff
        best_id = False
        cur_md_abs = 0
        best_prev_nearest_id = False

        nearest_id = prev_nearest
        while nearest_id < mz_array_l:
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
            elif cur_md > mass_diff:
                break

            nearest_id += 1
        if not best_prev_nearest_id:
            best_prev_nearest_id = prev_nearest
        return best_id, best_diff / i, best_prev_nearest_id

    def get_arrays(self, tmp1):
        tmp1_nearest_id_arr = np.array([x[0] for x in tmp1])
        tmp1_idx_arr = np.array([x[1] for x in tmp1])
        tmp1_diff_arr = np.array([x[2] for x in tmp1])
        return tmp1_nearest_id_arr, tmp1_idx_arr, tmp1_diff_arr

    def push_me_to_the_peak(self, next_peak, diff, min_length):

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
                        else None))
            if best_id:
                tmp1.append([best_id, idx, md_res])

        mask = [True] * (len(self.mz_array))

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
                            nearest = self.get_nearest_value(element_mz, mask)
                            nearest_id = self.newid(nearest, mask)
                            tmp1_nearest_id_arr[idx] = nearest_id

                            tmp1_diff_arr[idx] = abs(
                                self.mz_array[nearest_id] - element_mz) / \
                                element_mz
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
            self.mz_array[i] = np.mean(self.mass_array[i][-3:])

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

    def split_peaks(self, hillValleyFactor):
        set_to_del = set()
        new_hills = []
        for hill_idx, hill in enumerate(self.finished_hills):
            smothed_intensity = meanfilt(hill.intensity, 3)
            smothed_intensity = medfilt(smothed_intensity, 3)

            c_len = len(smothed_intensity) - 3
            idx = 3
            min_idx = False
            min_val = 1.0
            while idx <= c_len:
                l_r = float(smothed_intensity[idx]) / \
                    max(smothed_intensity[:idx])
                r_r = float(smothed_intensity[idx]) / \
                    max(smothed_intensity[idx:])
            #     print(l_r, r_r)
                if l_r < hillValleyFactor and r_r < hillValleyFactor:
                    mult_val = l_r * r_r
                    if mult_val < min_val:
                        min_val = mult_val
                        min_idx = idx
                idx += 1
            if min_idx:
                set_to_del.add(hill_idx)
                new_hills.append(ready_hill(
                                    intensity=hill.intensity[:min_idx],
                                    scan_id=hill.scan_id[:min_idx],
                                    mass=hill.mass[:min_idx],
                                    ion_mobility=(
                                        hill.ion_mobility[:min_idx] if not
                                        (hill.ion_mobility is None) else
                                        None)))
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

    def __init__(self, finished_hills, each, each_id, negative_mode):

        self.charge = each[1][0][1]
        self.shift = each[3]
        # self.mz = finished_hills[each[0]].mz
        self.mz = np.median(finished_hills[each[0]].mass)
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
        self.idict = finished_hills[each[0]].idict
        self.sqrt_of_i_sum_squares = math.sqrt(
            sum(v**2 for v in self.idict.values()))
        self.scan_set = finished_hills[each[0]].scan_set
        if not (finished_hills[each[0]].ion_mobility is None):
            self.ion_mobility = finished_hills[each[0]].opt_ion_mobility
        else:
            self.ion_mobility = None

        self.scan_id = finished_hills[each[0]].scan_id[self.id_for_scan]
        # self.scan_id = finished_hills[each[0]]
        self.RT = self.scan_numb
        self.sulfur = (1 if each[2] else 0)
        self.cos_corr = each[4][0]
        self.cos_corr_2 = each[4][1]
        self.corr_fill_zero = each[4][2]
        self.diff_for_output = each[4][3]
        self.intensity_1 = each[4][4]
        self.scan_id_1 = each[4][5]
        self.mz_std_1 = each[4][6]
        self.intensity_2 = each[4][7]
        self.scan_id_2 = each[4][8]
        self.mz_std_2 = each[4][9]
        self.id = each_id
