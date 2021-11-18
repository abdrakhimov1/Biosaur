from . import classes
import numpy as np
from scipy.stats import binom
from scipy.stats import scoreatpercentile
from scipy.optimize import curve_fit
import operator
import math
from multiprocessing import Queue, Process, cpu_count
import logging
import itertools
from copy import deepcopy
from collections import defaultdict
logging.basicConfig(format=u'%(filename)s[LINE:%(lineno)d]#\
%(levelname)-8s [%(asctime)s] %(message)s', level=logging.DEBUG)


# def check_its_ready(id_real, peak, check_degree):

#     mz_list = peak.mz_array
#     scan_list = peak.scan_id

#     for i in range(len(mz_list)):
#         if id_real - scan_list[i][-1] > check_degree:

#             tmp_ready_hill = classes.ready_hill(
#                 intensity=peak.intensity.pop(i),
#                 scan_id=peak.scan_id.pop(i),
#                 mass=peak.mass_array.pop(i),
#                 mz=peak.mz_array.pop(i))

#             peak.finished_hills.append(tmp_ready_hill)


def data_to_features(input_file, max_diff, min_length_hill, proccess_number, start_index, end_index):

    # data = mzml.read(input_file)

    # working_area = next(data)

    # working_area = input_file[0]

    # idx = working_area['intensity array'] >= 1
    # working_area['m/z array'] = working_area['m/z array'][idx]
    # working_area['intensity array'] = working_area['intensity array'][idx]

    # print(len(input_file[0]['m/z array']))
        # if idx > 10:
        #     break

    RT_dict = dict()

    mz_step = max_diff * 1e-6 * 2500

    k = 0
    # print(len(input_file))

    for i in input_file:
        idx = (i['m/z array'] >= start_index) & (i['m/z array'] < end_index)
        # (dists >= r) & (dists <= r+dr)
        new_mz_array = i['m/z array'][idx]
        new_intensity_array = i['intensity array'][idx]
        if 'mean inverse reduced ion mobility array' in i:
            new_ion_mobility_array = i['mean inverse reduced ion mobility array'][idx]
        else:
            new_ion_mobility_array = None
        if k == 0:
            peak1 = classes.peak(
                new_mz_array,
                new_intensity_array,
                i['index'],
                i['index'],
                new_ion_mobility_array,
                )
            RT_dict[i['index']] = float(
                i['scanList']['scan'][0]['scan start time'])

        if k > 0:
            next_peak_i = classes.next_peak(
                new_mz_array,
                new_intensity_array,
                i['index'],
                new_ion_mobility_array,
                )
            peak1.push_me_to_the_peak(next_peak_i, max_diff, min_length_hill, mz_step)
            RT_dict[i['index']] = float(
                i['scanList']['scan'][0]['scan start time'])
        # if k > 10:
            # break
            # pass
        k += 1
    peak1.push_left(min_length=min_length_hill)
    # peak1.medar = np.array(peak1.medar)
    # peak1.medar = np.cumprod(peak1.medar)

    # for idx in range(len(peak1.finished_hills)):
    #     tmp_mass = [mv * peak1.medar[sv-1] for mv, sv in zip(peak1.finished_hills[idx].mass, peak1.finished_hills[idx].scan_id)]
    #     peak1.finished_hills[idx].mass = tmp_mass
    #     peak1.finished_hills[idx].mz = np.median(tmp_mass)

    logging.info(
        u'Data converted to features with process /' +
        str(proccess_number + 1) + '/ --->')
    return peak1, RT_dict


def first_or_second(id1, id2, charge1, charge2, first, second, theoretiacal):

    if abs(theoretiacal - first) <= abs(theoretiacal - second):
        return id1, charge1
    else:
        return id2, charge2


def cos_correlation(theoretical_list, experimental_list):

    suit_len = min(len(theoretical_list), len(experimental_list))

    theoretical_list = theoretical_list[:suit_len]
    experimental_list = experimental_list[:suit_len]

    top = 0

    bottom = math.sqrt(sum([numb * numb for numb in theoretical_list])) * \
        math.sqrt(sum([numb * numb for numb in experimental_list]))

    for i1, i2 in zip(theoretical_list, experimental_list):
        top += i1 * i2

    return top / bottom


def cos_correlation_new(theoretical_list, experimental_list, shf):

    theor_total_sum = sum(theoretical_list)
    theoretical_list = theoretical_list[shf:]
    suit_len = min(len(theoretical_list), len(experimental_list))
    theoretical_list = theoretical_list[:suit_len]
    experimental_list = experimental_list[:suit_len]

    top = 0

    for i1, i2 in zip(theoretical_list, experimental_list):
        top += i1 * i2

    if not top:
        return 0, 0
    else:

        bottom = math.sqrt(sum([numb * numb for numb in theoretical_list])) * \
            math.sqrt(sum([numb * numb for numb in experimental_list]))

        averagineExplained = sum(theoretical_list) / theor_total_sum

        return top / bottom, averagineExplained


def cos_correlation_fill_zeroes(hill_1, hill_2):

    inter_set = hill_1.scan_set.intersection(hill_2.scan_set)
    if len(inter_set) >= 2:

        top = 0
        for i in inter_set:
            h1_val = hill_1.idict.get(i, 0)
            h2_val = hill_2.idict.get(i, 0)
            top += h1_val * h2_val
        
        
        bottom = hill_1.sqrt_of_i_sum_squares * hill_2.sqrt_of_i_sum_squares
        # bottom = math.sqrt(sum(v * v for key, v in hill_1.idict.items() if key in inter_set)) * math.sqrt(sum(v * v for key, v in hill_2.idict.items() if key in inter_set))
        return top / bottom

    else:
        return 0


def checking_cos_correlation_for_carbon_noshift(
        theoretical_list, experimental_list, thresh):

    # prev_corr = 0
    # size = 1
    best_value = 0
    best_shift = 0
    best_pos = 1
    best_cor = 0

    # for shf in range(4):
    for shf in range(1):
        # shf = 0
        pos = len(experimental_list)

        while pos != 1:

            averagineCorrelation, averagineExplained = cos_correlation_new(
                theoretical_list, experimental_list[:pos], shf)

            # if averagineExplained < 0.5:
            #     break
            if averagineExplained >= 0.5 and averagineCorrelation >= thresh:
                tmp_val = averagineCorrelation * averagineExplained
                if tmp_val > best_value:
                    best_value = tmp_val
                    best_cor = averagineCorrelation
                    best_shift = shf
                    best_pos = pos

                break

            pos -= 1

            # if correlation >= thresh:

            # if correlation >= prev_corr:

            #     prev_corr = correlation
            #     size = len(experimental_list)

            # else:
            # size = len(experimental_list)
            # return correlation, pos#, shift

        # experimental_list = experimental_list[:-1]

        if best_value:
            break

    return best_cor, best_pos, best_shift

def checking_cos_correlation_for_carbon(
        theoretical_list, experimental_list, thresh):

    # prev_corr = 0
    # size = 1
    best_value = 0
    best_shift = 0
    best_pos = 1
    best_cor = 0

    for shf in range(4):
        # shf = 0
        pos = len(experimental_list)

        while pos != 1:

            averagineCorrelation, averagineExplained = cos_correlation_new(
                theoretical_list, experimental_list[:pos], shf)

            # if averagineExplained < 0.5:
            #     break
            if averagineExplained >= 0.5 and averagineCorrelation >= thresh:
                tmp_val = averagineCorrelation# * averagineExplained
                if tmp_val > best_value:
                    best_value = tmp_val
                    best_cor = averagineCorrelation
                    best_shift = shf
                    best_pos = pos

                break

            pos -= 1

            # if correlation >= thresh:

            # if correlation >= prev_corr:

            #     prev_corr = correlation
            #     size = len(experimental_list)

            # else:
            # size = len(experimental_list)
            # return correlation, pos#, shift

        # experimental_list = experimental_list[:-1]

        if best_value:
            break

    return best_cor, best_pos, best_shift


def iter_hills(
        peak,
        min_charge,
        max_charge,
        min_intensity,
        mass_acc,
        start_index,
        end_index,
        min_length,
        proccess_number=1):

    ready = []
    averagine_mass = 111.1254
    averagine_C = 4.9384
    tmplist = list(range(10))
    prec_masses = []
    # prec_isotopes = []
    # prec_minisotopes = []
    # isotopes_int = []

    a = dict()

    mz_step = mass_acc * 1e-6 * 2500

    # s_list = [1, 2, 3]
    # s_dict = dict()

    # for i in s_list:
    #     int_arr = binom.pmf(tmplist, i, 0.0425)
    #     s_dict[i] = int_arr

    # tmplist_s = [1, 2, 3]

    # s_list = [1, 2, 3]
    # s_dict = dict()

    # for i in s_list:
    #     int_arr = binom.pmf(tmplist_s, i, 0.0425)
    #     s_dict[i] = int_arr

    # s_list = [s_dict[i] for i in tmplist_s]

    # s_list = [0.9575, 0.0425, 0.0425**2, 0.0425**3]

    for i in range(100, 20000, 100):
        int_arr = binom.pmf(
            tmplist,
            float(i) /
            averagine_mass *
            averagine_C,
            0.0107)
        prec_masses.append(i)
        # int_arr_norm = int_arr / int_arr.max()
        int_arr_norm = int_arr / int_arr.sum()
        # prec_is = np.where(int_arr_norm >= 0.01)[0]
#         isotopes_int.append(int_arr_norm[prec_is])
#         prec_minisotopes.append(prec_is.min())
#         prec_isotopes.append(prec_is - prec_minisotopes[-1])
        a[i] = int_arr_norm

    end_index = min(end_index, len(peak.finished_hills))
    size = end_index

    ready_set = set()
    
    charges = list(range(min_charge, max_charge + 1, 1)[::-1])

    numbers = []

    for k in range(10):
        numbers.append(k)

    for i in range(start_index, end_index, 1):

        if peak.finished_hills[i].scan_len >= min_length:

            peak_1_mz = peak.finished_hills[i].mz
            left_border_i = peak.finished_hills[i].scan_id[0]
            right_border_i = peak.finished_hills[i].scan_id[-1]
            mz_tol = mass_acc * 1e-6 * peak.finished_hills[i].mz

            # s_tmp_intensity = s_list

            # s_all_theoretical_int = [
            #     peak.finished_hills[i].max_intensity *
            #     s_tmp_intensity[z] /
            #     s_tmp_intensity[0] for z in numbers[:2]]

            for charge in charges:
                
                candidates = []
                s_candidates = []

                k = i
                ks = i
                for numb in numbers[1:]:

                    tmp_candidates = []
                    tmp_s_candidates = []

                    m_to_check = peak_1_mz + (1.00335 * numb / charge)
                    m_to_check_fast = int(m_to_check/mz_step)

                    for j in peak.get_potential_isotope_id(m_to_check_fast, i):

                        peak_2_mz = peak.finished_hills[j].mz

                        diff = peak_2_mz - m_to_check

                        if abs(diff) <= mz_tol and (peak.finished_hills[i].opt_ion_mobility is None or abs(peak.finished_hills[i].opt_ion_mobility-peak.finished_hills[j].opt_ion_mobility) <= 0.01):

                            cos_cor_test = cos_correlation_fill_zeroes(
                                                peak.finished_hills[i],
                                                peak.finished_hills[j])

                            if cos_cor_test >= 0.6:

                                tmp_candidates.append((j, charge, cos_cor_test, diff/m_to_check*1e6, 0))

                                if numb == 1:
                                    diff_for_output = diff / peak_2_mz


                    # if numb == 2:

                    #     for n_sulf in range(1, 4, 1):

                    #         m_to_check2 = peak_1_mz + (1.00335 * (numb - 2) / charge) + (1.9957958999999974 / charge)

                    #         sulf_int = s_dict[n_sulf][0]
                    #         # sulf_int = 0.0425 * tmp_intensity[numb-2]

                    #         m_to_check2 = (m_to_check2 * sulf_int + m_to_check * tmp_intensity[numb]) / (sulf_int+tmp_intensity[numb])
                    #         m_to_check2_fast = int(m_to_check2/0.02)
                    #         # print(m_to_check, m_to_check2)

                    #         for j in peak.get_potential_isotope_id(m_to_check2_fast, i):
                    #             if j not in ready_set and (peak.finished_hills[i].opt_ion_mobility is None or abs(peak.finished_hills[i].opt_ion_mobility-peak.finished_hills[j].opt_ion_mobility) <= 0.01):
                    #                 peak_2_mz = peak.finished_hills[j].mz
                    #                 diff = peak_2_mz - m_to_check2
                    #                 if abs(diff) <= mz_tol:
                    #                     cos_cor_test = cos_correlation_fill_zeroes(
                    #                                         peak.finished_hills[i],
                    #                                         peak.finished_hills[j])
                    #                     if cos_cor_test >= 0.6:
                    #                         tmp_candidates.append((j, charge, cos_cor_test, diff/m_to_check2*1e6, n_sulf))


                    # if numb == 2:

                    #     m_to_check = peak_1_mz + \
                    #         (1.9957958999999974 / 2 * numb / charge)
                    #     m_to_check_fast = int(m_to_check/0.02)

                    #     # for j in range(ks + 1, size, 1):
                    #     for j in peak.get_potential_isotope_id(m_to_check_fast, i):

                    #         if j not in ready_set and j != candidates[-1][0]:
                    #             peak_2_mz = peak.finished_hills[j].mz
                    #             diff = peak_2_mz - m_to_check

                    #             # if diff > mz_tol:
                    #             #     ks = j - 1
                    #             #     break

                    #             # if
                    #             # cos_correlation_fill_zeroes(peak.finished_hills[i],
                    #             # peak.finished_hills[j]) >= 0.7:
                    #             if abs(diff) <= mz_tol:
                    #                 s_cos_cor = cos_correlation_fill_zeroes(
                    #                         peak.finished_hills[i],
                    #                         peak.finished_hills[j])
                    #                 if s_cos_cor >= 0.6:

                    #                     tmp_s_candidates.append([j, charge, s_cos_cor, diff/m_to_check*1e6])

                    #                     # if len(s_candidates) < numb / 2:

                    #                     #     s_candidates.append(
                    #                     #         (j, charge))

                    #                     # else:

                    #                     #     intensity1 = (
                    #                     #         peak.finished_hills[
                    #                     #             s_candidates[-1][0]]
                    #                     #         .max_intensity)
                    #                     #     intensity2 = (
                    #                     #         peak.finished_hills[j]
                    #                     #         .max_intensity)
                    #                     #     s_th_i = s_all_theoretical_int[
                    #                     #         int(numb / 2)]
                    #                     #     s_candidates[-1] = (
                    #                     #         first_or_second(
                    #                     #             s_candidates[-1][0],
                    #                     #             j,
                    #                     #             s_candidates[-1][1],
                    #                     #             charge,
                    #                     #             intensity1,
                    #                     #             intensity2,
                    #                     #             s_th_i))

                    #                         # pass
                    if len(tmp_candidates):
                        # if len(tmp_candidates) > 1:
                        #     print(len(tmp_candidates))
                        candidates.append(tmp_candidates)

                        # if len(tmp_s_candidates):
                        #     s_candidates = tmp_s_candidates


                    if len(candidates) < numb:
                        break

                # if len(candidates) > 0:  # FIXME

                #     break

                if candidates:

                    neutral_mass = peak_1_mz * charge
                    tmp_intensity = a[int(100 * (neutral_mass // 100))]

                    all_theoretical_int = [
                        peak.finished_hills[i].max_intensity *
                        tmp_intensity[z] /
                        tmp_intensity[0] for z in numbers]

                    # s_candidates = []


                    # if len(s_candidates):

                    #     tmp_s_candidates = []

                    #     for iter_s_candidates in s_candidates:


                    #         s_all_exp_intensity = [peak.finished_hills[i].max_intensity]
                    #         # for k in iter_s_candidates:
                    #         s_all_exp_intensity.append(
                    #             peak.finished_hills[iter_s_candidates[0]].max_intensity)

                    #         s_c_cor = cos_correlation(
                    #                 s_all_theoretical_int,
                    #                 s_all_exp_intensity)
                    #         if s_c_cor > 0.6:
                    #             tmp_s_candidates.append(iter_s_candidates)
                    #             tmp_s_candidates[-1].append(s_c_cor)

                    #     if len(tmp_s_candidates):
                    #         s_candidates = sorted(tmp_s_candidates, key=lambda x: -x[3])
                    #     else:
                    #         s_candidates = []

                    for iter_candidates in itertools.product(*candidates):
                        # if len(iter_candidates) > 1:
                        #     basic_sulfur = iter_candidates[1][4]
                        # else:
                        #     basic_sulfur = 0
                        # # print(basic_sulfur)
                        # iter_candidates_new = []
                        # for z_idx, z in enumerate(iter_candidates):
                        #     if z_idx > 0:
                        #         if z[4] == basic_sulfur:
                        #             iter_candidates_new.append(z)
                        #         else:
                        #             break
                        #     else:
                        #         iter_candidates_new.append(z)
                        # iter_candidates = iter_candidates_new

                        all_exp_intensity = [peak.finished_hills[i].max_intensity]

                        for j in iter_candidates:
                            if j[1] != 0:
                                all_exp_intensity.append(
                                    peak.finished_hills[j[0]].max_intensity)
                            else:
                                all_exp_intensity.append(0)

                        (
                            cos_corr,
                            number_of_passed_isotopes,
                            shift) = checking_cos_correlation_for_carbon(
                            all_theoretical_int, all_exp_intensity, 0.6)

                        cos_corr_for_output = cos_correlation(
                                                all_theoretical_int[0:1],
                                                all_exp_intensity[0:1])

                        if cos_corr:  # прикрутить изменение параметра 0.6

                            # print(shift)

                            iter_candidates = iter_candidates[:number_of_passed_isotopes]

                            # добавить s_candidates
                            j2 = iter_candidates[0][0]
                            scan_id_2 = peak.finished_hills[j2].scan_id
                            mass_2 = peak.finished_hills[j2].mass
                            intensity_2 = peak.finished_hills[j2].intensity
                            ready.append([
                                i,
                                iter_candidates,
                                s_candidates,
                                shift,
                                [cos_corr,
                                    cos_corr_for_output,
                                    cos_cor_test,
                                    diff_for_output,
                                    peak.finished_hills[i].intensity,
                                    peak.finished_hills[i].scan_id,
                                    peak.finished_hills[i].mass,
                                    intensity_2,
                                    scan_id_2,
                                    mass_2],
                                    [all_theoretical_int, all_exp_intensity]])
                            # ready_set.add(i)
                            # for ic in candidates:
                            #     if ic[1] != 0:
                            #         ready_set.add(ic[0])


                                    # for ic in s_candidates:
                                    #     ready_set.add(ic[0])


    # ready = sorted(ready, key=lambda x: -len(x[1]))
    # ready_final = []
    # ready_set = set()


    # for pep_feature in ready:
    #     if pep_feature[0] not in ready_set:
    #         ready_final.append(pep_feature)
    #         ready_set.add(pep_feature[0])


    logging.info(
        u'All hills were iterated correctly with this process /' +
        str(proccess_number + 1) + '/ -->')
    return ready


def worker_data_to_features(
        data_for_analyse,
        qout,
        start_index,
        end_index,
        mass_accuracy,
        min_length_hill, hillValleyFactor, proccess_number):

    start_index = start_index * (1 - 1e-6 * 2 * mass_accuracy)
    end_index = end_index * (1 + 1e-6 * 2 * end_index)

    result_peak, result_RT_dict = data_to_features(
        data_for_analyse,
        mass_accuracy,
        min_length_hill,
        proccess_number,
        start_index,
        end_index
        )

    result_peak.split_peaks(hillValleyFactor, min_length_hill)

    set_to_del = set()
    for hill_idx, hill in enumerate(result_peak.finished_hills):
        if len(hill.mass) >= 40:
            if max(hill.intensity) < 2 * max(hill.intensity[0], hill.intensity[-1]):
                set_to_del.add(hill_idx)

    for idx in sorted(list(set_to_del))[::-1]:
        del result_peak.finished_hills[idx]

    result_peak.calc_accurate_mz()

    if result_peak:
        qout.put((result_peak, result_RT_dict))

    qout.put(None)


def boosting_firststep_with_processes(
        number_of_processes,
        data_for_analyse,
        mass_accuracy,
        min_length_hill,
        hillValleyFactor,
        data_start_index=0):

    
    for idx, v in enumerate(data_for_analyse):
        v['index'] = idx + 1 + data_start_index

    if number_of_processes == 0:

        try:
            number_of_processes = cpu_count()

        except NotImplementedError:
            number_of_processes = 1

    if number_of_processes == 1:

        result_peak, result_RT_dict = data_to_features(
            data_for_analyse, mass_accuracy, min_length_hill, 1, 0, 2500)


        result_peak.split_peaks(hillValleyFactor, min_length_hill)

        set_to_del = set()
        for hill_idx, hill in enumerate(result_peak.finished_hills):
            if len(hill.mass) >= 40:
                if max(hill.intensity) < 2 * max(hill.intensity[0], hill.intensity[-1]):
                    set_to_del.add(hill_idx)

        for idx in sorted(list(set_to_del))[::-1]:
            del result_peak.finished_hills[idx]

        result_peak.calc_accurate_mz()

    else:
        qout = Queue()

#         qin = list(islice(it, 500000))
#         if not len(qin):
#             break
# #           print 'Loaded 500000 items. Ending cycle.'
        procs = []

        data_for_analyse_len = len(data_for_analyse)
        # step = int(data_for_analyse_len / number_of_processes) + 1
        step = int(2500 / number_of_processes / 3) + 1
        # start_index = 0
        start_mz = 100

        for i in range(number_of_processes * 3):
            p = Process(
                target=worker_data_to_features,
                args=(
                    data_for_analyse,
                    qout,
                    start_mz,
                    step + start_mz,
                    mass_accuracy,
                    min_length_hill, hillValleyFactor, i))
            # print(start_index)
            p.start()
            procs.append(p)
            start_mz += step

        result_peak = False
        result_RT_dict = False#dict()

        # all_peaks = []
        for _ in range(number_of_processes * 3):
            for item in iter(qout.get, None):
                # all_peaks.append(item[0])
                # result_RT_dict.update(item[1])
                # print(len(item[0].finished_hills))
                if not result_peak:
                    # print(item[0].mz_array)
                    result_peak, result_RT_dict = item[0], item[1]

                else:

                    # print(item[0].mz_array)
                    result_peak.concat_peak_with(item[0])
                    result_RT_dict.update(item[1])
        # result_peak = concat_peaks(all_peaks)

            # print(len(result_peak.finished_hills))
        for p in procs:
            p.join()

    return result_peak, result_RT_dict

def concat_peaks(all_peaks):
    all_peaks = sorted(all_peaks, key=lambda x: x.intervals[0])
    result_peak = all_peaks[0]
    for peak in all_peaks[1:]:
        result_peak.concat_peak_new(peak)
    return result_peak


def worker_iter_hills(
        peak,
        qout,
        start_index,
        end_index,
        min_charge,
        max_charge,
        min_intensity,
        mass_accuracy,
        min_length,
        proccess_number
        ):

    result_q = iter_hills(
        peak,
        min_charge,
        max_charge,
        min_intensity,
        mass_accuracy,
        start_index,
        end_index,
        min_length,
        proccess_number)

    if result_q:
        qout.put(result_q)

    qout.put(None)


def noisygaus(x, a, x0, sigma, b):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b

def calibrate_mass(bwidth, mass_left, mass_right, true_md):

    bbins = np.arange(-mass_left, mass_right, bwidth)
    H1, b1 = np.histogram(true_md, bins=bbins)
    b1 = b1 + bwidth
    b1 = b1[:-1]

    popt, pcov = curve_fit(noisygaus, b1, H1, p0=[1, np.median(true_md), 1, 1])
    mass_shift, mass_sigma = popt[1], abs(popt[2])
    return mass_shift, mass_sigma, pcov[0][0]



def boosting_secondstep_with_processes(
        number_of_processes,
        peak,
        min_charge,
        max_charge,
        min_intensity,
        mass_accuracy,
        min_length):

    if number_of_processes == 0:

        try:
            number_of_processes = cpu_count()

        except NotImplementedError:
            number_of_processes = 1

    if number_of_processes == 1:

        ready = iter_hills(
            peak, min_charge, max_charge, min_intensity, mass_accuracy, 0, len(
                peak.finished_hills), min_length)

    else:
        qout = Queue()

#         qin = list(islice(it, 500000))
#         if not len(qin):
#             break
# #           print 'Loaded 500000 items. Ending cycle.'
        procs = []

        peak_len = len(peak.finished_hills)
        step = int(peak_len / number_of_processes) + 1
        start_index = 0

        for i in range(number_of_processes):
            p = Process(
                target=worker_iter_hills,
                args=(
                    peak,
                    qout,
                    start_index,
                    step + start_index,
                    min_charge,
                    max_charge,
                    min_intensity,
                    mass_accuracy,
                    min_length,
                    i))
            # print(start_index)
            p.start()
            procs.append(p)
            start_index += step

        ready = False

        for _ in range(number_of_processes):
            for item in iter(qout.get, None):
                if not ready:
                    ready = item
                else:
                    ready = ready + item
        for p in procs:
            p.join()

    ready = sorted(ready, key=lambda x: -len(x[1]))
    ready_final = []
    ready_set = set()

    # import pickle
    # pickle.dump(ready, open('ready.pickle', 'wb'))

    isotopes_mass_error_map = {}
    for ic in range(1, 10, 1):
        isotopes_mass_error_map[ic] = []
    for pep_feature in ready:
        for icc, cand in enumerate(pep_feature[1]):
            if icc != 1 or cand[4] == 0:
                isotopes_mass_error_map[icc+1].append(cand[3])
    for ic in range(1, 10, 1):
        if ic == 1 and len(isotopes_mass_error_map[ic]) >= 10:

            try:

                true_md = np.array(isotopes_mass_error_map[ic])

                mass_left = -min(isotopes_mass_error_map[ic])
                mass_right = max(isotopes_mass_error_map[ic])

                try:
                    mass_shift, mass_sigma, covvalue = calibrate_mass(0.01, mass_left, mass_right, true_md)
                except:
                    try:
                        mass_shift, mass_sigma, covvalue = calibrate_mass(0.05, mass_left, mass_right, true_md)
                    except:
                        mass_shift, mass_sigma, covvalue = calibrate_mass(0.25, mass_left, mass_right, true_md)
                if np.isinf(covvalue):
                    mass_shift, mass_sigma, covvalue = calibrate_mass(0.05, mass_left, mass_right, true_md)
                
                isotopes_mass_error_map[ic] = [mass_shift, mass_sigma]

            except:
                isotopes_mass_error_map[ic] = isotopes_mass_error_map[ic-1]

        else:
            if ic-1 in isotopes_mass_error_map:
                isotopes_mass_error_map[ic] = deepcopy(isotopes_mass_error_map[ic-1])
                isotopes_mass_error_map[ic][0] = isotopes_mass_error_map[ic][0] - 0.45
            else:
                isotopes_mass_error_map[ic] = [0, 10]
    print(isotopes_mass_error_map)

    for pfidx, pep_feature in enumerate(ready):
        allowed_idx = 1
        for icc, cand in enumerate(pep_feature[1]):
            if abs(cand[3] - isotopes_mass_error_map[icc+1][0])/isotopes_mass_error_map[icc+1][1] <= 5 or (icc == 1 and cand[4] > 0):
                allowed_idx += 1
            else:
                break


        all_theoretical_int, all_exp_intensity = pep_feature[5]
        all_theoretical_int = all_theoretical_int[:allowed_idx]
        all_exp_intensity = all_exp_intensity[:allowed_idx]

        ready[pfidx][1] = ready[pfidx][1][:allowed_idx]
        ready[pfidx][5] = [all_theoretical_int, all_exp_intensity]
    
        ready[pfidx].append(min(checking_cos_correlation_for_carbon(
                        all_theoretical_int, all_exp_intensity, 0.6)[0], 0.99999999))
        
    
    # ready = sorted(ready, key=lambda x: -x[-1])
    ready = sorted(ready, key=lambda x: -len(x[-2][0])-x[-1])
    # ready = sorted(ready, key=lambda x: -len(x[-2][0]))

    # for pep_feature in ready:
    #     if pep_feature[0] not in ready_set:
    #         if not any(cand[0] in ready_set for cand in pep_feature[1]):
    #             ready_final.append(pep_feature)
    #             ready_set.add(pep_feature[0])
    #             for cand in pep_feature[1]:
    #                 ready_set.add(cand[0])
    #             # for s_cand in pep_feature[2]:
    #             #     if s_cand[0] not in ready_set:
    #             #         ready_set.add(s_cand[0])
    #             #         break

    #         else:
    #             tmp = []
    #             for cand in pep_feature[1]:
    #                 if cand[0] not in ready_set:
    #                     tmp.append(cand)
    #                 else:
    #                     break
    #             if len(tmp):
    #                 pep_feature[1] = tmp
    #                 all_theoretical_int, all_exp_intensity = pep_feature[5]
    #                 all_theoretical_int = all_theoretical_int[:len(tmp)]
    #                 all_exp_intensity = all_exp_intensity[:len(tmp)]
    #                 (cos_corr,
    #                         number_of_passed_isotopes,
    #                         shift) = checking_cos_correlation_for_carbon(
    #                         all_theoretical_int, all_exp_intensity, 0.6)

    #                 if cos_corr:
    #                     ready_final.append(pep_feature)
    #                     ready_set.add(pep_feature[0])
    #                     for cand in pep_feature[1]:
    #                         ready_set.add(cand[0])
    #                     # for s_cand in pep_feature[2]:
    #                     #     if s_cand[0] not in ready_set:
    #                     #         ready_set.add(s_cand[0])
    #                     #         break
    max_l = len(ready)
    cur_l = 0

    ready_final = []
    ready_set = set()
    ready = sorted(ready, key=lambda x: -len(x[-2][0])-x[-1])
    cur_isotopes = len(ready[0][-2][0])

    cnt_mark = 0

    while cur_l < max_l:
        cnt_mark += 1
    #     if cnt_mark > 1000:
    #         break
        pep_feature = ready[cur_l]
        # print(cur_l, max_l, cur_isotopes, len(ready_final), -len(pep_feature[-2][0])-pep_feature[-1])
        n_iso = len(pep_feature[-2][0])
        if n_iso < cur_isotopes:
            ready = sorted(ready, key=lambda x: -len(x[-2][0]))
            cur_isotopes = n_iso
            cur_l = 0
            
        if pep_feature[0] not in ready_set:
            if not any(cand[0] in ready_set for cand in pep_feature[1]):
                ready_final.append(pep_feature)
                ready_set.add(pep_feature[0])
                for cand in pep_feature[1]:
                    ready_set.add(cand[0])
                for s_cand in pep_feature[2]:
                    ready_set.add(s_cand[0])
                del ready[cur_l]
                max_l -= 1
                cur_l -= 1

            else:
                tmp = []
                
    #             cur_isotopes = len(pep_feature[1])
                
                for cand in pep_feature[1]:
                    if cand[0] not in ready_set:
                        tmp.append(cand)
                    else:
                        break
                if len(tmp):
                    pep_feature[1] = tmp
                    all_theoretical_int, all_exp_intensity = pep_feature[5]
                    all_theoretical_int = all_theoretical_int[:len(tmp)]
                    all_exp_intensity = all_exp_intensity[:len(tmp)]
                    (cos_corr,
                            number_of_passed_isotopes,
                            shift) = checking_cos_correlation_for_carbon(
                            all_theoretical_int, all_exp_intensity, 0.6)
                    if cos_corr:
                        ready[cur_l] = [pep_feature[0],
                                        pep_feature[1],
                                        pep_feature[2],
                                        pep_feature[3],
                                        pep_feature[4],
                                        [all_theoretical_int, all_exp_intensity],
                                        cos_corr]

                    
                    else:
                        del ready[cur_l]
                        max_l -= 1
                        cur_l -= 1
                    
                    
                else:
                    del ready[cur_l]
                    max_l -= 1
                    cur_l -= 1
        else:
            del ready[cur_l]
            max_l -= 1
            cur_l -= 1


    #                 ready = ready[:cur_l] + sorted(ready[cur_l:], key=lambda x: -len(x[-2][0])-x[-1])

    #                 cur_l -= 1
        cur_l += 1
        

    return ready_final, isotopes_mass_error_map

#FIXME исправить функцию для подсчета по списку необходимых индексов 
def func_for_correlation_matrix(set_of_features):

    logging.info(u'Counting features correlation...')

    out_put_dict = defaultdict(list)

    set_length = len(set_of_features)
    each_id = 0
    # for each_id, each_feature in enumerate(set_of_features[:-1]):
    while each_id < set_length - 1:
        each_feature = set_of_features[each_id]
        if each_id % 50 == 0:
            logging.info(
                u'Calculated ' +
                str(each_id + 1) +
                '/' + str(set_length) +
                ' features.')

        other_id = each_id + 1
        while other_id < set_length:
            other_feature = set_of_features[other_id]

            if other_feature.scans[0] > each_feature.scans[-1]:
                break

            tmp_corr = cos_correlation_fill_zeroes(
                each_feature,
                other_feature)
            if tmp_corr > 0.5:
                out_put_dict[each_feature.id] += [{other_feature.id: tmp_corr}]
                out_put_dict[other_feature.id] += [{each_feature.id: tmp_corr}]
            other_id += 1

        each_id += 1

    return out_put_dict


def func_for_correlation_matrix_2(set_of_features, idx_list):

    logging.info(u'Counting features correlation...')

    out_put_dict = defaultdict(list)

    # for each_id, each_feature in enumerate(set_of_features[:-1]):
    for i in idx_list:
        each_feature = set_of_features[i]
        if i % 50 == 0:
            logging.info(
                u'Calculated ' +
                str(i + 1) +
                '/' + str(len(idx_list)) +
                ' features.')

        other_id = i + 1
        while other_id < idx_list[-1] + 1:
            other_feature = set_of_features[other_id]

            if other_feature.scans[0] > each_feature.scans[-1]:
                break

            tmp_corr = cos_correlation_fill_zeroes(
                each_feature,
                other_feature)
            if tmp_corr > 0.5:
                out_put_dict[each_feature.id] += [{other_feature.id: tmp_corr}]
                out_put_dict[other_feature.id] += [{each_feature.id: tmp_corr}]
            other_id += 1

    return out_put_dict


def worker_func_for_correlation_matrix(
        set_of_features,
        qout,
        idx_list
        ):

    result_q = func_for_correlation_matrix_2(
        set_of_features,
        idx_list,
            )

    if result_q:
        qout.put(result_q)

    qout.put(None)


def boosting_correlation_matrix_with_processes(
        number_of_processes, set_of_features):

    if number_of_processes == 0:

        try:
            number_of_processes = cpu_count()

        except NotImplementedError:
            number_of_processes = 1

    if number_of_processes == 1:

        result_q = func_for_correlation_matrix(
            set_of_features
        )

    else:
        qout = Queue()

#         qin = list(islice(it, 500000))
#         if not len(qin):
#             break
# #           print 'Loaded 500000 items. Ending cycle.'
        procs = []
#FIXME пофиксить(проверить) индексы в idx_list 
        set_len = len(set_of_features)
        step = int(set_len / 2 / number_of_processes) + 1
        start_index = 0

        for i in range(number_of_processes):

            idx_list = [x for x in range(i * step, (i + 1) * step)] + [x for x in range(set_len - (i+1) * step, set_len - i * step)]

            p = Process(
                target=worker_func_for_correlation_matrix,
                args=(
                    set_of_features,
                    qout,
                    idx_list))
            # print(start_index)
            p.start()
            procs.append(p)

        result_q = False

        for _ in range(number_of_processes):
            for item in iter(qout.get, None):
                # print(len(item[0].finished_hills))
                if not result_q:
                    # print(item[0].mz_array)
                    result_q = item

                else:
                    # print(item[0].mz_array)
                    for key in item:
                        result_q.update(item[key])
            # print(len(result_peak.finished_hills))
        for p in procs:
            p.join()

    return result_q
