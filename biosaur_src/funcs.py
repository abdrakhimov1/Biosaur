from . import classes
import numpy as np
from scipy.stats import binom
import math
from multiprocessing import Queue, Process, cpu_count
import logging
from collections import defaultdict
logging.basicConfig(format=u'%(filename)s[LINE:%(lineno)d]#\
%(levelname)-8s [%(asctime)s] %(message)s', level=logging.DEBUG)


def check_its_ready(id_real, peak, check_degree):

    mz_list = peak.mz_array
    scan_list = peak.scan_id

    for i in range(len(mz_list)):
        if id_real - scan_list[i][-1] > check_degree:

            tmp_ready_hill = classes.ready_hill(
                intensity=peak.intensity.pop(i),
                scan_id=peak.scan_id.pop(i),
                mass=peak.mass_array.pop(i),
                mz=peak.mz_array.pop(i))

            peak.finished_hills.append(tmp_ready_hill)


def data_to_features(input_file, max_diff, min_length, proccess_number):

    # data = mzml.read(input_file)

    # working_area = next(data)

    # working_area = input_file[0]

    # idx = working_area['intensity array'] >= 1
    # working_area['m/z array'] = working_area['m/z array'][idx]
    # working_area['intensity array'] = working_area['intensity array'][idx]

    # print(len(input_file[0]['m/z array']))
    if 'mean inverse reduced ion mobility array' in input_file[0]:
        for idx, i in enumerate(input_file):
            # print(idx, len(i['m/z array']))
            peak_ion_mobility_object = False
            for mz, intensity, ion_mobility in zip(
                    i['m/z array'],
                    i['intensity array'],
                    i['mean inverse reduced ion mobility array']):
                if intensity >= 300:
                    if not peak_ion_mobility_object:
                        peak_ion_mobility_object = classes.peak_ion_mobility(
                            mz, intensity, ion_mobility)
                    else:
                        peak_ion_mobility_object.push_me_to_the_peak(
                            mz, intensity, ion_mobility, max_diff)
            input_file[idx]['m/z array'] = np.array(
                peak_ion_mobility_object.mz_array)
            input_file[idx]['intensity array'] = np.array(
                peak_ion_mobility_object.intensity_max)
            tmp_string = 'mean inverse reduced ion mobility array'
            input_file[idx][tmp_string] = np.array(
                peak_ion_mobility_object.ion_mobility_opt)
        # if idx > 10:
        #     break
    # print(len(peak_ion_mobility_object.mz_array))
    # print(len(i['m/z array']))
    # print(len(input_file[0]['m/z array']))

    # print(input_file[0]['m/z array'])
    # print(input_file[0]['intensity array'])
    # print('HERE')

    RT_dict = dict()

    k = 0
    # print(len(input_file))

    for i in input_file:
        # print(i)
        if k == 0:
            peak1 = classes.peak(
                i['m/z array'],
                i['intensity array'],
                i['index'],
                i['index'],
                i.get(
                    'mean inverse reduced ion mobility array',
                    None))
            RT_dict[i['index']] = float(
                i['scanList']['scan'][0]['scan start time'])

        if k > 0:
            next_peak_i = classes.next_peak(
                i['m/z array'],
                i['intensity array'],
                i['index'],
                i.get(
                    'mean inverse reduced ion mobility array',
                    None))
            peak1.push_me_to_the_peak(next_peak_i, max_diff, min_length)
            RT_dict[i['index']] = float(
                i['scanList']['scan'][0]['scan start time'])
        # if k > 10:
            # break
            # pass
        k += 1
    # print(peak1.mz_array)
    peak1.push_left(min_length=min_length)
    logging.info(
        u'Data converted to features with process /' +
        str(proccess_number + 1) + '/ --->')
    # print(peak1.mz_array)
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

    bottom = math.sqrt(sum([numb * numb for numb in theoretical_list])) * \
        math.sqrt(sum([numb * numb for numb in experimental_list]))

    for i1, i2 in zip(theoretical_list, experimental_list):
        top += i1 * i2

    return top / bottom, sum(theoretical_list) / theor_total_sum


def cos_correlation_fill_zeroes(hill_1, hill_2):

    inter_set = hill_1.scan_set.intersection(hill_2.scan_set)
    if len(inter_set):

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
            if averagineCorrelation >= thresh:
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

    return best_cor, best_pos, best_shift


def iter_hills(
        peak,
        min_charge,
        max_charge,
        min_intensity,
        mass_acc,
        start_index,
        end_index,
        proccess_number):

    ready = []
    averagine_mass = 111.1254
    averagine_C = 4.9384
    tmplist = list(range(10))
    prec_masses = []
    # prec_isotopes = []
    # prec_minisotopes = []
    # isotopes_int = []

    a = dict()

    # s_list = [1, 2, 3]
    # s_dict = dict()

    # for i in s_list:
    #     int_arr = binom.pmf(tmplist, i, 0.0425)
    #     s_dict[i] = int_arr

    s_list = [0.9575, 0.0425, 0.0425**2, 0.0425**3]

    for i in range(100, 20000, 100):
        int_arr = binom.pmf(
            tmplist,
            float(i) /
            averagine_mass *
            averagine_C,
            0.0107)
        prec_masses.append(i)
        int_arr_norm = int_arr / int_arr.max()
        # prec_is = np.where(int_arr_norm >= 0.01)[0]
#         isotopes_int.append(int_arr_norm[prec_is])
#         prec_minisotopes.append(prec_is.min())
#         prec_isotopes.append(prec_is - prec_minisotopes[-1])
        a[i] = int_arr_norm

    end_index = min(end_index, len(peak.finished_hills))
    size = end_index

    ready_set = set()

    for i in range(start_index, end_index, 1):

        if i not in ready_set:

            candidates = []
            s_candidates = []
            numbers = []

            for k in range(10):
                numbers.append(k)

            charges = list(range(min_charge, max_charge + 1, 1)[::-1])

            peak_1_mz = peak.finished_hills[i].mz
            left_border_i = peak.finished_hills[i].scan_id[0]
            right_border_i = peak.finished_hills[i].scan_id[-1]
            # middle_index_i = peak.finished_hills[i].scan_of_max_intensity
            # i_len = peak.finished_hills[i].scan_len
            mz_tol = mass_acc * 1e-6 * peak.finished_hills[i].mz
            # mz_tol = 5 * peak.finished_hills[i].mz_std

            for charge in charges:

                neutral_mass = peak_1_mz * charge
                tmp_intensity = a[int(100 * (neutral_mass // 100))]
                s_tmp_intensity = s_list

                all_theoretical_int = [
                    peak.finished_hills[i].max_intensity *
                    tmp_intensity[z] /
                    tmp_intensity[0] for z in numbers]
                s_all_theoretical_int = [
                    peak.finished_hills[i].max_intensity *
                    s_tmp_intensity[z] /
                    s_tmp_intensity[0] for z in numbers[:2]]

                k = i
                ks = i
                for numb in numbers[1:]:

                    m_to_check = peak_1_mz + (1.00335 * numb / charge)

                    for j in range(k + 1, size, 1):

                        if j not in ready_set:

                            peak_2_mz = peak.finished_hills[j].mz
                            left_border_j = (
                                            peak.finished_hills[j].scan_id[0]
                                            - 1)
                            right_border_j = (
                                            peak.finished_hills[j].scan_id[-1]
                                            + 1)
                            # j_len = peak.finished_hills[j].scan_len

                            diff = peak_2_mz - m_to_check
                            if numb == 1:
                                diff_for_output = diff / peak_2_mz

                            if diff > mz_tol:
                                # print(k)
                                k = j - 1
                                # print(k)
                                break
                            if abs(diff) <= mz_tol:

                                cos_cor_test = cos_correlation_fill_zeroes(
                                                    peak.finished_hills[i],
                                                    peak.finished_hills[j])

                                if cos_cor_test >= 0.7:

                                    if len(candidates) < numb:

                                        candidates.append((j, charge))

                                    else:

                                        j_prev = candidates[-1][0]
                                        intensity2 = (
                                                    peak.finished_hills[j]
                                                    .max_intensity)

                                        if abs(
                                            peak.finished_hills[j].mz -
                                            (peak.finished_hills[i]
                                                .mz)) < abs(
                                            (peak.finished_hills[j_prev]
                                                .mz) -
                                                peak.finished_hills[i].mz):
                                            candidates[-1] = (j, charge)

                    # lc = len(candidates)
                    # if lc < numb:
                    #     if lc and candidates[-1][1] != 0:
                    #         candidates.append((0, 0))
                    #     else:
                    #         if lc:
                    #             candidates = candidates[:-1]
                    #         break

                    # if numb % 2 == 0:
                    if numb == 2:

                        m_to_check = peak_1_mz + \
                            (1.9957958999999974 / 2 * numb / charge)

                        for j in range(ks + 1, size, 1):

                            if j not in ready_set and j != candidates[-1][0]:
                                peak_2_mz = peak.finished_hills[j].mz
                                left_border_j = (
                                    peak.finished_hills[j]
                                    .scan_id[0]) - 1
                                # j_len = peak.finished_hills[j].scan_len
                                right_border_j = (
                                    peak.finished_hills[j]
                                    .scan_id[-1]) + 1
                                diff = peak_2_mz - m_to_check

                                if diff > mz_tol:
                                    ks = j - 1
                                    break

                                # if
                                # cos_correlation_fill_zeroes(peak.finished_hills[i],
                                # peak.finished_hills[j]) >= 0.7:
                                if abs(diff) <= mz_tol:
                                    if 1 or (
                                            left_border_i -
                                            1 <= left_border_j and
                                            right_border_i +
                                            1 >= right_border_j):
                                        # if abs(middle_index_i -
                                        # middle_index_j) <= max(1, 0.5
                                        # *max(j_len, i_len)):

                                        if cos_correlation_fill_zeroes(
                                                peak.finished_hills[i],
                                                peak.finished_hills[j]) >= 0.7:

                                            if len(s_candidates) < numb / 2:

                                                s_candidates.append(
                                                    (j, charge))

                                            else:

                                                intensity1 = (
                                                    peak.finished_hills[
                                                        s_candidates[-1][0]]
                                                    .max_intensity)
                                                intensity2 = (
                                                    peak.finished_hills[j]
                                                    .max_intensity)
                                                s_th_i = s_all_theoretical_int[
                                                    int(numb / 2)]
                                                s_candidates[-1] = (
                                                    first_or_second(
                                                        s_candidates[-1][0],
                                                        j,
                                                        s_candidates[-1][1],
                                                        charge,
                                                        intensity1,
                                                        intensity2,
                                                        s_th_i))

                                                pass

                    if len(candidates) < numb:
                        break

                if len(candidates) > 0:  # FIXME

                    break

            if candidates:

                all_exp_intensity = [peak.finished_hills[i].max_intensity]
                s_all_exp_intensity = [peak.finished_hills[i].max_intensity]

                for j in candidates:
                    if j[1] != 0:
                        all_exp_intensity.append(
                            peak.finished_hills[j[0]].max_intensity)
                    else:
                        all_exp_intensity.append(0)

                for k in s_candidates:
                    s_all_exp_intensity.append(
                        peak.finished_hills[k[0]].max_intensity)

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

                    candidates = candidates[:number_of_passed_isotopes]

                    # добавить s_candidates
                    j2 = candidates[0][0]
                    scan_id_2 = peak.finished_hills[j2].scan_id
                    mz_std_2 = peak.finished_hills[j2].mz_std
                    intensity_2 = peak.finished_hills[j2].intensity
                    ready.append([
                        i,
                        candidates,
                        [],
                        shift,
                        [cos_corr,
                            cos_corr_for_output,
                            cos_cor_test,
                            diff_for_output,
                            peak.finished_hills[i].intensity,
                            peak.finished_hills[i].scan_id,
                            peak.finished_hills[i].mz_std,
                            intensity_2,
                            scan_id_2,
                            mz_std_2]])

                    ready_set.add(i)
                    for i in candidates:
                        if i[1] != 0:
                            ready_set.add(i[0])

                    if cos_correlation(
                            s_all_theoretical_int,
                            s_all_exp_intensity) > 0.6:

                        ready[-1][2] = s_candidates

                        for i in s_candidates:
                            ready_set.add(i[0])
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
        min_length, proccess_number):

    result_peak, result_RT_dict = data_to_features(
        data_for_analyse[start_index:end_index],
        mass_accuracy,
        min_length,
        proccess_number
        )

    if result_peak:
        qout.put((result_peak, result_RT_dict))

    qout.put(None)


def boosting_firststep_with_processes(
        number_of_processes,
        data_for_analyse,
        mass_accuracy,
        min_length,
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
            data_for_analyse, mass_accuracy, min_length)

    else:
        qout = Queue()

#         qin = list(islice(it, 500000))
#         if not len(qin):
#             break
# #           print 'Loaded 500000 items. Ending cycle.'
        procs = []

        data_for_analyse_len = len(data_for_analyse)
        step = int(data_for_analyse_len / number_of_processes) + 1
        start_index = 0

        for i in range(number_of_processes):
            p = Process(
                target=worker_data_to_features,
                args=(
                    data_for_analyse,
                    qout,
                    start_index,
                    step + start_index,
                    mass_accuracy,
                    min_length, i))
            # print(start_index)
            p.start()
            procs.append(p)
            start_index += step

        result_peak = False
        result_RT_dict = False

        for _ in range(number_of_processes):
            for item in iter(qout.get, None):
                # print(len(item[0].finished_hills))
                if not result_peak:
                    # print(item[0].mz_array)
                    result_peak, result_RT_dict = item[0], item[1]

                else:

                    # print(item[0].mz_array)
                    result_peak.concat_peak_with(item[0])
                    result_RT_dict.update(item[1])
            # print(len(result_peak.finished_hills))
        for p in procs:
            p.join()

    return result_peak, result_RT_dict


def worker_iter_hills(
        peak,
        qout,
        start_index,
        end_index,
        min_charge,
        max_charge,
        min_intensity,
        mass_accuracy,
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
        proccess_number)

    if result_q:
        qout.put(result_q)

    qout.put(None)


def boosting_secondstep_with_processes(
        number_of_processes,
        peak,
        min_charge,
        max_charge,
        min_intensity,
        mass_accuracy):

    if number_of_processes == 0:

        try:
            number_of_processes = cpu_count()

        except NotImplementedError:
            number_of_processes = 1

    if number_of_processes == 1:

        result_q = iter_hills(
            peak, min_charge, max_charge, min_intensity, mass_accuracy, 0, len(
                peak.finished_hills))

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
                    i))
            # print(start_index)
            p.start()
            procs.append(p)
            start_index += step

        result_q = False

        for _ in range(number_of_processes):
            for item in iter(qout.get, None):
                # print(len(item[0].finished_hills))
                if not result_q:
                    # print(item[0].mz_array)
                    result_q = item

                else:
                    # print(item[0].mz_array)
                    result_q = result_q + item
            # print(len(result_peak.finished_hills))
        for p in procs:
            p.join()

    return result_q

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
