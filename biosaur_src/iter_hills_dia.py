import itertools
from multiprocessing import Queue, Process, cpu_count

from scipy.stats import binom
from .funcs import cos_correlation_fill_zeroes, cos_correlation, checking_cos_correlation_for_carbon, calibrate_mass
import numpy as np


def iter_hills_dia(
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

    # s_list = [1, 2, 3]
    # s_dict = dict()

    # for i in s_list:
    #     int_arr = binom.pmf(tmplist, i, 0.0425)
    #     s_dict[i] = int_arr

    tmplist_s = [1, 2, 3]

    s_list = [1, 2, 3]
    s_dict = dict()

    for i in s_list:
        int_arr = binom.pmf(tmplist_s, i, 0.0425)
        s_dict[i] = int_arr

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

            for charge in charges:

                candidates = []
                s_candidates = []

                neutral_mass = peak_1_mz * charge
                tmp_intensity = a[int(100 * (neutral_mass // 100))]
                s_tmp_intensity = s_list

                s_all_theoretical_int = [
                    peak.finished_hills[i].max_intensity *
                    s_tmp_intensity[z] /
                    s_tmp_intensity[0] for z in numbers[:2]]

                k = i
                ks = i
                for numb in numbers[1:]:

                    tmp_candidates = []
                    tmp_s_candidates = []


                    m_to_check = peak_1_mz + (1.00335 * numb / charge)
                    m_to_check_fast = int(m_to_check/0.02)

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


                    if numb == 2:

                        m_to_check = peak_1_mz + \
                            (1.9957958999999974 / 2 * numb / charge)
                        m_to_check_fast = int(m_to_check/0.02)

                        # for j in range(ks + 1, size, 1):
                        for j in peak.get_potential_isotope_id(m_to_check_fast, i):

                            if j not in ready_set and j != candidates[-1][0]:
                                peak_2_mz = peak.finished_hills[j].mz
                                diff = peak_2_mz - m_to_check

                                # if diff > mz_tol:
                                #     ks = j - 1
                                #     break

                                # if
                                # cos_correlation_fill_zeroes(peak.finished_hills[i],
                                # peak.finished_hills[j]) >= 0.7:
                                if abs(diff) <= mz_tol:
                                    s_cos_cor = cos_correlation_fill_zeroes(
                                            peak.finished_hills[i],
                                            peak.finished_hills[j])
                                    if s_cos_cor >= 0.6:

                                        tmp_s_candidates.append([j, charge, s_cos_cor, diff/m_to_check*1e6])

                                        # if len(s_candidates) < numb / 2:

                                        #     s_candidates.append(
                                        #         (j, charge))

                                        # else:

                                        #     intensity1 = (
                                        #         peak.finished_hills[
                                        #             s_candidates[-1][0]]
                                        #         .max_intensity)
                                        #     intensity2 = (
                                        #         peak.finished_hills[j]
                                        #         .max_intensity)
                                        #     s_th_i = s_all_theoretical_int[
                                        #         int(numb / 2)]
                                        #     s_candidates[-1] = (
                                        #         first_or_second(
                                        #             s_candidates[-1][0],
                                        #             j,
                                        #             s_candidates[-1][1],
                                        #             charge,
                                        #             intensity1,
                                        #             intensity2,
                                        #             s_th_i))

                                            # pass
                    if len(tmp_candidates):
                        # if len(tmp_candidates) > 1:
                        #     print(len(tmp_candidates))
                        candidates.append(tmp_candidates)

                        if len(tmp_s_candidates):
                            s_candidates = tmp_s_candidates


                    if len(candidates) < numb:
                        break

                if len(s_candidates):

                    tmp_s_candidates = []

                    for iter_s_candidates in s_candidates:


                        s_all_exp_intensity = [peak.finished_hills[i].max_intensity]
                        # for k in iter_s_candidates:
                        s_all_exp_intensity.append(
                            peak.finished_hills[iter_s_candidates[0]].max_intensity)

                        s_c_cor = cos_correlation(
                                s_all_theoretical_int,
                                s_all_exp_intensity)
                        if s_c_cor > 0.6:
                            tmp_s_candidates.append(iter_s_candidates)
                            tmp_s_candidates[-1].append(s_c_cor)

                    if len(tmp_s_candidates):
                        s_candidates = sorted(tmp_s_candidates, key=lambda x: -x[3])
                    else:
                        s_candidates = []

                for iter_candidates in itertools.product(*candidates):

                    all_theoretical_int = [
                        peak.finished_hills[i].max_intensity *
                        tmp_intensity[z] /
                        tmp_intensity[0] for z in numbers]

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

                    if len(iter_candidates):
                        j2 = iter_candidates[0][0]
                        scan_id_2 = peak.finished_hills[j2].scan_id
                        mz_std_2 = peak.finished_hills[j2].mz_std
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
                                peak.finished_hills[i].mz_std,
                                intensity_2,
                                scan_id_2,
                                mz_std_2],
                                [all_theoretical_int, all_exp_intensity]])
    print(len(ready))
    return ready


def worker_iter_hills_dia(
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

    result_q = iter_hills_dia(
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


def boosting_secondstep_dia_with_processes(
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

        ready = iter_hills_dia(
            peak, min_charge, max_charge, min_intensity, mass_accuracy, 0, len(
                peak.finished_hills), min_length)

    else:
        qout = Queue()
        procs = []

        peak_len = len(peak.finished_hills)
        step = int(peak_len / number_of_processes) + 1
        start_index = 0

        for i in range(number_of_processes):
            p = Process(
                target=worker_iter_hills_dia,
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
                isotopes_mass_error_map[icc + 1].append(cand[3])
    for ic in range(1, 10, 1):
        if len(isotopes_mass_error_map[ic]) >= 1000:

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
                isotopes_mass_error_map[ic] = isotopes_mass_error_map[ic - 1]

        else:
            isotopes_mass_error_map[ic] = isotopes_mass_error_map[ic - 1]
    # print(isotopes_mass_error_map)

    for pfidx, pep_feature in enumerate(ready):
        allowed_idx = 1
        for icc, cand in enumerate(pep_feature[1]):
            if abs(cand[3] - isotopes_mass_error_map[icc + 1][0]) / isotopes_mass_error_map[icc + 1][1] <= 5 or (
                    icc == 1 and cand[4] > 0):
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
    ready = sorted(ready, key=lambda x: -len(x[-2][0]) - x[-1])
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
    ready = sorted(ready, key=lambda x: -len(x[-2][0]) - x[-1])
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