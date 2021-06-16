import time
from os import path
from pyteomics import mzml
import pandas as pd
from pyteomics import pepxml
import numpy as np
from . import funcs
from . import classes
import logging
from . import utills
import os
logging.basicConfig(format=u'%(filename)s[LINE:%(lineno)d]#\
%(levelname)-8s [%(asctime)s]  %(message)s', level=logging.DEBUG)


def process_files(args):
    start_time = time.time()
    logging.info(u'Reading scans...')
    input_mzml_path = args['input_mzml_path'][0]
    number_of_processes = int(args['number_of_processes'])
    mass_accuracy = args['mass_accuracy']
    min_length = args['min_length']
    min_length_hill = args['min_length_hill']
    min_charge = args['min_charge']
    max_charge = args['max_charge']
    min_intensity = args['min_intensity']
    pep_xml_file_path = args['pep_xml_file_path']

    if args['negative_mode']:
        negative_mode = True
    else:
        negative_mode = False

    if args['output_file']:
        output_file = args['output_file']
    else:
        output_file = path.splitext(input_mzml_path)[0]\
             + path.extsep + 'features.tsv'
    hillValleyFactor = args['hill_valley_factor']

    correlation_map = args['correlation_map']

    # data_for_analyse = list(z for z in mzml.read(
    #     input_mzml_path) if z['ms level'] == 1)

    data_for_analyse = []
    for z in mzml.read(input_mzml_path):
        if z['ms level'] == 1:

            # if 5.3 <= float(z['scanList']['scan'][0]['scan start time']) <= 5.5:
            if 1:

                idx = z['intensity array'] >= min_intensity
                z['intensity array'] = z['intensity array'][idx]
                z['m/z array'] = z['m/z array'][idx]
                if 'mean inverse reduced ion mobility array' in z:
                    z['mean inverse reduced ion mobility array'] = z['mean inverse reduced ion mobility array'][idx]

                idx = np.argsort(z['m/z array'])
                z['m/z array'] = z['m/z array'][idx]
                z['intensity array'] = z['intensity array'][idx]
                if 'mean inverse reduced ion mobility array' in z:
                    z['mean inverse reduced ion mobility array'] = z['mean inverse reduced ion mobility array'][idx]


                data_for_analyse.append(z)
                # if len(data_for_analyse) > 100:
                #     break

    logging.info(u'Number of MS1 scans: ' + str(len(data_for_analyse)))
    tmp_str = 'maximum amount of'
    logging.info(
        u'Converting your data, using ' +
        str(number_of_processes if number_of_processes != 0 else tmp_str) +
        ' processes...')

    out_file = open(output_file, 'w')
    if correlation_map:
        out_file.write('massCalib\
\trtApex\
\tintensityApex\
\tcharge\
\tnIsotopes\
\tnScans\
\tsulfur\
\tcos_corr_1\
\tcos_corr_2\
\tdiff_for_output\
\tcorr_fill_zero\
\tintensity_1\
\tscan_id_1\
\tmz_std_1\
\tintensity_2\
\tscan_id_2\
\tmz_std_2\
\tmz\
\trtStart\
\trtEnd\
\tcorrMap\
\tid\
\tion_mobility\
\tFAIMS\
\ttargeted_mode\
\n')
    else:
        out_file.write('massCalib\
\trtApex\
\tintensityApex\
\tcharge\
\tnIsotopes\
\tnScans\
\tsulfur\
\tcos_corr_1\
\tcos_corr_2\
\tdiff_for_output\
\tcorr_fill_zero\
\tintensity_1\
\tscan_id_1\
\tmz_std_1\
\tintensity_2\
\tscan_id_2\
\tmz_std_2\
\tmz\
\trtStart\
\trtEnd\
\tid\
\tion_mobility\
\tFAIMS\
\ttargeted_mode\
\n')

    out_file.close()

    if 'FAIMS compensation voltage' in data_for_analyse[0]:
        logging.info(u'\nFAIMS was detected in data')
        faims_set = set()
        for z in data_for_analyse:
            if z['FAIMS compensation voltage'] not in faims_set:
                faims_set.add(z['FAIMS compensation voltage'])
            else:
                break
        logging.info(u'Detected FAIMS values: %s' % str(faims_set))
    else:
        faims_set = set([None, ])

    data_start_index = 0

    for faims_val in faims_set:

        if faims_val is None:
            data_for_analyse_tmp = data_for_analyse
            faims_val = 0
        else:
            data_for_analyse_tmp = []
            for z in data_for_analyse:
                if z['FAIMS compensation voltage'] == faims_val:
                    data_for_analyse_tmp.append(z)


        if 'mean inverse reduced ion mobility array' in data_for_analyse_tmp[0]:
            for idx, i in enumerate(data_for_analyse_tmp):
                # print(idx, len(i['m/z array']))
                peak_ion_mobility_object = False
                for mz, intensity, ion_mobility in zip(
                        i['m/z array'],
                        i['intensity array'],
                        i['mean inverse reduced ion mobility array']):
                    # if intensity >= 300:
                    if not peak_ion_mobility_object:
                        peak_ion_mobility_object = classes.peak_ion_mobility(
                            mz, intensity, ion_mobility)
                    else:
                        peak_ion_mobility_object.push_me_to_the_peak_ion_mob(
                            mz, intensity, ion_mobility, mass_accuracy)
                if peak_ion_mobility_object is False:
                    data_for_analyse_tmp[idx]['m/z array'] = np.array([])
                    data_for_analyse_tmp[idx]['intensity array'] = np.array([])
                    tmp_string = 'mean inverse reduced ion mobility array'
                    data_for_analyse_tmp[idx][tmp_string] = np.array([])
                else:
                    # peak_ion_mobility_object.mz_array = [np.mean(z) for z in peak_ion_mobility_object.mass_array]
                    peak_ion_mobility_object.ion_mobility_opt = [np.mean(z) for z in peak_ion_mobility_object.ion_mobility_array]
                    peak_ion_mobility_object.intensity_max = [np.sum(z) for z in peak_ion_mobility_object.intensity_array]
                    data_for_analyse_tmp[idx]['m/z array'] = np.array(
                        peak_ion_mobility_object.mz_array)
                    data_for_analyse_tmp[idx]['intensity array'] = np.array(
                        peak_ion_mobility_object.intensity_max)
                    tmp_string = 'mean inverse reduced ion mobility array'
                    data_for_analyse_tmp[idx][tmp_string] = np.array(
                        peak_ion_mobility_object.ion_mobility_opt)


        test_peak, test_RT_dict = funcs.boosting_firststep_with_processes(
            number_of_processes, data_for_analyse_tmp, mass_accuracy,
            min_length_hill, hillValleyFactor, data_start_index=data_start_index)
        

        data_start_index += len(data_for_analyse_tmp)

        # funcs.cutting_filter_on(test_peak)

        # print(
        #     "Timer: " +
        #     str(round((data_to_features_time - start_time) / 60, 1))
        #     + " minutes.")

        # test_peak.crosslink(mass_accuracy)
        # test_peak.cutting_down(0.5)
        logging.info(u'All data converted to hills...')
        logging.info('Processing hills...')
        logging.info(
            'Your hills proccesing with ' +
            str(number_of_processes if number_of_processes != 0 else tmp_str) +
            ' processes...')

        # test_peak.crosslink_simple(mass_accuracy)
        # print(
        #     "Timer: " +
        #     str(round((time.time() - start_time) / 60, 1)) + " minutes.")
        # test_peak.split_peaks(hillValleyFactor)
        # print(
        #     "Timer: " +
        #     str(round((time.time() - start_time) / 60, 1)) + " minutes.")
        # test_peak.split_peaks(hillValleyFactor)


        
        logging.info(
            str(len(test_peak.finished_hills)) +
            u' hills were detected...')

        # test_peak.split_peaks2(hillValleyFactor)
        logging.info(
            str(len(test_peak.finished_hills)) +
            u' hills were detected...')

        test_peak.sort_finished_hills()


        logging.info('Start recalc_fast_array_for_finished_hills...')

        mz_step = mass_accuracy * 1e-6 * 2500
        test_peak.recalc_fast_array_for_finished_hills(mz_step)

        # output = open('/home/mark/first_step.pkl', 'wb')
        # import pickle
        # pickle.dump(test_peak, output)
        # output.close()

        # pickle.dump(test_RT_dict, open('/home/mark/test_RT_dict.pkl', 'wb'))
        

        logging.info('Start boosting_secondstep_with_processes...')

        tmp, isotopes_mass_error_map = funcs.boosting_secondstep_with_processes(
            number_of_processes,
            test_peak,
            min_charge,
            max_charge,
            min_intensity,
            mass_accuracy,
            min_length)
        # tmp = funcs.iter_hills(test_peak, 1 , 5, 10, mass_accuracy)
        # output = open('second_step.pkl', 'wb')
        # pickle.dump(tmp, output)
        # output.close()

        # print(
        #     "Timer: " +
        #     str(round((iter_hills_time - start_time) / 60, 1)) + " minutes.")

        features = []
        
        

        for each_id, each in enumerate(tmp):
            features.append(
                classes.feature(
                    test_peak.finished_hills,
                    each,
                    each_id,
                    negative_mode, isotopes_mass_error_map, mass_accuracy))
        
        
        
        
        if pep_xml_file_path != '0':
#             targeted_pepxml_file = pepxml.read(pep_xml_file_path)
            targeted_dataframe = utills.prepare_dataframe(pep_xml_file_path)[0]
            targeted_mode_dict = dict()
            for index, i in targeted_dataframe.iterrows():
                targeted_mode_dict[i['spectrum']] = {
                    'RT' : i['RT exp'], 
                    'expect_score' : i['expect'],
                    'mz' :(i['calc_neutral_pep_mass'] + i['assumed_charge'] * 1.0072) / i['assumed_charge']}
                
            
            for idx, f in enumerate(features):
                keys_to_del = []
                for key, value in targeted_mode_dict.items():
                    if abs(f.mz - value['mz']) < f.mz_tol:
                        
                        if test_RT_dict[f.scans[0]] < value['RT']:

                            if value['RT'] < test_RT_dict[f.scans[-1]]:

                                f.targeted((key, value['expect_score']))
                                keys_to_del.append(key)
                  
                for each in keys_to_del:
                    del targeted_mode_dict[each]
                print(len(targeted_mode_dict))

            new_features = []

            for i in features:
                if not len(i.ms2_scan) == 0:
                    new_features.append(i)

            features = new_features

        if os.path.exists(output_file):
            dft = pd.read_table(output_file)
            max_id_val = dft['id'].max()
            if np.isnan(max_id_val):
                new_id = 1
            else:
                new_id = max_id_val + 1
        else:
            new_id = 1
        for ftr in features:
            ftr.id = new_id
            new_id += 1
        
        
        
        
        
        
        # print(
        #     "Timer: " +
        #     str(round((features_time - start_time) / 60, 1)) + " minutes.")

        out_file = open(output_file, 'a')

        if correlation_map:
            features = sorted(features, key=lambda x: x.scans[0])
            #FIXME исправить количество процессов для подсчета матрицы кореляции (сейчас 1 процесс)
            tmp_dict = funcs.boosting_correlation_matrix_with_processes(1, features)
            for x in features:
                out_file.write('\t'.join([str(z) for z in [
                    x.neutral_mass,
                    test_RT_dict[x.scan_id],
                    x.intensity,
                    x.charge,
                    x.isotopes_numb + 1,
                    x.scan_numb,
                    x.sulfur,
                    x.cos_corr,
                    x.cos_corr_2,
                    x.diff_for_output,
                    x.corr_fill_zero,
                    x.intensity_1,
                    x.scan_id_1,
                    x.mz_std_1,
                    x.intensity_2,
                    x.scan_id_2,
                    x.mz_std_2,
                    x.mz,
                    test_RT_dict[x.scans[0]],
                    test_RT_dict[x.scans[-1]],
                    tmp_dict[x.id],
                    x.id,
                    (
                        x.ion_mobility if not
                        (x.ion_mobility is None)
                        else 0),
                    faims_val,
                    x.ms2_scan]]) + '\n')
            out_file.close()
        else:
            for x in features:
                out_file.write('\t'.join([str(z) for z in [
                    x.neutral_mass,
                    test_RT_dict[x.scan_id],
                    x.intensity,
                    x.charge,
                    x.isotopes_numb + 1,
                    x.scan_numb,
                    x.sulfur,
                    x.cos_corr,
                    x.cos_corr_2,
                    x.diff_for_output,
                    x.corr_fill_zero,
                    x.intensity_1,
                    x.scan_id_1,
                    x.mz_std_1,
                    x.intensity_2,
                    x.scan_id_2,
                    x.mz_std_2,
                    x.mz,
                    test_RT_dict[x.scans[0]],
                    test_RT_dict[x.scans[-1]],
                    x.id,
                    (
                        x.ion_mobility if not
                        (x.ion_mobility is None)
                        else 0),
                    faims_val,
                    x.ms2_scan]]) + '\n')
            out_file.close()

        if pep_xml_file_path != '0':    
            bio = pd.read_table(output_file)
            ms_ms_df = bio.sort_values(['cos_corr_1', 'nIsotopes'], ascending=[False, False]).drop_duplicates(subset='targeted_mode', keep="last")
            ms_ms_df['ms_2'] = ms_ms_df['targeted_mode'].apply(lambda x: x.split(',')[0][3:-1:])
            ms_ms_df['ms_2_expect_val'] = ms_ms_df['targeted_mode'].apply(lambda x: x.split(',')[1][0:-2])
            bio = ms_ms_df.drop(['targeted_mode'], axis=1)
            bio.to_csv(output_file)
        
        total_time = time.time()
        print('=========================================================== \n')
        logging.info("Ready!")
        logging.info(str(len(features)) + u' features were detected.')
        logging.info(u'All your features were saved in ' + output_file)
        # print("Ready!")
        # print(
        #     "Total time: " +
        #     str(round((total_time - start_time) / 60, 1)) + " minutes.")
        logging.info(
            "Total time: " +
            str(round((total_time - start_time) / 60, 1)) +
            " minutes.")
