import time
from os import path
from pyteomics import mzml
from . import funcs
from . import classes
import logging
logging.basicConfig(format=u'%(filename)s[LINE:%(lineno)d]#\
%(levelname)-8s [%(asctime)s]  %(message)s', level=logging.DEBUG)


def process_files(args):
    start_time = time.time()
    logging.info(u'Reading scans...')
    input_mzml_path = args['input_mzml_path'][0]
    number_of_processes = int(args['number_of_processes'])
    mass_accuracy = args['mass_accuracy']
    min_length = args['min_length']
    min_charge = args['min_charge']
    max_charge = args['max_charge']
    min_intensity = args['min_intensity']

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

    data_for_analyse = list(z for z in mzml.read(
        input_mzml_path) if z['ms level'] == 1)
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
\n')

    out_file.close()

    if not args['faims']:
        faims_set = set([None, ])
        if 'FAIMS compensation voltage' in data_for_analyse[0]:
            logging.warning(u'\nWARNING: FAIMS detected in data,\
                 but option --faims was not enabled!\n')

    else:
        faims_set = set()
        for z in data_for_analyse:
            if z['FAIMS compensation voltage'] not in faims_set:
                faims_set.add(z['FAIMS compensation voltage'])
            else:
                break
        logging.info(u'Detected FAIMS values: ', faims_set)

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

        test_peak, test_RT_dict = funcs.boosting_firststep_with_processes(
            number_of_processes, data_for_analyse_tmp, mass_accuracy,
            min_length, data_start_index=data_start_index)

        data_start_index += len(data_for_analyse_tmp)

        # funcs.cutting_filter_on(test_peak)

        # print(
        #     "Timer: " +
        #     str(round((data_to_features_time - start_time) / 60, 1))
        #     + " minutes.")

        # test_peak.crosslink(mass_accuracy)
        # test_peak.cutting_down(0.5)
        logging.info(u'All data converted to hills...')
        logging.info(
            str(len(test_peak.finished_hills)) +
            u' hills were detected...')
        logging.info('Processing hills...')
        logging.info(
            'Your hills proccesing with ' +
            str(number_of_processes if number_of_processes != 0 else tmp_str) +
            ' processes...')
        test_peak.crosslink_simple(mass_accuracy)
        # print(
        #     "Timer: " +
        #     str(round((time.time() - start_time) / 60, 1)) + " minutes.")
        test_peak.split_peaks(hillValleyFactor)
        # print(
        #     "Timer: " +
        #     str(round((time.time() - start_time) / 60, 1)) + " minutes.")
        # test_peak.split_peaks(hillValleyFactor)

        test_peak.sort_finished_hills()

        # output = open('first_step.pkl', 'wb')
        # pickle.dump(test_peak, output)
        # output.close()

        tmp = funcs.boosting_secondstep_with_processes(
            number_of_processes,
            test_peak,
            min_charge,
            max_charge,
            min_intensity,
            mass_accuracy)
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
                    negative_mode))

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
                    faims_val]]) + '\n')
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
                    faims_val]]) + '\n')
            out_file.close()

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
