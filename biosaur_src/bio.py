from . import funcs
from . import classes
import time


def process_files(args):
    start_time = time.time()
    print('Starting program...')
    input_mzml_path = args['input_mzml_path'][0]
    number_of_processes = int(args['number_of_processes'])
    mass_accuracy = args['mass_accuracy']
    min_length = args['min_length']
    min_charge = args['min_charge']
    max_charge = args['max_charge']
    min_intensity = args['min_intensity']
    output_file = args['output_file']
    hillValleyFactor = args['hill_valley_factor']

    test_peak, test_RT_dict = funcs.boosting_firststep_with_processes(
        number_of_processes, input_mzml_path, mass_accuracy, min_length)

    # funcs.cutting_filter_on(test_peak)

    data_to_features_time = time.time()

    print(
        "Timer: " +
        str(round((data_to_features_time - start_time) / 60, 1)) + " minutes.")

    # test_peak.crosslink(mass_accuracy)
    # test_peak.cutting_down(0.5)
    print(len(test_peak.finished_hills))
    test_peak.crosslink_simple(mass_accuracy)
    print(
        "Timer: " +
        str(round((time.time() - start_time) / 60, 1)) + " minutes.")
    test_peak.split_peaks(hillValleyFactor)
    print(
        "Timer: " +
        str(round((time.time() - start_time) / 60, 1)) + " minutes.")
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
    iter_hills_time = time.time()

    # output = open('second_step.pkl', 'wb')
    # pickle.dump(tmp, output)
    # output.close()

    print(
        "Timer: " +
        str(round((iter_hills_time - start_time) / 60, 1)) + " minutes.")

    features = []

    for each in tmp:
        features.append(classes.feature(test_peak.finished_hills, each))

    features_time = time.time()

    print(
        "Timer: " +
        str(round((features_time - start_time) / 60, 1)) + " minutes.")

    out_file = open(output_file, 'w')
    out_file.write(
        'massCalib\
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
        \tion_mobility\
        \n')

    # output = open('step3.pkl', 'wb')
    # pickle.dump(features, output)
    # output.close()

    for x in features:
        out_file.write('\t'.join([str(z) for z in [x.neutral_mass,
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
                                                   (
                                                       x.ion_mobility if not
                                                       (x.ion_mobility is None)
                                                       else 0)]]) + '\n')
    out_file.close()

    total_time = time.time()
    print("Ready!")
    print(
        "Total time: " +
        str(round((total_time - start_time) / 60, 1)) + " minutes.")
