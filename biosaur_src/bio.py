import funcs
from copy import copy
from pyteomics import mzml
import pandas as pd
from scipy import spatial
from scipy.spatial import KDTree
import numpy as np
from scipy.stats import binom
import classes
import time
import pickle



start_time = time.time()
print('Starting program...')
number_of_processes = 0
mass_accuracy = 8
min_length = 3
min_charge = 1
max_charge = 6
min_intensity = 1
hillValleyFactor = 0.8
# input_mzml_path = 'QEHF1_10924_JB.mzML'
#input_mzml_path = 'QEHF1_09766_JB.mzML'
input_mzml_path = 'QEHF1_10930_JB.mzML'
#input_mzml_path = 'plasma.mzML'
#test_peak, test_RT_dict = funcs.boosting_with_processes('plasma.mzML', mass_accuracy, 3)
test_peak, test_RT_dict = funcs.boosting_firststep_with_processes(number_of_processes, input_mzml_path, mass_accuracy, min_length)

print(np.histogram([z.scan_len for z in test_peak.finished_hills], bins = [2,3,4,5,10,20, 50, 100, 200, 500, 1000, 2000]))



data_to_features_time = time.time()

print("Timer: " + str(round((data_to_features_time - start_time) / 60, 1)) + " minutes.")

# test_peak.crosslink(mass_accuracy)
print(len(test_peak.finished_hills))
test_peak.crosslink_simple(mass_accuracy)
print(len(test_peak.finished_hills))
test_peak.split_peaks(hillValleyFactor)
print(len(test_peak.finished_hills))
test_peak.split_peaks(hillValleyFactor)
print(len(test_peak.finished_hills))

print(np.histogram([z.scan_len for z in test_peak.finished_hills], bins = [2,3,4,5,10,20, 50, 100, 200, 500, 1000, 2000]))

# # test_peak.cutting_down(0.5)
test_peak.sort_finished_hills()



# output = open('first_step.pkl', 'wb')
# pickle.dump(test_peak, output)
# output.close()

tmp = funcs.boosting_secondstep_with_processes(number_of_processes, test_peak, min_charge, max_charge, min_intensity, mass_accuracy)
#tmp = funcs.iter_hills(test_peak, 1 , 5, 10, mass_accuracy)
iter_hills_time = time.time()

# output = open('second_step.pkl', 'wb')
# pickle.dump(tmp, output)
# output.close()

print("Timer: " + str(round((iter_hills_time - start_time)/60, 1)) + " minutes.")

features = []

for each in tmp:
    features.append(classes.feature(test_peak.finished_hills, each))

features_time = time.time()

print("Timer: " + str(round((features_time - start_time)/60, 1)) + " minutes.")

out_file = open('test_custom.features.tsv', 'w')
out_file.write('massCalib\trtApex\tintensityApex\tcharge\tnIsotopes\tnScans\tsulfur\tshift\n')

# output = open('step3.pkl', 'wb')
# pickle.dump(features, output)
# output.close()

for x in features:
    out_file.write('\t'.join([str(z) for z in [x.neutral_mass, test_RT_dict[x.scan_id], x.intensity, x.charge, x.isotopes_numb + 1, x.scan_numb, x.sulfur, x.shift]]) + '\n')
out_file.close()

total_time = time.time()
print("Ready!")
print("Total time: " + str(round((total_time - start_time)/60, 1)) + " minutes.")
