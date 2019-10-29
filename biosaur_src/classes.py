from .funcs import *
from copy import copy
from pyteomics import mzml
import pandas as pd
from scipy import spatial
from scipy.spatial import KDTree
import numpy as np
from scipy.stats import binom

class ready_hill:
    
    def __init__(self, intensity, scan_id, mass):

        self.mz = np.mean(mass)
        
        self.intensity = intensity
        self.scan_id = scan_id
        self.mass = mass
        tmp = max(range(len(self.intensity)), key=self.intensity.__getitem__)
        self.scan_of_max_intensity = self.scan_id[tmp]
        self.max_intensity = self.intensity[tmp]
        self.scan_len = len(self.scan_id)


    def set_new_hill(self, start_index, end_index):
        
        new_hill = ready_hill(self.intensity[start_index:end_index], self.scan_id[start_index:end_index], self.mass[start_index:end_index])

        return new_hill

class next_peak:
    
    def __init__(self, next_mz_array, next_intensity_array, next_scan_id):
        
        self.next_mz_array = next_mz_array
        self.next_intensity_array = next_intensity_array
        self.next_scan_id = next_scan_id


        
class peak:
    def __init__(self, mz_array, intensity, scan_id):
        
        self.mz_array = copy(mz_array)
        
        self.scan_id = []

        for _ in range(len(mz_array)):
            self.scan_id.append([scan_id, ])

        
        
        self.intensity = []
        for i in intensity:
            self.intensity.append([i, ])
        
        self.mass_array = []
        for i in mz_array:
            self.mass_array.append([i, ])
            
        self.finished_hills = []
        self.crosslinked_hills = []
    
    def concat_peak_with(self, second_peak):
        
        self.mz_array = self.mz_array + second_peak.mz_array
        self.intensity = self.intensity + second_peak.intensity
        self.mass_array = self.mass_array + second_peak.mass_array
        self.finished_hills = self.finished_hills + second_peak.finished_hills
        self.crosslinked_hills = self.crosslinked_hills + second_peak.crosslinked_hills


    

    def crosslink(self, mass_accuracy):
        
        crosslink_counter = 0
        crosslink_counter2 = 0
        self.finished_hills = sorted(self.finished_hills, key = lambda x: x.scan_id[0])

        # for idx, hill in enumerate(self.finished_hills):

        #     for hill2 in self.finished_hills[idx+1:]:

        #         if hill.scan_id[-1] == hill2.scan_id[0]:
                    
        #             #self.crosslinked_hills.append(ready_hill(hill.mz + hill2.mz, hill.intensity + hill2.intensity, hill.scan_id + hill2.scan_id, hill.mass + hill2.mass)) # добавить сшивку, счетчик сшитых, проверку на массу и косинусную корреляцию
        #             self.crosslinked_hills.append(ready_hill(hill.mz + hill2.mz, hill.intensity + hill2.intensity, hill.scan_id + hill2.scan_id, hill.mass + hill2.mass))
        #             crosslink_counter += 1

        #         else:

        #             self.crosslinked_hills.append(hill)

        i = 0 
        ini_len = len(self.finished_hills)   

        while i < ini_len:
            
            hill = self.finished_hills[i]
            j = i + 1

            while j < ini_len:

                hill2 = self.finished_hills[j]

                #if hill.scan_id[-1] == hill2.scan_id[0]:
                if abs(hill.scan_id[-1] - hill2.scan_id[0]) <= 1:
                    #crosslink_counter2 += 1
                    if abs(hill.mz - hill2.mz)/hill.mz <= mass_accuracy * 1e-6:

                        self.finished_hills[i] = ready_hill(intensity=hill.intensity + hill2.intensity, scan_id=hill.scan_id + hill2.scan_id, mass=hill.mass + hill2.mass)
                        del self.finished_hills[j]
                        ini_len -= 1
                        crosslink_counter += 1
                elif hill2.scan_id[0] > hill.scan_id[-1] +1:
                    break

                j += 1

            i += 1


        #print(crosslink_counter)
        #print(crosslink_counter2)

    

    def sort_finished_hills(self):
        self.finished_hills = sorted(self.finished_hills, key = lambda x : x.mz)

    
    
    def check_its_ready(self, id_real, check_degree, min_length):

        mask_to_del = [True] * self.mz_array.size
        for i in range(self.mz_array.size)[::-1]:
            if id_real - self.scan_id[i][-1] > check_degree:

                tmp_ready_hill = ready_hill(intensity = self.intensity.pop(i), 
                                            scan_id = self.scan_id.pop(i), 
                                            mass = self.mass_array.pop(i), 
                                            )
                mask_to_del[i] = False
                
                if len(tmp_ready_hill.scan_id) >= min_length:
                    self.finished_hills.append(tmp_ready_hill)
                    
        self.mz_array = self.mz_array[mask_to_del]
        

    def push_left(self, min_length):
        mask_to_del = [True] * self.mz_array.size
        for i in range(self.mz_array.size)[::-1]:

            tmp_ready_hill = ready_hill(intensity = self.intensity.pop(i), 
                                        scan_id = self.scan_id.pop(i), 
                                        mass = self.mass_array.pop(i), 
                                        )
            mask_to_del[i] = False

            if len(tmp_ready_hill.scan_id) >= min_length:
                self.finished_hills.append(tmp_ready_hill)
                    
        self.mz_array = self.mz_array[mask_to_del] 
        
    def get_nearest_value(self, value, mask):
        #return min(self.mz_array[mask], key=lambda x: abs(x - value))
        #return min(np.abs(self.mz_array[mask]-value))
        return np.argmin(np.abs(self.mz_array[mask]-value))
    
    def newid(self, nearest, mask):
        
        #cnt = 0
        #i1 = 0
        
        # while True:
        #     if not mask[i1]:
        #         cnt += 1
        #     i1 += 1
            
        #     if i1 >= nearest:
        #         break
                   
        cnt2 = nearest -1 - sum(mask[:nearest-1])
        #print(cnt-cnt2)
            
        return cnt2 + nearest
        
    
    def push_me_to_the_peak(self, next_peak, diff, min_length):
        
        next_mz_array = next_peak.next_mz_array
        next_intensity_array = next_peak.next_intensity_array
        next_scan_id = next_peak.next_scan_id
        
        self.check_its_ready(id_real=next_scan_id, check_degree = 2, min_length = min_length)
        
        mask = [True]  * (len(self.mz_array))
        tmp1 = []
        tmp2 = []

        #FIXME
        #изменить метод добавления в пик. Применять среднее значение по предыдущим объектам
        for idx, i in enumerate(next_mz_array):
            nearest = self.get_nearest_value(i, mask)
            nearest_id = self.newid(nearest, mask)
            #nearest_id = np.nonzero(self.mz_array == nearest)[0][0]
            tmp_diff = abs(self.mz_array[nearest_id] - i) / i
            if tmp_diff <= diff * 1e-6:
                tmp1.append([nearest_id, idx, tmp_diff])


        
        tmp1_nearest_id_arr = np.array([ x[0] for x in tmp1])
        tmp1_idx_arr = np.array([x[1] for x in tmp1])
        tmp1_diff_arr = np.array([ x[2] for x in tmp1])
        
        sort_list = np.argsort(tmp1_diff_arr) #try different kinds
        tmp1_nearest_id_arr =  tmp1_nearest_id_arr[sort_list]
        tmp1_idx_arr = tmp1_idx_arr[sort_list]
        tmp1_diff_arr = tmp1_diff_arr[sort_list]
        
        saved_index = set()

        while tmp1:
            
            tmp_id = tmp1_idx_arr[0]
            
            if tmp1_diff_arr.size == 0:
                break

            if tmp1_diff_arr[0] > diff * 1e-6:
                break

            tmp2.append((tmp1_nearest_id_arr[0], tmp1_idx_arr[0]))

            #tmp_id = [x[2] for x in tmp1].index(min(x[2] for x in tmp1))
            # tmp2.append(tmp1[tmp_id])
            
            saved_index.add(tmp1_idx_arr[0])

            mask[tmp2[-1][0]] = False
            if sum(mask):
                # del tmp1[tmp_id]
                tmp1_nearest_id_arr = tmp1_nearest_id_arr[1:]

                tmp1_idx_arr = tmp1_idx_arr[1:]

                tmp1_diff_arr = tmp1_diff_arr[1:]

                if tmp1_diff_arr.size == 0:
                    break

                if tmp1_idx_arr[0] in saved_index:

                    for idx, element in enumerate(tmp1_idx_arr):
                        

                        if element in saved_index:

                        # if element[0] == tmp2[-1][0]:

                            nearest = self.get_nearest_value(element, mask)
                            nearest_id = self.newid(nearest, mask)
                            #element[0] = np.nonzero(self.mz_array == nearest)[0][0]
                            tmp1_nearest_id_arr[idx] = nearest_id
                            
                            tmp1_diff_arr[idx] = abs(self.mz_array[nearest_id] - element)/element
                        else:
                            break
                    sort_list = np.argsort(tmp1_diff_arr, kind='quicksort') #try different kinds
                    tmp1_nearest_id_arr =  tmp1_nearest_id_arr[sort_list]
                    tmp1_idx_arr = tmp1_idx_arr[sort_list]
                    tmp1_diff_arr = tmp1_diff_arr[sort_list]

            else:
                break
        

        #print(self.mass_array[1])


        for i, idx in tmp2:
            #FIXME
            #self.mz_array[i] = (self.mz_array[i] + next_mz_array[idx])/2
            self.scan_id[i].append(next_scan_id)
            self.intensity[i].append(next_intensity_array[idx])
            self.mass_array[i].append(next_mz_array[idx])
            self.mz_array[i] = np.mean(self.mass_array[i][-3:])

            
            
        added = set(x[1] for x in tmp2)
        mask2 = [(False if i in added else True) for i in range(len(next_mz_array))]
        next_mz_array_size = next_mz_array[mask2].size
        self.mz_array = np.append(self.mz_array, next_mz_array[mask2])
        
        for i in range(next_mz_array_size):
            self.scan_id.append([next_scan_id, ])
            self.intensity.append([next_intensity_array[mask2][i], ])
            self.mass_array.append([next_mz_array[mask2][i], ])

    def cutting_down(self, intensity_propotion):

    #result = []

        for idx, peak in enumerate(self.finished_hills):
            
            max_intensity_propotion = peak.max_intensity * intensity_propotion
            #FIXME try "and" 

            if (peak.intensity[0] >= max_intensity_propotion and peak.intensity[-1] >= max_intensity_propotion):

                del self.finished_hills[idx]
        


    #self.finished_hills = result

class feature:
    
    def __init__(self, finished_hills, each):

        self.charge = each[1][0][1]
        #self.mz = finished_hills[each[0]].mz
        self.mz = np.median(finished_hills[each[0]].mass)
        self.neutral_mass = self.mz * self.charge - 1.0073 * self.charge

        self.isotopes_numb = len(each[1])
        
        #self.scan_numb = finished_hills[each[0]].scan_id[-1] - finished_hills[each[0]].scan_id[0]
        self.scan_numb = len(finished_hills[each[0]].scan_id)

        self.id_for_scan = finished_hills[each[0]].intensity.index(max(finished_hills[each[0]].intensity))
        self.intensity = finished_hills[each[0]].max_intensity
        
        self.scan_id = finished_hills[each[0]].scan_id[self.id_for_scan]
        #self.scan_id = finished_hills[each[0]]
        self.RT = self.scan_numb
        self.sulfur = (1 if each[2] else 0)
