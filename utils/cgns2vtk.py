#!/usr/bin/python

import os
import glob

# Descriptions for each file type
file_type = {0: 'volumetric',
             1: 'burning',
             2: 'nonburning',
             3: 'noninteracting',
             4: 'burn_surf',
             5: 'iburn_surf',
             6: 'volumetric_solid'}

# String pattern for each file
file_string = {0: 'fluid_*.cgns',  # volumetric
               1: 'ifluid_b_*.cgns',  # burning
               2: 'ifluid_nb_*.cgns',  # non-burning
               3: 'ifluid_ni_*.cgns',  # non-interacting
               4: 'burn*.cgns',
               5: 'iburn*.cgns',
               6: 'solid*.cgns'}

# Is surface?
surf_bool = {0: 0,
             1: 1,
             2: 1,
             3: 1,
             4: 0,
             5: 1,
             6: 0}

# Prefix string pattern for each file
file_prefix_string = {0: 'fluid',  # volumetric
                      1: 'ifluid_b',  # burning
                      2: 'ifluid_nb',  # non-burning
                      3: 'ifluid_ni',  # non-interacting
                      4: 'burn',
                      5: 'iburn_all',
                      6: 'solid'}


# Function for sorting filenames by time
def get_time(fname):
    dec_str = fname.split('.')[1]
    dec_str = dec_str.split('_')[0]

    exp_str = fname.split('.')[0][-2:]

    t = float('.' + dec_str) * float(10 ** int(exp_str))
    first, second = str(t).split('.')
    t = float(first + '.' + second[0])
    return t


# For each type of file
for itype in range(len(file_string)):

    # Create stitched result subdirectory
    os.system('rm -rf ' + file_type[itype])
    os.system('mkdir ' + file_type[itype])

    # Get and sort filenames by time
    file_list = glob.glob(file_string[itype])
    file_list = sorted(file_list, key=get_time)

    # Get number of partitions for each file type
    surf = surf_bool[itype]
    file_list_single = []
    if len(file_list) != 0:
        if surf:
            base_t = file_list[0].split('_')[2]
        else:
            base_t = file_list[0].split('_')[1]
        for file in file_list:
            if base_t in file:
                print(file)
                file_list_single.append(file)

    if len(file_list_single) != 0:
        npart = len(file_list_single)
        ntime = int(len(file_list) / npart)

        # For each time step
        for itime in range(ntime):
            file = file_list[itime * npart]

            surf = surf_bool[itype]

            if surf:
                base_t = file.split('_')[2]
            else:
                base_t = file.split('_')[1]
            prefix = file_prefix_string[itype]

            os.system('rocStitchMesh ' + '.' + ' ' + prefix + ' ' + base_t + ' ' + str(surf))

            list_of_files = glob.glob('./*.vtu')
            if len(list_of_files) != 0:
                latest_file = max(list_of_files, key=os.path.getctime)
                fname = file_type[itype] + '_' + str(itime)
                os.system('mv ' + latest_file + ' ' + file_type[itype] + '/' + fname + '.vtu')
