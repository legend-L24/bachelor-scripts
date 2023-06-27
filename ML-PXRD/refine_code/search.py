import numpy as np 
import os 



def split_by_space(string):
    if string[-1] == '\n':
        string = string[:-1]
    string = string.split(' ')
    new_list = []
    for each in string:
        if each != '':
            new_list.append(each)
    return new_list


def lst2cif(lstfile,ciffile):
    with open(lstfile,'r') as lstf:
        lines = lstf.readlines()

        Atoms_index = []
        for ii in range(len(lines)):
            if lines[ii][:7] == ' Atoms:':
                Atoms_index.append(ii)
            elif lines[ii][:9] == ' Density:':
                end_index = ii
            elif lines[ii][:15] == ' New unit cell:':
                cell_unit = split_by_space(lines[ii+2])[1:7]
                print(cell_unit)
        #print(Atoms_index)
        start_index = Atoms_index[1]

        coord_raw = lines[start_index+3:end_index-1]
        coordinate = [] 

        print(len(coord_raw))
        for ii in range(int(len(coord_raw)/3)):
            atom = []
            atom += split_by_space(coord_raw[3*ii])[:1]
            atom += split_by_space(coord_raw[3*ii+1])[1:4]
            coordinate.append(atom)

    
    with open(ciffile,'w') as ciff:
        ciff.write('data_global\n')
        ciff.write('_space_group_IT_number 205\n')
        ciff.write("_symmetry_space_group_name_H-M    'P a -3'\n")
        ciff.write('_cell_length_a %s\n'%(cell_unit[0]))
        ciff.write('_cell_length_b %s\n'%(cell_unit[1]))
        ciff.write('_cell_length_c %s\n'%(cell_unit[2]))
        ciff.write('_cell_angle_alpha %s\n'%(cell_unit[3]))
        ciff.write('_cell_angle_beta %s\n'%(cell_unit[4]))
        ciff.write('_cell_angle_gamma %s\n'%(cell_unit[5]))
        ciff.write('loop_\n_space_group_symop_operation_xyz\n')
        ciff.write("\t'x,y,z' '1/2-x,-y,1/2+z'        '-x,1/2+y,1/2-z'\n")
        ciff.write("\t'1/2+x,1/2-y,-z'        'z,x,y' '1/2+z,1/2-x,-y'\n")
        ciff.write("\t'1/2-z,-x,1/2+y'        '-z,1/2+x,1/2-y'        'y,z,x'\n")
        ciff.write("\t'-y,1/2+z,1/2-x'        '1/2+y,1/2-z,-x'        '1/2-y,-z,1/2+x'\n")
        ciff.write("\t'-x,-y,-z'      '1/2+x,y,1/2-z' 'x,1/2-y,1/2+z'\n")
        ciff.write("\t'1/2-x,1/2+y,z' '-z,-x,-y'      '1/2-z,1/2+x,y'\n")
        ciff.write("\t'1/2+z,x,1/2-y' 'z,1/2-x,1/2+y' '-y,-z,-x'\n")
        ciff.write("\t'y,1/2-z,1/2+x' '1/2-y,1/2+z,x' '1/2+y,z,1/2-x'\n")
        ciff.write('loop_\n')
        ciff.write('_atom_site_label\n')
        ciff.write('_atom_site_fract_x\n')
        ciff.write('_atom_site_fract_y\n')
        ciff.write('_atom_site_fract_z\n')
        for atom in coordinate:
            ciff.write('\t%s\t%s\t%s\t%s\n'%(atom[0],atom[1],atom[2],atom[3]))

def sort_and_pick(numstr_list):
    ciffile = numstr_list[0]
    ciffile = ciffile[:-1].split('/')[-1]
    minimum = 1000000
    index = -1
    num_list = []
    for each in numstr_list[1:]:
        num_list.append(float(each[:-1]))
    #print('!!!!!!!')    
    #print(num_list)
    for ii in range(len(num_list)):
        if ii==0:
            minimum = num_list[ii]
            index = ii
        else:
            if num_list[ii] < minimum:
                minimum = num_list[ii]
                index = ii

    return ciffile, minimum, index



wdir = os.getcwd()
all_list = os.listdir(wdir)
if not os.path.exists('selected_and_refined'):
    os.mkdir('selected_and_refined')

error_list = []

for each in all_list: 
    if each[:5] == 'Error':
        error_list.append(each)

temp_storage = []
for each in error_list:
    with open(each,'r') as f:
        data = f.readlines()
    ciffile, minimum, index = sort_and_pick(data)

    
    print(ciffile,minimum,index)
    if (index==-1):
        print('no calculation')
    if minimum < 20.0:
        lstfile = 'job%s/refine_trial_%d.lst'%(each[5:],index) 
        ciffile = os.path.join('selected_and_refined',ciffile+str(minimum))

        

        try:
            lst2cif(lstfile,ciffile)
            print(ciffile,minimum)
            print('reached 20.0\% accuracy')
        except:
            print('good but no result')
    print()
