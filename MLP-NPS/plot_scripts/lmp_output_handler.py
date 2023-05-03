#!/usr/bin/env python

import os
import numpy as np
from file_iterator import read_iter
from ase.atoms import Atoms

class lmp_output_handler(object):
    """ DOC """
    def __init__(self):
        self.path = os.getcwd()
    
    def _read_file(self,infilename):
        f = open(infilename,'r')
        lines = f.readlines()
        f.close()
        return(lines)

    def read_per_atom_dump(self,infile,sort=True,select=None):
        #TODO: fixed structure adjust to below if below is faster
        lines = self._read_file(infile)
        sample_step, data = [], np.array([])
        for i in range(0,len(lines)):
            if (lines[i].find('ITEM: TIMESTEP') != -1):
                sample_step.append(int(lines[i+1]))
            if (lines[i].find('ITEM: NUMBER OF ATOMS') != -1):
                natoms = int(lines[i+1])
            if (lines[i].find('ITEM: ATOM') != -1):
                nval = len(lines[i].split())-2
                step = []
                for j in range(0,natoms):
                    inner = []
                    for k in range(0,nval):
                        inner.append(float(lines[i+j+1].split()[k]))
                    step.append(inner)
                step = np.array(step)
                if (select != None):
                    step = step[np.where(step[:,select[0]]==select[1])[0],:]
                if (sort == True):
                    step = self.__sort_by_id(step)
                if (data.size == 0):
                    data = step
                    NoSel = data[:,0].size
                else:
                    data = np.vstack((data,step))
        return(np.array(sample_step),data,NoSel)

    def __sort_by_id(self,array):
        ind = np.argsort(array[:,0])
        new_array = np.zeros((array[:,0].size,array[0,:].size))
        new_array[:,:] = array[ind,:]
        return(new_array)

    def lmp_dump_timestep(self,infile):
        line_obj = read_iter(infile)
        l_unit = 9+int(line_obj.get_lines(3,4)[0]) #per atoms (atoms length)
        select = range(0,line_obj.nlines/l_unit) #TIMESTEP blocks
        tstep = []
        for i in select:
            lines = line_obj.get_lines(i*l_unit,(i+1)*l_unit) #grep step lines
            tstep.append(int(lines[1]))
        return(tstep)

    def lmp_per_atom_dump_2_nparray_list(self,infile,sort=True):
        ''' general function to read per-atom dump files of LAMMPS
            input:  infile: filename of dump file to be read
                    sort:   sort atom arrays for the first output
                            argument (makes sense if that is atom id)
            output: tstep:  list of timesteps found in LAMMPS output
                    out:    list containting data of each timestep 
                            as a numpy array
        '''
        line_obj = read_iter(infile)
        natom = int(line_obj.get_lines(3,4)[0]) #per atoms (atoms length)
        l_unit = natom+9 #TIMESTEP block length
        select = range(0,line_obj.nlines/l_unit) #TIMESTEP blocks
        out,tstep = [],[]
        for i in select:
            step = []
            lines = line_obj.get_lines(i*l_unit,(i+1)*l_unit) #grep step lines
            tstep.append(int(lines[1]))
            for j in range(natom):
                l = j+9
                line = [float(x) for x in lines[l].split()]
                step.append(line)
            step = np.array(step)
            if (sort):
                step = self.__sort_by_id(step)
            out.append(step)
        return(tstep,out)

    def lammps_dump_file_2_ase(self,template,select=None,species={},remove_type=[]):
        #TODO: put wrapper around per-atom-dump file; do timing!
       #lines = self._read_file(template)
       #l_unit2 = int(lines[3])+9
        line_obj = read_iter(template)
        l_unit = int(line_obj.get_lines(3,4)[0])+9
        if select == None: #standard, get all snapshots
            select = range(0,int(line_obj.nlines/l_unit))
        out = []
        for i in select:
            box, geometry = self._lines_2_pos(line_obj.get_lines(i*l_unit,(i+1)*l_unit))
            geometry = self.__remove_type_from_geom(geometry,remove_type)
            out.append(self._transform_to_ase(box,geometry,element_dict=species))
        return(out)

    def _lines_2_pos(self,lines):
        """ puts lines into data for single snapshot """
        nparticles = int(lines[3]) # 3 is a fixed offset
        #get box - 5 is a fixed offset
        box = []
        for i in range(0,3):
            inner = [float(lines[5+i].split()[j]) for j in range(0,len(lines[5+i].split()))]
            if (len(inner) == 2):
                inner.append(0.0)
            box.append(inner)
        box_ase = self.__box2vector_ase(np.array(box))
        #get coordinates - 9 is a fixed offset
        geometry = np.zeros((nparticles,5))
        for i in range(0,nparticles):
            line = [float(x) for x in lines[i+9].split()]
            geometry[int(line[0]-1),:] = np.array(line)
        if ((geometry[:,0] == 0).nonzero()[0].size != 0):
            raise Exception("not all atoms found")
        return(box_ase,geometry)
    
    def __box2vector_ase(self,box):
        #[xlo xhi xy],[ylo yhi xz],[zlo zhi yz]'''
        xhilo = (box[0,1] - box[0,0]) - (box[0,2]**2)**0.5 - (box[1,2]**2)**0.5
        yhilo = (box[1,1] - box[1,0]) - (box[2,2]**2)**0.5
        zhilo = (box[2,1] - box[2,0])
        vector = np.zeros((3,3))
        vector = [[xhilo, 0.0, 0.0],[box[0,2],yhilo, 0.0],[box[1,2],box[2,2],zhilo]]
        return(vector)

    def __remove_type_from_geom(self,geometry,atyp=[]):
        for i in range(0,len(atyp)):
            geometry = geometry[np.where(geometry[:,1] != atyp[i])[0],:]
        return(geometry)

    def _transform_to_ase(self,lmp_box,geometry,element_dict={}):
        if (len(element_dict)==0):
            [element_dict.update({x:x}) for x in range(1,100)] #if no dict type == element
        e_numbers = []
        [e_numbers.append(element_dict[x]) for x in geometry[:,1]]
        a = Atoms(numbers=e_numbers,positions=geometry[:,2:5],cell=lmp_box,pbc=True)
        return(a)

    def read_lmp_time_average_global(self,infile):
        ''' read file and return each timestep as a np.array in a list '''
        #TODO: add timestep dependence - now it only reads first output
        lines = self._read_file(infile)
        offset, nentry, nrow, dat, count = 4, int(lines[3].split()[1]), int(len(lines[2].split())-2), [], 4
        while count < len(lines):
            inner = []
            for i in range(count,count+nentry):
                inner.append([float(x) for x in lines[i].split()[-nrow:]])
            dat.append(np.array(inner))
            count += nentry+1
        return(dat)

### outer funcs

def lammps_dump_file_2_ase(template,species={},remove_type=[],select=None):
    ''' function to wrap around lmp_output_handler object and fct 
        lammps_dump_file_2_ase
        input:  template  = name of lammps dump file (str)
                species   = dictionary sorting lammps atom type to 
                            an atomic number, 
                            i.e. {lmp_atomtype:atomic_number} (dict)
                remove_type = types to not read (i.e. shell particles) (list)
                select = snap shots to read (None == all snapshots)
        output: atoms  = ase atoms object
    
    '''
    lmp_out = lmp_output_handler()
    atoms = lmp_out.lammps_dump_file_2_ase(template,select,species,remove_type)
    return(atoms)


