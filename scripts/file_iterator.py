#!/usr/bin/env python

import os
import numpy as np


class iter(object):
    """ DOC """
    def __init__(self,infile):
        self.nlines = self.__get_nlines(infile)
        self._open_file(infile)

    def _open_file(self,string):
        self.infile_string = string
        self.infile = open(string,'r')
        self.iterpos = 0

    def _close_file(self):
        self.infile.close()

    def get_lines(self,start,stop):
        if (self.iterpos > start):
            self._close_file()
            self._open_file(self.infile_string)
        while (self.iterpos != start):
            trash = self.infile.readline()
            self.iterpos += 1
        trash = None
        lines = []
        while (self.iterpos != stop):
            lines.append(self.infile.readline())
            self.iterpos += 1
        return(lines)
    
    def __get_nlines(self,string):
        with open(string, "r") as ins:
            nlines = 0
            for line in ins:
                nlines += 1
        ins.close()
        return(nlines)

    def _get_snap_list(self,string):
        nlines = self._get_nlines(string)
        self._open_file(string)
        line = self.__iter_file(3,4)
        self._close_file()
        self._open_file(string)
        #lines = self.__get_in_string(string)
        #self.nparticles = int(lines[3])
        self.nparticles = int(line[0])
        timesteps = []
        skip = 9 + self.nparticles
        n = 0
        while (n < nlines):
            line = self.__iter_file(n+1,n+2)
            timesteps.append(int(line[0]))
            #timesteps.append(int(lines[n+1]))
            n += skip
        self.timesteps = np.array(timesteps)
        #del lines
        #print sys.getsizeof(self.timesteps)
        #lines = None
        self._close_file()

def read_iter(infile):
    return(iter(infile))

