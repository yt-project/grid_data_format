#!/projects/skillman/local/bin/python

import os
import weakref
import numpy as na
import h5py as h5

from collections import \
    defaultdict
from string import \
    strip, \
    rstrip
from stat import \
    ST_CTIME

fields = []
translation_dict = {}
translation_dict['density'] = 'density'
translation_dict['total_energy'] = 'specific_energy'
translation_dict['velocity_x'] = 'velocity_x'
translation_dict['velocity_y'] = 'velocity_y'
translation_dict['velocity_z'] = 'velocity_z'
translation_dict['cell_centered_B_x'] = 'mag_field_x'
translation_dict['cell_centered_B_y'] = 'mag_field_y'
translation_dict['cell_centered_B_z'] = 'mag_field_z'

def parse_line(line, grid):
#    print line
    # grid is a dictionary
    splitup = line.strip().split()
    if "vtk" in splitup:
        grid['vtk_version'] = splitup[-1]
    elif "Really" in splitup:
        grid['time'] = splitup[-1]
    elif "DIMENSIONS" in splitup:
        grid['dimensions'] = na.array(splitup[-3:]).astype('int')
    elif "ORIGIN" in splitup:
        grid['left_edge'] = na.array(splitup[-3:]).astype('float64')
    elif "SPACING" in splitup:
        grid['dds'] = na.array(splitup[-3:]).astype('float64')
    elif "CELL_DATA" in splitup:
        grid["ncells"] = int(splitup[-1])
    elif "SCALARS" in splitup:
        field = splitup[1]
        grid['read_field'] = field
        grid['read_type'] = 'scalar'
    elif "VECTORS" in splitup:
        field = splitup[1]
        grid['read_field'] = field
        grid['read_type'] = 'vector'
        
def read_grid(filename):
    """ Read Athena legacy vtk file from single cpu """
    f = open(filename,'rb')
    grid = {}
    grid['read_field'] = None
    grid['read_type'] = None
    table_read=False
    line = f.readline()
    while line is not '':
        while grid['read_field'] is None:
            parse_line(line, grid)
            if grid['read_type'] is 'vector':
                break
            if table_read is False:             
                line = f.readline()
            if 'TABLE' in line.strip().split():
                table_read = True
            if len(line) == 0: break
        #    print line

        if len(line) == 0: break
        if na.prod(grid['dimensions']) != grid['ncells']:
            grid['dimensions'] -= 1
        if na.prod(grid['dimensions']) != grid['ncells']:
            print 'product of dimensions %i not equal to number of cells %i' % \
                  (na.prod(grid['dimensions']), grid['ncells'])
            raise TypeError
            
        if grid['read_type'] is 'scalar':
            grid[grid['read_field']] = \
                na.fromfile(f, dtype='>f4', count=grid['ncells']).reshape(grid['dimensions'],order='F')
            fields.append(grid['read_field'])
        elif grid['read_type'] is 'vector':
            data = na.fromfile(f, dtype='>f4', count=3*grid['ncells'])
            grid[grid['read_field']+'_x'] = data[0::3].reshape(grid['dimensions'],order='F')
            grid[grid['read_field']+'_y'] = data[1::3].reshape(grid['dimensions'],order='F')
            grid[grid['read_field']+'_z'] = data[2::3].reshape(grid['dimensions'],order='F')
            fields.append(grid['read_field']+'_x')
            fields.append(grid['read_field']+'_y')
            fields.append(grid['read_field']+'_z')
        else:
            raise TypeError
        grid['read_field'] = None
        grid['read_type'] = None
        line = f.readline()
        if len(line) == 0: break
    grid['right_edge'] = grid['left_edge']+grid['dds']*(grid['dimensions'])
    return grid

def write_to_hdf5(fn, grid):
    f = h5.File(fn,'a')

    ## --------- Begin level nodes --------- ##
    g = f.create_group('gridded_data_format')
    g.attrs['format_version']=na.float32(1.0)
    g.attrs['data_software']='athena'
    data_g = f.create_group('data')
    field_g = f.create_group('field_types')
    part_g = f.create_group('particle_types')
    pars_g = f.create_group('simulation_parameters')

    dle = grid['left_edge'] # True only in this case of one grid for the domain
    gles = na.array([grid['left_edge']])
    gdims = na.array([grid['dimensions']])
    glis = ((gles - dle)/grid['dds']).astype('int64')
    gris = glis + gdims

    # grid_left_index
    gli = f.create_dataset('grid_left_index',data=glis)
    # grid_dimensions
    gdim = f.create_dataset('grid_dimensions',data=gdims)

    levels = na.array([0]).astype('int64') # unigrid example
    # grid_level
    level = f.create_dataset('grid_level',data=levels)

    ## ----------QUESTIONABLE NEXT LINE--------- ##
    # This data needs two dimensions for now. 
    n_particles = na.array([[0]]).astype('int64')
    #grid_particle_count
    part_count = f.create_dataset('grid_particle_count',data=n_particles)

    # Assume -1 means no parent.
    parent_ids = na.array([-1]).astype('int64')
    # grid_parent_id
    pids = f.create_dataset('grid_parent_id',data=parent_ids)

    ## --------- Done with top level nodes --------- ##
    
    f.create_group('hierarchy')

    ## --------- Store Grid Data --------- ##
    
    g0 = data_g.create_group('grid_%010i'%0)
    for field in fields:
        name = field
        if field in translation_dict.keys():
            name = translation_dict[name]
        if not name in g0.keys(): 
            g0.create_dataset(name,data=grid[field])
    
    ## --------- Store Particle Data --------- ##

    # Nothing to do

    ## --------- Attribute Tables --------- ##

    pars_g.attrs['refine_by'] = na.int64(1)
    pars_g.attrs['dimensionality'] = na.int64(3)
    pars_g.attrs['domain_dimensions'] = grid['dimensions']
    try:
	pars_g.attrs['current_time'] = grid['time']
    except:
	pars_g.attrs['current_time'] = 0.0
    pars_g.attrs['domain_left_edge'] = grid['left_edge'] # For Now
    pars_g.attrs['domain_right_edge'] = grid['right_edge'] # For Now
    pars_g.attrs['unique_identifier'] = 'athenatest'
    pars_g.attrs['cosmological_simulation'] = na.int64(0)
    pars_g.attrs['num_ghost_zones'] = na.int64(0)
    pars_g.attrs['field_ordering'] = na.int64(0)
    pars_g.attrs['boundary_conditions'] = na.int64([0]*6) # For Now

    # Extra pars:
    pars_g.attrs['n_cells'] = grid['ncells']
    pars_g.attrs['vtk_version'] = grid['vtk_version']

    # Add Field Attributes
    for name in g0.keys():
        this_field = field_g.create_group(name)
        this_field.attrs['field_to_cgs'] = na.float64('1.0') # For Now

    # Add particle types
    # Nothing to do here

    # Add particle field attributes
    f.close()

import sys
if __name__ == '__main__':
    n = sys.argv[-1]
    n = n.split('.')
    fn = '%s.%04i'%(n[0],int(n[1]))
    grid = read_grid(fn+'.vtk')
    write_to_hdf5(fn+'.gdf',grid)
    
