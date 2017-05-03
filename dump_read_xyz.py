"""
    An example script utilising the Config class to read in a configuration 
    and output in .xyz format

    USAGE: 

    $ python2.7 dump_read_xyz.py -f <path-to-lammps-output-file> 
                                 -tse <start-timestep> <end-timestep> 
                                 -re <read-every-n-timesteps>  
"""

from config_dump import *

import argparse

def main()

    ## Parse command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", type=str, 
                        help = "Name of LAMMPS dump file to be read", 
                        required = True)
    parser.add_argument("-tse", type=int, nargs = 2, 
                        help = "Start and end timesteps for averaging",
                        required = True)
    parser.add_argument("-re", type=int, 
                        help = "Read every n timesteps",
                        required = True)
    args = parser.parse_args() 

    fname = args.filename
    print fname
    type_index = {'O':1, 'H':2, 'C':3, 'H_C':4, 'C_C':5}
    mass_index = {1:16.0, 2:1.0, 3:12.0, 4:1.0, 5:12.0}


    t_start = args.tse[0]
    t_end = args.tse[1]
    every_n = args.re
    new_conf = Config()
    new_conf.new_type_list(type_index, mass_index)
    new_conf.open_file(fname)
    #new_conf.generate_density_profile("y", 200, "TMC")
    #new_conf.generate_density_profile("y", 200, "MPD")
    mol_read_flag = False
    timeflag = 0
    while not new_conf.read_config() and not timeflag:

        if (new_conf.timestep >= t_start and new_conf.timestep <= t_end 
            and not new_conf.timestep%every_n):  
            if not mol_read_flag:
                new_conf.construct_molecule_list()
#               mol_read_flag = True
            new_conf.remove_pbc_mollist(new_conf.mollist)
            new_conf.print_xyz("graphene_NEMD")
        if new_conf.timestep > t_end: timeflag = 1

if __name__ == "__main__":
    main()
