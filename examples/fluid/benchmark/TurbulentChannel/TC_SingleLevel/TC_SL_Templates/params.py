import os
apesFolder = os.getenv('HOME')+'/apes/'
# Production directory - default production directory is 'prod'
# comment out if you don't want user defined production directory
prod_dir = 'prod'
mail_address = None #'kannan.masilamani@dlr.de'
smtp_server = { 'tunnel' : {'host':'robin.inf.tu-dresden.de'} }

# name of the shepherd log file
shepherd_out = 'shepherd.log'

# name of the log and rror file of the clone and build function
clone_build_out = 'clone_build.log'
clone_build_err = 'clone_build_error.log'
grep_performance = False
mus_rev = None

run_seeder = False
submit_job = False
create_plot = False

import math
pi = math.pi

shepherd_jobs = []

shepherd_jobs.append(dict(executable=apesFolder+'seeder/build/seeder',
                          template='seeder.template',
                          extension='lua',
                          run_exec = run_seeder,
                          run_command = '',
                          params = [['y_plus',25, 50]],
                          additional_params = dict(MESH='mesh/',
                                                   l_h=2*pi,
                                                   w_h=2*pi),
                          create_subdir = ['mesh'],
                          prefix = 'sdr',
                          label = 'seeder'))

shepherd_jobs.append(dict(executable=apesFolder+'musubi/build/musubi',
                          template='musubi.template',
                          extension='lua',
                          run_exec = False,
                          run_command = 'mpirun --oversubscribe -np 8',
                          #params = [["wall_func", 'musker', 'reichardt'],
                          #params = [["wall_func", 'werner_wengle'],
                          #          ["les",'smagorinsky', 'vreman', 'wale'],
                          params = [["wall_func", 'musker'],
                                    ["les", 'vreman'],
                          ],
                          additional_params = dict(mesh='../mesh/', 
                                                   relaxation ='mrt',
                                                   layout='d3q19',
                                                   omega_lim=0.01),
                          create_subdir = ['tracking','restart' ],
                          prefix = 'mus',
                          depend = ['seeder'],
                          create_depend_path = True,
                          label = 'musubi',
                          mail = False))

shepherd_jobs.append(dict(executable=None,
                          template = 'printSimParam.template',
                          extension='lua',
                          run_exec = False,
                          run_command = '',
                          depend = ['seeder','musubi'],
                          run_last = False,
                          create_depend_path = True,
                          create_depend_params = True,
                          label = 'print',
                          mail = False))

shepherd_jobs.append(dict(executable=None,
                          template = 'plot_track.template',
                          extension='py',
                          run_exec = create_plot,
                          run_command = 'python',
                          depend = ['musubi'],
                          run_last = False,
                          create_depend_path = True,
                          create_depend_params = True,
                          label = 'gleaner',
                          mail = False))

shepherd_jobs.append(dict(executable=None,
                     template='cara.template',
                     extension='sh',
                     run_exec = submit_job,
                     run_command = 'sbatch',
                     additional_params = dict(nprocs=256),
                     create_dir = False,
                     create_depend_path = True,
                     create_depend_params = True,
                     depend = ['musubi'],
                     label = 'cara'))

