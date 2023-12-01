import os
apesFolder = os.getenv('HOME')+'/apes/'
# Production directory - default production directory is 'prod'
# comment out if you don't want user defined production directory
prod_dir = 'prod'
loglevel = 'INFO'
shepherd_jobs = []
run_seeder = False
run_seeder = True

run_musubi = False
run_musubi = True

run_hvs = True
run_hvs = False

run_glr = True
run_glr = False

run_coeff = True
run_coeff = False

submit_job = True
submit_job = False


shepherd_jobs.append(dict(executable=apesFolder+'seeder/build/seeder',
                          template='seeder.template',
                          extension='lua',
                          run_exec = run_seeder,
                          run_command = '',
                          params = [
                                    #['nHeight', 64, 128, 256, 512]
                                    ['nHeight', 64]
                                   ] ,
                          additional_params = dict(MESH='mesh/',
                                                   stl_file=os.getcwd()+'/cylinder.stl'),
                          create_subdir = ['mesh'],
                          prefix = 'sdr',
                          label = 'seeder'))

shepherd_jobs.append(dict(executable=apesFolder+'harvester/build/harvester',
                          template='harvester_mesh.lua',
                          extension='lua',
                          run_exec = False,
                          run_command = '',
                          create_subdir = ['harvest'],
                          depend = ['seeder'],
                          prefix = 'hvs',
                          label = 'hvs'))

shepherd_jobs.append(dict(executable=apesFolder+'musubi/build/musubi',
                          template='musubi.template',
                          extension='lua',
                          run_exec = run_musubi,
                          run_command = 'mpirun -np 4',
                          params = [[ "Re", 20]],
                          additional_params = dict(MESH='../mesh/'),
                          create_subdir = ['tracking','restart'],
                          prefix = 'mus',
                          depend = ['seeder'],
                          create_depend_path = True,
                          create_depend_params = True,
                          label = 'musubi'))

shepherd_jobs.append(dict(executable=apesFolder+'harvester/scripts/harvest_series.py',
                          template = 'params_hvs.py',
                          extension = 'py',
                          run_exec = run_hvs,
                          run_command = '',
                          input_option = '--config',
                          depend = ['musubi'],
                          #run_last = True,
                          #create_depend_path = True,
                          #create_depend_params = True,#["musubi_phase_shift","musubi_theta_eq"],
                          label = 'hvs'))

shepherd_jobs.append(dict(executable=apesFolder+'harvester/gleaner/gleaner.py',
                          template = 'params_plot.template',
                          extension = 'py',
                          run_exec = run_glr,
                          run_command = '',
                          input_option = '',
                          depend = ['seeder','musubi'],
                          create_depend_path = True,
                          create_depend_params = True,
                          label = 'glr'))

shepherd_jobs.append(dict(executable=None,
                          template = 'get_coeff_max.template',
                          extension = 'py',
                          run_exec = run_coeff,
                          run_command = 'python3',
                          input_option = '',
                          depend = ['musubi'],
                          create_depend_path = True,
                          create_depend_params = True,
                          label = 'coeff'))

shepherd_jobs.append(dict(executable=None,
                     template='horus.template',
                     extension='sh',
                     run_exec = submit_job,
                     run_command = 'sbatch',
                     params = [
                               ["nNodes", 2]
                              ],
                     create_dir = False,
                     create_depend_path = True,
                     depend = ['musubi'],
                     label = 'horus'))

