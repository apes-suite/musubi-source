import os
apesFolder = os.getenv('HOME')+'/apes/'
prod_dir = 'convergence_comp_MfrEq_PressEq'
loglevel = 'info'
mail_address = None

run_seeder = True
run_musubi = True
run_gleaner = True
run_gleaner2 = True
grep_performance = False

shepherd_jobs = []

shepherd_jobs.append(dict(executable=apesFolder+'/seeder/build/seeder',
                          template='seeder.template',
                          extension='lua',
                          run_exec = run_seeder,
                          run_command = '',
                          params = [[ 'nHeight',10,20,40,80 ]],
                          abort_failure = True,
                          create_subdir = ['mesh'],
                          prefix = 'sdr',
                          label = 'seeder'))

shepherd_jobs.append(dict(executable=apesFolder+'/musubi/build/musubi',
                          template='musubi.template',
                          extension='lua',
                          run_exec = run_musubi,
                          run_command = 'mpirun -np 12',
                          params = [[ 'omega', 1.1, 1.6, 1.9 ]],
                          additional_params = dict(MESH='../mesh/',),
                          abort_failure = False,
                          create_subdir = ['tracking','restart'],
                          prefix = 'mus',
                          depend = ['seeder'],
                          create_depend_path = True,
                          create_depend_params = True,
                          label = 'musubi'))

shepherd_jobs.append(dict(executable=None,
                          template='printParams.template',
                          extension='lua',
                          run_exec = False,
                          run_command = '',
                          abort_failure = True,
                          depend = ['seeder','musubi'],
                          create_depend_path = True,
                          create_depend_params = True,
                          label = 'print'))

shepherd_jobs.append(dict(executable=None,
                          template = 'plot_track.template',
                          extension='py',
                          run_exec = run_gleaner,
                          run_command = 'python3.7',
                          depend = ['musubi'],
                          abort_failure = False,
                          run_last = True,
                          create_depend_path = True,
                          create_depend_params = True,
                          label = 'plot_track'))

shepherd_jobs.append(dict(executable=None,
                          template = 'plot_l2error.template',
                          extension='py',
                          run_exec = run_gleaner2,
                          run_command = 'python3.7',
                          depend = ['seeder','musubi'],
                          abort_failure = True,
                          run_last = True,
                          create_depend_path = True,
                          create_depend_params = True,
                          label = 'plot_l2error'))
