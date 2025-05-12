#! /usr/bin/env python
# encoding: utf-8
# Harald Klimach 2011
#
# Compilation will always the executable and execute the unit tests,
# if only compilations should be done, without executions, use --notests
import sys
import os
import glob

APPNAME = 'musubi'
VERSION = '1'

top = '.'
out = 'build'

def options(opt):
    opt.load('compiler_cxx')
    opt.add_option('--no_harvesting', action='store_true',
                   default=False,
                   help = 'Do not include harvesting files for compilation.',
                   dest = 'no_harvesting')
    opt.add_option('--stream', choices=['PUSH', 'PULL'],
                   help='Streaming approach to use, either PUSH or PULL')

    mus_opt = opt.add_option_group('External libraries to include')
    mus_opt.add_option('--with_ext_tdf', action='store_true',
                       help='add external library to compute thermodynamic factors.')
    mus_opt.add_option('--libdmapp', action='store_true',
                       help='to use the highly optimized GHAL-based DMAPP \
                       collective algorithms, if available.')

def subconf(conf):
    conf.setenv('')
    if not conf.options.no_harvesting:
        conf.env.build_hvs = True
    else:
        conf.env.build_hvs = False

    conf.env['with_ext_tdf'] = False
    if conf.options.with_ext_tdf:
        conf.env['with_ext_tdf'] = True

    if conf.options.stream:
        for key, confenv in conf.all_envs.items():
            confenv.COCOFLAGS.append('-D{0}'.format(conf.options.stream))

    if not conf.env.LIB_PRECICE and conf.env['with_ext_tdf']:
        conf.setenv('')
        conf.setenv('cenv')
        # load c++ compiler for external thermodynamicFactor
        conf.load('compiler_cxx')
        # KM:todo set debug CXX flags only for debug variant
        #conf.env['CXXFLAGS'] = ['-g']
        conf.check(lib='stdc++', uselib_store='STDCXX', mandatory=True)
        # Need to have the libraries in the Fortran environment for linking
        conf.all_envs[''].LIB_STDCXX = conf.env.LIB_STDCXX
        conf.all_envs[''].LIBPATH_STDCXX = conf.env.LIBPATH_STDCXX

    if conf.options.libdmapp:
        conf.setenv('')
        conf.env['libdmapp'] = True
        conf.check_fc(lib='dmapp', uselib_store='DMAPP', mandatory=False)

def configure(conf):
    from waflib import Logs
    conf.setenv('')

    # Compile flags are set in treelm!
    # Check in tem/wscript for the flags that are set, for the various
    # available targets.
    subconf(conf)

    # Avoid some warnings:
    if not conf.options.nowarn:
        for key, fenv in conf.all_envs.items():
            if fenv.FC_NAME == 'GFORTRAN':
                fenv.FCFLAGS.append('-Wno-unused-dummy-argument')
                fenv.FCFLAGS.append('-fmax-stack-var-size=131072')
            if fenv.FC_NAME == 'IFORT':
                stderrflag = len(fenv.FCFLAGS)
                for flag in range(len(fenv.FCFLAGS)-1):
                    if (fenv.FCFLAGS[flag] == '-warn') and (fenv.FCFLAGS[flag+1] == 'stderrors'):
                        stderrflag = flag
                if stderrflag < len(fenv.FCFLAGS)-1:
                    del(fenv.FCFLAGS[stderrflag+1])
                    del(fenv.FCFLAGS[stderrflag])

    Logs.warn('Musubi modified flags:')
    Logs.warn('Default flags: {0}'.format(' '.join(conf.all_envs[''].FCFLAGS)))
    Logs.warn('Debug   flags: {0}'.format(' '.join(conf.all_envs['debug'].FCFLAGS)))

    conf.setenv('ford', conf.env)
    conf.env.ford_mainpage = 'mus_mainpage.md'
    conf.env.fordurl_mus = 'https://apes-suite.github.io/musubi/'

def build(bld):
    from waflib.extras.utest_results import utests

    pp_sources = bld.path.ant_glob('source/compute/*.fpp')
    pp_sources += bld.path.ant_glob('source/bc/*.fpp')
    pp_sources += bld.path.ant_glob('source/intp/*.fpp')
    pp_sources += bld.path.ant_glob('source/scheme/*.fpp')
    pp_sources += bld.path.ant_glob('source/geom/*.fpp')
    pp_sources += bld.path.ant_glob('source/derived/*.fpp')
    pp_sources += bld.path.ant_glob('source/variables/*.fpp')
    pp_sources += bld.path.ant_glob('source/turbulence/*.fpp')
    pp_sources += bld.path.ant_glob('source/particles/*.fpp')
    pp_sources += bld.path.ant_glob('source/*.fpp')
    mus_sources = bld.path.ant_glob('source/*.f90',
                                       excl=['source/musubi.f90','source/*eNRTL*.f90'])
    mus_sources += bld.path.ant_glob('source/init/*.f90')
    mus_sources += bld.path.ant_glob('source/bc/*.f90')
    mus_sources += bld.path.ant_glob('source/compute/*.f90')
    mus_sources += bld.path.ant_glob('source/intp/*.f90')
    mus_sources += bld.path.ant_glob('source/scheme/*.f90')
    mus_sources += bld.path.ant_glob('source/geom/*.f90')
    mus_sources += bld.path.ant_glob('source/derived/*.f90')
    mus_sources += bld.path.ant_glob('source/variables/*.f90')
    mus_sources += bld.path.ant_glob('source/turbulence/*.f90')
    mus_sources += bld.path.ant_glob('source/particles/*.f90')
    mus_sources += pp_sources

    if bld.cmd != 'docu':

      # source files to compute thermodynamic factors
      extLib_objs=[]
      if bld.env.with_ext_tdf:
          ext_tf_source = bld.path.ant_glob('external/thermFactor/Density.cpp')
          ext_tf_source += bld.path.ant_glob('external/thermFactor/eNRTL.cpp')
          ext_tf_source += bld.path.ant_glob('external/thermFactor/initialize.cpp')
          ext_tf_source += bld.path.ant_glob('external/thermFactor/Matrixoperation.cpp')
          mus_sources += bld.path.ant_glob('source/mus_eNRTL_module.f90')

          bld(
            features = 'cxx',
            source = ext_tf_source,
            target = 'ext_tdf')
          extLib_objs.append('ext_tdf')
          extLib_objs.append('STDCXX')
      else:
          mus_sources += bld.path.ant_glob('source/mus_eNRTL_dummy.f90')

      if bld.env.libdmapp:
          extLib_objs.append('DMAPP')

      bld(
        features = 'coco fc',
        source   = mus_sources,
        target   = 'mus_objs')

      utest_sources = bld.path.ant_glob('utests/*_module.f90')
      utest_sources += bld.path.ant_glob('utests/*_module.fpp')

      bld(
        features = 'coco fc',
        source   = utest_sources,
        use      = ['tem_objs', 'aotus', 'mus_objs'],
        target   = 'utest_objs')

      bld(
        features = 'coco fc fcprogram',
        name     = 'musubi',
        source   = 'source/musubi.f90',
        use      = ['blas_objs', 'lapack_objs', 'tem_objs', bld.env.mpi_mem_c_obj,
                    'aotus', 'mus_objs']+extLib_objs,
        target   = bld.path.parent.find_or_declare('musubi'))

      if bld.env.build_hvs and not bld.options.no_harvesting:
          mus_hvs_sources = bld.path.ant_glob('source/mus_harvesting/*.fpp')
          mus_hvs_sources += bld.path.ant_glob('source/mus_harvesting/*.f90',
                                              excl='source/mus_harvesting/mus_harvesting.f90')
          bld(
              features = 'fc',
              source   = mus_hvs_sources,
              target   = 'mus_hvs_objs')
          bld(
              features = 'fc fcprogram',
              name     = 'mus_harvesting',
              source   = 'source/mus_harvesting/mus_harvesting.f90',
              use      = ['tem_objs', 'aotus', 'mus_objs', bld.env.mpi_mem_c_obj,
                          'mus_hvs_objs', 'base64']+extLib_objs,
              target = bld.path.parent.find_or_declare('mus_harvesting'))

      utests(bld, ['lapack_objs', 'aotus', 'tem_objs', bld.env.mpi_mem_c_obj,
                   'mus_objs','utest_objs']+extLib_objs, preprocessor='coco')

    else:
      from waflib.extras.make_fordoc import gendoc

      mpp = bld(
        features = 'includes coco',
        source   = pp_sources)

      mpp.post()
      mus_preprocessed = []
      for ppm in mpp.tasks:
        for f in ppm.outputs:
          mus_preprocessed.append(f)

      if not bld.env.fordonline:
        mus_preprocessed.append(bld.env.fordext_aotus)
        mus_preprocessed.append(bld.env.fordext_tem)

      tgt = bld.path.get_bld().make_node('docu/modules.json')
      bld.env.fordext_mus = tgt

      bld( rule = gendoc,
           source = mus_preprocessed,
           src_paths = [bld.path.find_node('source').abspath()],
           target = tgt,
           extern = ['aoturl = {0}'.format(bld.env.fordext_aotus),
                     'temurl = {0}'.format(bld.env.fordext_tem)
                    ],
           extern_urls = ['aoturl = {0}'.format(bld.env.fordurl_aotus),
                          'temurl = {0}'.format(bld.env.fordurl_tem)
                         ],
           mainpage = os.path.join(bld.top_dir, 'mus', 'mus_mainpage.md')
      )


# clean output files
# add different extension format to remove
# by extending outputfiles list using outfiles.extend(glob.glob('*.yourfileext')
def cleanoutput(ctx):
    outputfiles=[]
    outputfiles.extend(glob.glob('*.vtk'))
    #print outputfiles
    for output in outputfiles:
        os.remove(output)

