#!/bin/bash

# target machine and folder to transfer the tar file
target=juqueen:apes/musubi
# name of tar file
tar_file=musubi_app.tar

# Make a tar file using necessaru source files
tar cvf ${tar_file} ./source/
tar rvf ${tar_file} ./treelm/source/
tar rvf ${tar_file} ./treelm/external/
tar rvf ${tar_file} ./treelm/utests/
tar rvf ${tar_file} ./treelm/config.lua
tar rvf ${tar_file} ./treelm/waf ./treelm/wscript
tar rvf ${tar_file} ./treelm/aotus/LuaFortran/
tar rvf ${tar_file} ./treelm/aotus/external/
tar rvf ${tar_file} ./treelm/aotus/sample/
tar rvf ${tar_file} ./treelm/aotus/source/
tar rvf ${tar_file} ./treelm/aotus/utests/
tar rvf ${tar_file} ./treelm/aotus/waf
tar rvf ${tar_file} ./treelm/aotus/wscript
tar rvf ${tar_file} ./utests/
tar rvf ${tar_file} default.coco
tar rvf ${tar_file} musubi.lua
tar rvf ${tar_file} wscript
tar rvf ${tar_file} waf
tar rvf ${tar_file} testsuite/performance/parallel/TGV_3D_single_level/

# copy to target machine
scp -r  ${tar_file} ${target}
