#! /bin/bash
mkdir -p mesh
seeder
mkdir -p tracking
mkdir -p restart
bash test_stability.sh
bash test_order.sh
