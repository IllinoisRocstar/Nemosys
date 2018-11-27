#!/bin/bash
set -x
# setting source dir as determined by cmake (@...@)
source_dir=@CMAKE_SOURCE_DIR@

# finding directory with test data
test_dir=$(find ${source_dir} -name test_data)

# specifying directory with vol2meshTransfer exec as determined
# by cmake (@...@)
exec_dir=@CMAKE_BINARY_DIR@/bin

# check for diff_tol executable
diff_tol=$(find ${source_dir} -name diff_tol)
if [ ! -e ${diff_tol} ]; then
  echo "diff_tol executable does not exist"
  echo "Compile diff_tol.cpp to generate"
  exit 1
fi

# passing name of test case folder as command line arg
folder=${test_dir}/$1

# SET TOLERANCE
TOL="1e-6"

cd $folder

if [ $1 == "Refine_Quality_Transfer_test" ]; then
  rm -f driver_test.vtu refined_case0001.vtu meshQual.txt
  ref_mesh=$(ls *REF*.vtu)
  inp_file=$(ls *.json)
  out_mesh=refined_case0001.vtu
  ref_text=$(ls *REF*.txt)
  out_text=meshQual.txt
  log_file=Refine_Quality_Transfer_test.txt
  rm -f ../${log_file}
  echo "Testing ${folder} ..." | tee -a ../${log_file}
  ${exec_dir}/nemosysRun ${inp_file} | tee -a ../${log_file}
  if [ ! -e ${out_mesh} ] || [ ! -e ${out_text} ]; then
    echo "Either ${out_mesh} or ${out_text} did not result from execution!"
    exit 1
  else
    ${diff_vtu} ${out_mesh} ${ref_mesh} ${TOL} && 
    if [ $? -ne 0 ]; then
      echo "Test Failed: ${out_mesh} differs from ${ref_mesh}" | tee -a ../${log_file}
      printf "\n" >> ../${log_file}
      exit 1
    else
      echo "Test Passed: ${out_mesh} and ${ref_mesh} are the same" | tee -a ../${log_file}
    fi
    ${diff_tol} ${out_text} ${ref_text} ${TOL} &&
    if [ $? -ne 0 ]; then
      echo "Test Failed: ${out_text} differs from ${ref_text}" | tee -a ../${log_file}
      printf "\n" >> ../${log_file}
      exit 1
    else
      echo "Test Passed: ${out_text} and ${ref_text} are the same" | tee -a ../${log_file}
    fi
    exit 0
  fi
  printf "\n" >> ../${log_file}
elif [ $1 == "MeshGen_UnifRefine_test" ]; then
  rm -f $(ls | egrep -v 'stl|REF|json')
  ref_mesh=$(ls *REF*.vtu)
  inp_file=$(ls *.json)
  out_mesh=refined_uniform_hinge.vtu
  log_file=MeshGen_UnifRefine_test.txt
  rm -f ../${log_file}
  echo "Testing ${folder} ..." | tee -a ../${log_file}
  ${exec_dir}/nemosysRun ${inp_file} | tee -a ../${log_file}
  if [ ! -e ${out_mesh} ]; then
    echo "${out_mesh} did not result from execution!"
    exit 1
  else
    ${diff_vtu} ${out_mesh} ${ref_mesh} ${TOL} && 
    if [ $? -ne 0 ]; then
      echo "Test Failed: ${out_mesh} differs from ${ref_mesh}" | tee -a ../${log_file}
      printf "\n" >> ../${log_file}
      exit 1
    else
      echo "Test Passed: ${out_mesh} and ${ref_mesh} are the same" | tee -a ../${log_file}
    fi
    exit 0
  fi
  printf "\n" >> ../${log_file}
fi
exit 1



