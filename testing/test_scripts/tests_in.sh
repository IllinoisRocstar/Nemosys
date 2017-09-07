#!/bin/bash

# setting source dir as determined by cmake (@...@)
source_dir=@CMAKE_SOURCE_DIR@

# finding directory with test data
test_dir=$(find ${source_dir} -name test_data)

# specifying directory with vol2meshTransfer exec as determined
# by cmake (@...@)
exec_dir=@CMAKE_BINARY_DIR@/bin

# check for diff_tol executable
diff_exec=$(find ${source_dir} -name diff_tol)
if [ ! -e ${diff_exec} ]
then
	echo "diff_tol executable does not exist"
	echo "Compile diff_tol.cpp to generate"
	exit 1
fi

# passing name of test case folder as command line arg
folder=${test_dir}/$1

cd $folder
rm -f $(ls *.out | grep -v REF) # removing existing output
ref_file=$(ls *REF*)
inp_file=$(ls *.inp)
out_file=$(ls *REF* | sed 's/\-REF//g')
log_file=vol2planeTransfer_testresults.txt
rm -f ../${log_file}
echo "Testing ${folder} ..." | tee -a ../${log_file}

# SET TOLERANCE
TOL="1e-6"

${exec_dir}/vol2planeTransfer ${inp_file}

if [ ! -e ${out_file} ]
then
	echo "No ${out_file} results from execution!"
	exit 1
else
	# compare output file with gold file
	${diff_exec} ${out_file} ${ref_file} ${TOL}
	if [ $? -ne 0 ]
	then
		echo "Test Failed: ${out_file} differs from ${ref_file}" | tee -a ../${log_file}
		printf "\n" >> ../${log_file}
		exit 1
	else
		echo "Test Passed: ${out_file} and ${ref_file} are the same" | tee -a ../${log_file}
	fi
fi
printf "\n" >> ../${log_file}
exit 0
