#!/usr/bin/env python

# This python script runs a number of single iteration regression tests and diffs output with the expected output...
# Run the RegTest folder and pass the path to the relevant cactus executable as an argument on the command line (ex. runreg.py ../../../cactus)
# Note: if no differences exist, the regression test output file is deleted for convenience.
#
# The cactus binary must exist on your path.

import pytest

import os
import sys
import filecmp
import difflib
import subprocess


def _print_file_comparison(fn1, fn2):
    same = filecmp.cmp(fn1, fn2)

    if not same:
        print('Summary of differences between ' + fn1 + ' and ' + fn2 + ':')
        
        f = open(fn1,'r');
        f1 = f.readlines()
        f.close()
        
        f = open(fn2,'r');
        f2 = f.readlines()
        f.close()
        
        diff_out = difflib.unified_diff(f1, f2, fromfile=fn1, tofile=fn2)
        for line in diff_out:
            sys.stdout.write(line)

    else:
        print('No differences')
        os.remove(fn1)

    return same


@pytest.mark.parametrize("regression_name", ["RegTest1", "RegTest2", "RegTest3"])
def test_regtest(regression_name):

    # Run regression test
    print('Running regression test %s' % regression_name)
    command = ['cactus', regression_name + '.in']

    try:
        subprocess.check_output(command)
    except subprocess.CalledProcessError as e:
        print(e.output)

    # clean up standard output files which are meaningless for this calculation
    os.remove('%s_Param.csv' % regression_name)
    os.remove('%s_RevData.csv' % regression_name)
    os.remove('%s_TimeData.csv' % regression_name)

    # diff output
    fn1 = '%s_RegData.out' % regression_name
    fn2 = '%s_RegData_Ex.out' % regression_name

    same = _print_file_comparison(fn1, fn2)

    assert same, 'Differences between param output and gold standard.'
