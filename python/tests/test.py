import os
import sys
import unittest

sys.path.append('.')

from lmfit2 import read_file, lmfit2_mp, write_file, main

from numpy.testing import assert_array_equal

class TestLmfit2(unittest.TestCase):

    def test_fitting_multiprocessing(self):
        print("Running multiprocessing test on 2 threads.")

        import multiprocessing as mp

        pool = mp.Pool(2)

        input_filename = 'tests/test.rawacf'
        output_filename = 'tests/output_mp.lmfit2'
        expected_filename = 'tests/expected.lmfit2'

        fitted_records = list()

        #first, read the file we want to fit
        records = read_file(input_filename)

        temp = [pool.apply_async(lmfit2_mp,args=(x,records[x])) for x in range(len(records))]

        fitted_records = [p.get() for p in temp]
        fitted_records.sort()
        fitted_records = [r[1] for r in fitted_records]

        write_file(output_filename, fitted_records)

        output_file = read_file(output_filename)
        expected_file = read_file(expected_filename)

        for item1,item2 in zip(output_file,expected_file):
            item1_key = item1.keys()
            item1_key.sort()
            item2_key = item2.keys()
            item2_key.sort()

            self.assertEqual(item1_key,item2_key)

            for key in item1_key:
                assert_array_equal(item1[key],item2[key])


    def test_fitting_singleprocessing(self):
        print("Running single processing test.")

        input_filename = 'tests/test.rawacf'
        output_filename = 'tests/output_sp.lmfit2'
        expected_filename = 'tests/expected.lmfit2'

        main(input_filename,output_filename)


        output_file = read_file(output_filename)
        expected_file = read_file(expected_filename)

        for item1,item2 in zip(output_file,expected_file):
            item1_key = item1.keys()
            item1_key.sort()
            item2_key = item2.keys()
            item2_key.sort()

            self.assertEqual(item1_key,item2_key)

            for key in item1_key:
                assert_array_equal(item1[key],item2[key])


if __name__ == '__main__':
    unittest.main()

    os.remove('tests/output_mp.lmfit2')
    os.remove('tests/output_sp.lmfit2')
