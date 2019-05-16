rm -rf print_testres.txt
./test_laprouts3d
./test_lfmm3d
./test_lfmm3d_vec
mv print_testres.txt ../../print_testreslap.txt
rm -rf fort.13
