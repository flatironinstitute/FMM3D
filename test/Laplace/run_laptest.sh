rm -rf print_testres.txt
./test_laprouts3d
./test_rfmm3dpart
./test_rfmm3dpart_vec
mv print_testres.txt ../../print_testreslap.txt
rm -rf fort.13
