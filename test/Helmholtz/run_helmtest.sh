rm -rf print_testres.txt
./test_helmrouts3d
./test_hfmm3d
./test_hfmm3d_zkbig
./test_hfmm3d_vec
./test_hfmm3d_mps
mv print_testres.txt ../../print_testreshelm.txt
rm -rf fort.13
