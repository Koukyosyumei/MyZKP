nvcc -std=c++17 test_fr.cu -o test_fr
nvcc -std=c++17 test_sumcheck.cu -o test_sumcheck
./test_fr
./test_sumcheck
