for i in {0..4}; do echo "--> ${i} <--" > files{i}.out; done

cat ${inputdir}/*.out > one_file.out
