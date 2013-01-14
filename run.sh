cmake .
echo "This test suite will fail if you do not have Hermes, or Valgrind, etc.."
echo "CMake generation done."
echo "Run memory-leaks tests? (Long) [y/n]"
read ans
if [ "$ans" = "y" ]; then
  memory_start_time= `date +%s`
	echo "Processing memory leaks tests..."
  cd memory-leaks
  cd 01-memory-simple
  echo "Building the test..."
  make
  echo "Test built, running..."
  valgrind --leak-check=full --log-file=../01-memory-simple-ValgrindLogfile ./01-memory-simple 3 3
  echo "Valgrind output '01-memory-simple-ValgrindLogfile' available in memory-leaks/"
  cd ../02-memory-adapt
  echo "Building the test..."
  make
  echo "Test built, running..."
  valgrind --leak-check=full --log-file=../02-memory-adapt-ValgrindLogfile ./02-memory-adapt
  echo "Valgrind output '02-memory-adapt-ValgrindLogfile' available in memory-leaks/"
  cd ../03-memory-transient-adapt
  echo "Building the test..."
  make
  echo "Test built, running..."
  valgrind --leak-check=full --log-file=../03-memory-transient-adapt-ValgrindLogfile ./03-memory-transient-adapt
  echo "Valgrind output '03-memory-transient-adapt-ValgrindLogfile' available in memory-leaks/"
  echo Memory leaks tests runtime - $(expr `date +%s` - $memory_start_time) s
  echo "Memory leaks tests - Done."
  cd ../..
fi
echo "Run performance tests? (Very long) [y/n]"
read ans
if [ "$ans" = "y" ]; then
  perf_start_time= `date +%s`

	echo "Processing performance tests..."
  echo "In the meantime:"
  echo "     Valgrind(memcheck):   http://valgrind.org/docs/manual/mc-manual.html"
  echo "     Callgrind:   http://valgrind.org/docs/manual/cl-manual.html"
  echo "     Cachegrind:   http://valgrind.org/docs/manual/cg-manual.html"
  echo "     DHAT:   http://valgrind.org/docs/manual/dh-manual.html"
  echo "     Massif:   http://valgrind.org/docs/manual/ms-manual.html"
  echo ""
  cd performance
  cd 01-performance-simple
  echo "Building the test..."
  make
  echo "Test built, running..."
  rm -f cachegrind.out.*
  valgrind --log-file=temp --tool=cachegrind ./01-performance-simple 3 3
  cg_annotate cachegrind.out.* > ../01-performance-simple-CachegrindLogfile
  echo "Cachegrind output '01-performance-simple-CachegrindLogfile' available in performance/"
  valgrind --log-file=../01-performance-simple-DHATLogfile --tool=exp-dhat ./01-performance-simple 3 3
  echo "DHAT output '01-performance-simple-DHATLogfile' available in performance/"
  rm -f massif.out.*
  valgrind --log-file=temp --tool=massif ./01-performance-simple 3 3
  ms_print massif.out.* > ../01-performance-simple-MassifgrindLogfile
  echo "Massif output '01-performance-simple-MassifgrindLogfile' available in performance/"


  cd ../02-performance-adapt
  echo "Building the test..."
  make
  echo "Test built, running..."
  rm -f cachegrind.out.*
  valgrind --log-file=temp --tool=cachegrind ./02-performance-adapt
  cg_annotate cachegrind.out.* > ../02-performance-adapt-CachegrindLogfile
  echo "Cachegrind output '02-performance-adapt-CachegrindLogfile' available in performance/"
  valgrind --log-file=../02-performance-adapt-DHATLogfile --tool=exp-dhat ./02-performance-adapt
  echo "DHAT output '02-performance-adapt-DHATLogfile' available in performance/"
  rm -f massif.out.*
  valgrind --log-file=temp --tool=massif ./02-performance-adapt
  ms_print massif.out.* > ../02-performance-adapt-MassifgrindLogfile
  echo "Massif output '02-performance-adapt-MassifgrindLogfile' available in performance/"

  cd ../03-performance-transient-adapt
  echo "Building the test..."
  make
  echo "Test built, running..."
  rm -f cachegrind.out.*
  valgrind --log-file=temp --tool=cachegrind ./03-performance-transient-adapt
  cg_annotate cachegrind.out.* > ../03-performance-transient-adapt-CachegrindLogfile
  echo "Cachegrind output '03-performance-transient-adapt-CachegrindLogfile' available in performance/"
  valgrind --log-file=../03-performance-transient-adapt-DHATLogfile --tool=exp-dhat ./03-performance-transient-adapt
  echo "DHAT output '03-performance-transient-adapt-DHATLogfile' available in performance/"
  rm -f massif.out.*
  valgrind --log-file=temp --tool=massif ./03-performance-transient-adapt
  ms_print massif.out.* > ../03-performance-transient-adapt-MassifgrindLogfile
  echo "Massif output '03-performance-transient-adapt-MassifgrindLogfile' available in performance/"
  echo Performance leaks tests runtime - $(expr `date +%s` - $memory_start_time) s
  echo "Performance leaks tests - Done."
  cd ../..
fi
echo "Run visualization tests? (Short) [y/n]"
read ans
if [ "$ans" = "y" ]; then
	echo "Processing visualization tests..."
  cd visualization/views
  make
  ./01-all
  echo "Outputs saved in visualization/views/*.bmp"
  echo "Visualization tests - Done."
  cd ../../..
fi
echo "Run inner-funcionality tests? (Very short) [y/n]"
read ans
if [ "$ans" = "y" ]; then
	echo "Processing internal tests..."
  cd inner-functionality
  make -j4
  ctest -j4
  echo "Internal-functionality tests - Done."
  cd ..
fi
echo "Run calculations tests? [y/n]"
read ans
if [ "$ans" = "y" ]; then
	echo "Processing calculations tests..."
  cd calculations
  make -j4
  ctest -j4
  echo "Calculations tests - Done."

fi
echo "Quitting..."
