# TEST_CASE=tests/ti16.txt
# TEST_CASE=tests/utopia.txt
TEST_CASE=tests/Simplicity_Test_Case.txt

mkdir fast
glpsol -m src/osemosys_fast.txt -d $TEST_CASE -o fast.sol --xcheck
rm results/SelectedResults.csv
cp -f results/*.csv fast

mkdir short
glpsol -m src/osemosys_short.txt -d $TEST_CASE -o short.sol --xcheck
rm results/SelectedResults.csv
cp -f results/*.csv short

mkdir long
glpsol -m src/osemosys.txt -d $TEST_CASE -o long.sol --xcheck
rm results/SelectedResults.csv
cp -f results/*.csv long

diff fast long -wy > fast_long.diff
diff short long -wy > short_long.diff
diff fast short -wy > fast_short.diff
