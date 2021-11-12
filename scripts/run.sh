make clean
make
cp 00DescriptionFilePatrBasicos.txt ../results
echo --------- Table 2 -----------
source Obtaining_table_2.sh
mv patrBasicos.csv ../results/for_table2.csv

echo --------- Tables 3 to 8 -----------
source Obtaining_fromTable_3_to_8_Algorithm_2.sh
mv patrBasicos.csv ../results/for_tables_3to8_Alg_2.csv
source Obtaining_fromTable_3_to_8_Algorithm_3.sh
mv patrBasicos.csv ../results/for_tables_3to8_Alg_3.csv
source Obtaining_fromTable_3_to_8_Algorithm_4.sh
mv patrBasicos.csv ../results/for_tables_3to8_Alg_4.csv

echo --------- Tables 9 y 10 -----------
source Obtaining_Table_9_and_10_Algorithm_5_with_d_0.sh
mv patrBasicos.csv ../results/for_tables_9to10_Alg_5_d_0.csv
source Obtaining_Table_9_and_10_Algorithm_5_with_d_1.sh
mv patrBasicos.csv ../results/for_tables_9to10_Alg_5_d_1.csv
source Obtaining_Table_9_and_10_Algorithm_5_with_d_2.sh
mv patrBasicos.csv ../results/for_tables_9to10_Alg_5_d_2.csv
source Obtaining_Table_9_and_10_Algorithm_5_with_d_3.sh
mv patrBasicos.csv ../results/for_tables_9to10_Alg_5_d_3.csv
source Obtaining_Table_9_and_10_Algorithm_5_with_d_5.sh
mv patrBasicos.csv ../results/for_tables_9to10_Alg_5_d_5.csv

echo --------- Table 11 -----------
source Obtaining_Table_11_Algorithm_6.sh
mv patrBasicos.csv ../results/for_table_11.csv