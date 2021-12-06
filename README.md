# Ανάπτυξη Λογισμικού για Αλγοριθμικά Προβλήματα - 2021-2022
## 2η Προγραμματιστική Εργασία - Time-Series-Clustering-Search

g++ -lm -g -O3 -o bin/search search/src/* lib/*.cpp lsh/src/hashTable.cpp -I lib -I search/inc -I lsh/inc
./bin/search -i input/nasdaq2017_LQ.csv -q input/nasdaq2017_LQ_query.csv -o outputFiles/output -algorithm LSH