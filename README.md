# Ανάπτυξη Λογισμικού για Αλγοριθμικά Προβλήματα - 2021-2022
### 2η Προγραμματιστική Εργασία
### Κατακερματισμός και αναζήτηση για Χρονοσειρές στη C++

### Our team:
**Ioannis  Georgopoulos**, sdi1800026

**Vasileios Markopoulos**, sdi1800108

This project includes a program that, given a set of timeseries, finds the nearest neighbour of each, and a program that clusters these timeseries in discrete groups. The timeseries can be represented either as **vectors**, or as **curves**. The search is performed by utilizing probabilistic algorithms *LSH* and *Hypercube Projection* for vectors, and *Frechet*'s variation of the *LSH* algorithm for curves. Clustering can use exhaustive search (*Lloyd's algorithm*), or reverse range search with the probabilistic algorithms mentioned above.

This project has been developed using [Git](https://github.com/vmarkop/Time-Series-Clustering-Search).

### Timeseries Representation
Representing timeseries as *vectors* belonging to the Euclidean Space R^d is straight-forward, by using C++ vectors with d elements, where each element is a floating point coordinate.

On the other hand, *curves* are a bit trickier to convey. In order to maintain functionality with some core functions of the project, 2-D curves in this program are also projected into C++ vectors. For instance, a 2-D curve containing 730 2-D points, each with coordinates (x,y), would be projected into a vector of size 1460, in the form `(x1,y1,x2,y2,...,x1460,y1460)`.

Finally, 1-D curves for the Continuous timeseries are represented with vectors as well, without taking the y coordinate into consideration.

## Project Structure
- **bin**
  - exec files
- **cluster**
  - inc
  - src
- **hypercube**
  - inc
  - src
- **input**
- **lib**
  - Fred
  - utils
- **lsh**
  - inc
  - src
- makefile
- **outputFiles**
- **search**
  - inc
  - src
- **unit**
- README

The **bin** folder contains the executable files created by the makefile.

The **cluster**, and **search** folders contain the main functions for their respective programs, as well as headers.

The **hypercube** and **lsh** folders contain functions implemented in Project1 which are used in both programs. The LSH folder, in particular, holds the source and header files for the new Frechet algorithm with Discrete and Continuous metrics.

The **input**  folder ideally contains input files used by the application, to reduce clutter.

The **outputFiles** folder ideally contains output files produced by the application, to reduce clutter.

The **lib** folder contains two lib folders:
- The **Fred** folder contains files from the Fred library, used solely to calculate Continuous Frechet Distance.
- The **utils** folder contains source and header files for utilities (functions, classes and structs) used in both programs.

The **unit** folder contains the source file for the unit tests.

## Compilation
From the main project folder, the included *makefile* can be used.
- *make [search | cluster]* produces the executable given in the *bin* folder.
- *make unittest* produces the executable for the unit test in the *bin* folder.
- *make* (or *make all*) produces all three files in the *bin* folder.
- *make clean* deletes files in the bin folder, as well as output files produced in the *outputFiles* folder.

## Using the programs
Commands to use the 2 programs with default parameters:
- search: `./bin/search.out -i input/nasd_input.csv -q input/nasd_query.csv -o outputFiles/search_res.txt -algorithm <Alg> -metric <Metric>`

  - Alg can be: *LSH* or *Hypercube* or *Frechet*.
  - If Alg==Frechet, then Metric is required to be *discrete* or *continuous*.

- cluster:
`./bin/cluster.out -i input/nasd_input.csv -c input/cluster.conf -o outputFiles/cluster_res.txt -assignment <Assgn> -update <Updt> -silhouette`
  - Assignment can be: *Classic* or *LSH* or *Hypercube* or *LSH_Frechet*.
  - Update can be: *Mean_Frechet* or *Mean_Vector*. **NOTE!** don't forget the underscore!
  - Valid Assignment-Update method combinations are: Classic-Frechet, Classic-Vector, LSH-Vector, Hypercube-Vector, Frechet-Frechet.

## Input File Structure
For the input file, every input vector is represented by a line. The first word is the vector's ID (can be string or int), and it is followed by numbers (integers or floating point numbers) separated by tabs, which represent the vector's value in each dimension. The vector's dimension is implied by the amount of numbers provided in the line (excluding the first ID string, if it was also an int). The query file follows the same pattern.

## Unit Testing
A small Unit Testing program has been implemented using the gtest library by Google.
- The Google Tests library must be installed before Unit Testing can commence.
- If not already installed, you may follow these steps:

```
sudo apt-get install libgtest-dev
sudo apt-get install cmake
cd /usr/src/gtest
sudo cmake CMakeLists.txt
sudo make
cd lib
sudo cp *.a /usr/lib
```
- After installing the library, compile the test using
`make unittest`,and run using `./bin/unit_test.out`.
