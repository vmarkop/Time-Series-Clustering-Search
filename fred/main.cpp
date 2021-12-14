
#include <iostream>
#include "curve.hpp"
#include <vector>
#include "frechet.hpp"
#include <string>
#include <fstream>
int initialize_db(std::vector<Curve> &DB, std::ifstream &file);


int main() {
    const curve_size_t m = 120;
    const dimensions_t dimensions = 1;
    Curve c1(m, dimensions, "firstCurwa");
    Point p(1);
    p.set(1, 10);
    std::cout << p.repr() << std::endl;
    c1.push_back(p);
    Curve c2(m, dimensions, "secondCurwa");


    std::vector<Curve> DB;
    
    std::ifstream file("query_small_id");
    initialize_db(DB, file);


    Frechet::Continuous::Distance d = Frechet::Continuous::distance(DB[10], DB[1]);

    std::cout << "Distance is: " << d.value << std::endl;
    file.close();
}


int initialize_db(std::vector<Curve> &DB, std::ifstream &file) {
    std::string str;
    Curve *p2;
    int d;
    const curve_size_t m = 120;
    const dimensions_t dimensions = 1;
    while (std::getline(file, str)) {

        std::string s = str;
        std::string delimiter = " ";

        std::vector<Point> v;
        int i = 0;
        size_t pos = 0;
        std::string token;
        std::string item_id;
        while ((pos = s.find(delimiter)) != std::string::npos) {
            token = s.substr(0, pos);
            s.erase(0, pos + delimiter.length());

            if (i != 0) { // skip the item name
                Point p(1);
                p[0] = stoi(token);
                std::cout << p[0] << std::endl;
                v.push_back(p);
            }
            else {
                item_id = token;
            }
            i++;

        }
        d = i-1;
        // v.push_back(stoi(s));
        if (d == 0) break;

        Curve K(d-1, dimensions, "firstCurwa");
        for (Point k : v) {
            K.push_back(k);
        }
        DB.push_back(K);

    }
    return d;
}