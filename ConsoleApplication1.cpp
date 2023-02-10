#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>

using namespace std;

int main()
{
    ifstream fin("C:\\Users\\omjag\\Downloads\\connectivity.txt");
    double a, b, c, d, e, f, g, h, i, j, k, l, m;
    double elem_last;

    std::vector<double> elem;
    std::vector<double> con1;
    std::vector<double> con2;
    std::vector<double> con3;

    while (fin >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l)
    {
        elem.push_back(a);
        con1.push_back(b);
        con2.push_back(c);
        con3.push_back(d);
        //std::cout << a << " " << b << " " << c << " " << d << endl;

        //std::cout << "3 " << b - 1 << " " << c - 1 << " " << d - 1 << endl;
    }

    cout << elem[elem.size() -1 ];
    elem.push_back(con1);
    
    //elem.insert(elem.end(), con1.begin(), con1.end());
   

    cout << " " << elem[elem.size() - 1];

    return 0;
}
