// ConsoleApplication6.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;

int main()
{
    ifstream fin("C:\\Users\\omjag\\Downloads\\connectivity.txt");
    double a, b, c, d, e, f, g, h, i, j, k, l;

    std::vector<double> elem;
    std::vector<double> con1;
    std::vector<double> con2;
    std::vector<double> con3;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;



    while (fin >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l)
    {
        elem.push_back(a);
        con1.push_back(b);
        con2.push_back(c);
        con3.push_back(d);
        //std::cout << a << " " << b << " " << c << " " << d << endl;
    }

    ifstream finn1("C:\\Users\\omjag\\Downloads\\coordinates.txt");
    double cor1, cor2, cor3;

    while (finn1 >> cor1 >> cor2 >> cor3)
    {
        x.push_back(cor2);
        y.push_back(cor3);
        //std::cout << y << " " << z << " 0" << endl;
    }

    double npoin,nelem, ip1, ip2, ip3, loc1_x, loc1_y, loc2_x, loc2_y, loc3_x, loc3_y, D;

    npoin = x.size();
    nelem = elem.size();

    //Initailization of global stifness matrix
    vector<vector<double>> lhspo(npoin, vector<double>(npoin, 0));

    //cout << lhspo.size() << "\n";
    //cout << x.size() << "\n";
    //cout << con1.size();

    //loc1_x = x[con1[10]];
    //cout << loc1_x;
    //loc1_x = x[con1[]];


    for (int i = 0; i < nelem; i++)
    {
        ip1 = con1[i];
        ip2 = con2[i];
        ip3 = con3[i];
//        cout << i << " \t";

//        loc1_x = x[con1[i]-1];
//        cout << loc1_x << "\n";


        //local nodes data
        loc1_x = x[con1[i]-1];
        loc1_y = y[con1[i]-1];

        loc2_x = x[con2[i]-1];
        loc2_y = y[con2[i]-1];

        loc3_x = x[con3[i]-1];
        loc3_y = y[con3[i]-1];

        D = loc1_x * (loc2_y - loc3_y) - loc1_y * (loc2_x - loc3_x) + (loc2_x * loc3_y - loc2_y * loc3_x);

        ip1 = ip1 - 1;
        ip2 = ip2 - 1;
        ip3 = ip3 - 1;



        //Element Matrix
        lhspo[ip1][ip1] = lhspo[ip1][ip1] + (((loc2_y - loc3_y) * (loc2_y - loc3_y) + (loc2_x - loc3_x) * (loc2_x - loc3_x)) / (2 * D));
        lhspo[ip1][ip2] = lhspo[ip1][ip2] + (((loc2_y - loc3_y) * (loc3_y - loc1_y) + (loc2_x - loc3_x) * (loc3_x - loc1_x)) / (2 * D));
        lhspo[ip1][ip3] = lhspo[ip1][ip3] + (((loc2_y - loc3_y) * (loc1_y - loc2_y) + (loc2_x - loc3_x) * (loc1_x - loc2_x)) / (2 * D));
        lhspo[ip2][ip1] = lhspo[ip2][ip1] + (((loc3_y - loc2_y) * (loc2_y - loc3_y) + (loc3_x - loc2_x) * (loc2_x - loc3_x)) / (2 * D));
        lhspo[ip2][ip2] = lhspo[ip2][ip2] + (((loc3_y - loc2_y) * (loc3_y - loc1_y) + (loc3_x - loc2_x) * (loc3_x - loc1_x)) / (2 * D));
        lhspo[ip2][ip3] = lhspo[ip2][ip3] + (((loc3_y - loc2_y) * (loc1_y - loc2_y) + (loc3_x - loc2_x) * (loc1_x - loc2_x)) / (2 * D));
        lhspo[ip3][ip1] = lhspo[ip3][ip1] + (((loc1_y - loc2_y) * (loc2_y - loc3_y) + (loc1_x - loc2_x) * (loc2_x - loc3_x)) / (2 * D));
        lhspo[ip3][ip2] = lhspo[ip3][ip2] + (((loc1_y - loc2_y) * (loc3_y - loc1_y) + (loc1_x - loc2_x) * (loc3_x - loc1_x)) / (2 * D));
        lhspo[ip3][ip3] = lhspo[ip3][ip3] + (((loc1_y - loc2_y) * (loc1_y - loc2_y) + (loc1_x - loc2_x) * (loc1_x - loc2_x)) / (2 * D));

    }


    
    ifstream finn2("C:\\Users\\omjag\\Downloads\\Boundary_Faces.txt");

    double aa, bb, cc, dd, ee, ff, nbface, len, xx, yy;
    std::vector<double> bface;
    std::vector<double> bcon1;
    std::vector<double> bcon2;
    std::vector<double> vb;


    
    while (finn2 >> aa >> bb >> cc >> dd >> ee >> ff)
    {
        bface.push_back(aa);
        bcon1.push_back(bb);
        bcon2.push_back(cc);
        vb.push_back(ee);
        //std::cout << y << " " << z << " 0" << endl;
    }
    
    nbface = bface.size();
    vector<vector<double>> rhspo(npoin, vector<double>(1, 0));
    

    //cout << nbface;
    //cout << rhspo[0][0];
    
    for (int i = 0; i < nbface; i++)
    {

        //local nodes data
        loc1_x = x[bcon1[i] - 1];
        loc1_y = y[bcon1[i] - 1];

        loc2_x = x[bcon2[i] - 1];
        loc2_y = y[bcon2[i] - 1];

        len = sqrt(pow((loc2_y-loc1_y),2) + pow((loc2_x - loc1_x), 2));

        xx = bcon1[i] - 1;
        yy = bcon2[i] - 1;

        
        //Element Matrix
        if (vb[i] == 1)
        {
        rhspo[xx][0] = rhspo[xx][0] + ((vb[i] * len) / 2);
        rhspo[yy][0] = rhspo[yy][0] + ((vb[i] * len) / 2);
        //cout << yy << " " << bcon2[i] << "\n";
        //cout << vb[i];
        }
    }
    

    // Cholesky Decomposition
    vector<vector<double>> L(npoin, vector<double>(npoin, 0));

    for (int k = 0; k < npoin; k++)
    {
        for (int i = 0; i < npoin; i++)
        {

            if (k == i)
            {
                double var1 = 0, var2;
                for (int j = 0; j < i; j++)
                {
                    var2 = L[k][j];
                    var1 = var1 + pow(var2, 2);
                }
                L[k][k] = L[k][k] + sqrt(lhspo[k][k] - var1);
            }
            else
            {
                double var1 = 0;
                for (int j = 0; j < i; j++)
                {
                    var1 = var1 + L[i][j] * L[k][j];
                }
                if (L[i][i] == 0)
                {
                    L[k][i] = 0;
                }
                else
                {
                    L[k][i] = L[k][i] + (lhspo[k][i] - var1) / L[i][i];
                }

            }


        }
    }





 /*
    for (int i = 0; i < npoin-1; i++)
    {
        for (int j = 0; j < npoin-1; j++)
        {
            //lhspo[i][j] = 0;
            cout << lhspo[i][j] << "\t";
        }
        cout << endl;
    }
  */


    // Printing Out A Matrix

    ofstream
    fw("C:\\Users\\omjag\\Downloads\\A_matrix.txt");
    std::ofstream::out;

    if (fw.is_open())
    {
        for (int i = 0; i < npoin; i++)
        {
            for (int j = 0; j < npoin; j++)
            {
                //lhspo[i][j] = 0;
                fw << lhspo[i][j] << "\t";
            }
            fw << "\n";
        }
        fw.close();
    }

    else cout << "Problem with opening file";


    // Printing Out B Matrix
    ofstream
        fw1("C:\\Users\\omjag\\Downloads\\B_matrix.txt");
    std::ofstream::out;

    if (fw1.is_open())
    {
        for (int i = 0; i < npoin; i++)
        {
            fw1 << rhspo[i][0] << "\t";
            fw1 << "\n";
        }
        fw1.close();
    }

    else cout << "Problem with opening file";


}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
