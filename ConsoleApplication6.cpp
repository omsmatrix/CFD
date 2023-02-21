// ConsoleApplication6.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
using namespace std;

int main()
{
    ifstream fin("C:\\Users\\omjag\\Downloads\\Cy_2_Connectivity.txt");
    double a, b, c, d, e, f, g, h, i, j, k, l;

    std::vector<double> elem;
    std::vector<double> con1;
    std::vector<double> con2;
    std::vector<double> con3;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;



    while (fin >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k)
    {
        elem.push_back(a);
        con1.push_back(b);
        con2.push_back(c);
        con3.push_back(d);
        //std::cout << a << " " << b << " " << c << " " << d << endl;
    }

    ifstream finn1("C:\\Users\\omjag\\Downloads\\Cy_2_Corrdinates.txt");
    double cor1, cor2, cor3;

    while (finn1 >> cor1 >> cor2 >> cor3)
    {
        x.push_back(cor2);
        y.push_back(cor3);
        //std::cout << y << " " << z << " 0" << endl;
    }

    double npoin, nelem, ip1, ip2, ip3, loc1_x, loc1_y, loc2_x, loc2_y, loc3_x, loc3_y, D;

    npoin = x.size();
    nelem = elem.size();

    //Initailization of global stifness matrix
    vector<vector<double>> lhspo(npoin, vector<double>(npoin, 0));
    vector<double> D_mat(npoin, 0);
    vector<double> phi(npoin, 0);
    vector<double> phi_old(npoin, 0);
    vector<double> phi_new(npoin, 0);
    vector<double> diff(npoin, 0);

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
        loc1_x = x[con1[i] - 1];
        loc1_y = y[con1[i] - 1];

        loc2_x = x[con2[i] - 1];
        loc2_y = y[con2[i] - 1];

        loc3_x = x[con3[i] - 1];
        loc3_y = y[con3[i] - 1];

        D = loc1_x * (loc2_y - loc3_y) - loc1_y * (loc2_x - loc3_x) + (loc2_x * loc3_y - loc2_y * loc3_x);

        ip1 = ip1 - 1;
        ip2 = ip2 - 1;
        ip3 = ip3 - 1;


        //Element Matrix
        lhspo[ip1][ip1] = lhspo[ip1][ip1] + (((loc2_y - loc3_y) * (loc2_y - loc3_y) + (loc2_x - loc3_x) * (loc2_x - loc3_x)) / (2 * D));
        lhspo[ip1][ip2] = lhspo[ip1][ip2] + (((loc2_y - loc3_y) * (loc3_y - loc1_y) + (loc2_x - loc3_x) * (loc3_x - loc1_x)) / (2 * D));
        lhspo[ip1][ip3] = lhspo[ip1][ip3] + (((loc2_y - loc3_y) * (loc1_y - loc2_y) + (loc2_x - loc3_x) * (loc1_x - loc2_x)) / (2 * D));
        lhspo[ip2][ip1] = lhspo[ip2][ip1] + (((loc3_y - loc1_y) * (loc2_y - loc3_y) + (loc3_x - loc1_x) * (loc2_x - loc3_x)) / (2 * D));
        lhspo[ip2][ip2] = lhspo[ip2][ip2] + (((loc3_y - loc1_y) * (loc3_y - loc1_y) + (loc3_x - loc1_x) * (loc3_x - loc1_x)) / (2 * D));
        lhspo[ip2][ip3] = lhspo[ip2][ip3] + (((loc3_y - loc1_y) * (loc1_y - loc2_y) + (loc3_x - loc1_x) * (loc1_x - loc2_x)) / (2 * D));
        lhspo[ip3][ip1] = lhspo[ip3][ip1] + (((loc1_y - loc2_y) * (loc2_y - loc3_y) + (loc1_x - loc2_x) * (loc2_x - loc3_x)) / (2 * D));
        lhspo[ip3][ip2] = lhspo[ip3][ip2] + (((loc1_y - loc2_y) * (loc3_y - loc1_y) + (loc1_x - loc2_x) * (loc3_x - loc1_x)) / (2 * D));
        lhspo[ip3][ip3] = lhspo[ip3][ip3] + (((loc1_y - loc2_y) * (loc1_y - loc2_y) + (loc1_x - loc2_x) * (loc1_x - loc2_x)) / (2 * D));

    }

    // Saving Diagonal Matrix

    for (int i = 0; i < npoin; i++)
    {
        for (int j = 0; j < npoin; j++)
        {
            if (i == j)
            {
                D_mat[i] = lhspo[i][j];
                //cout << D_mat[i] << "\n";
            }
        }

    }



    ifstream finn2("C:\\Users\\omjag\\Downloads\\Cy_2_Boundary_Faces.txt");

    double aa, bb, cc, dd, ee, ff, nbface, len, xx, yy, m;
    std::vector<double> bface;
    std::vector<double> bcon1;
    std::vector<double> bcon2;
    std::vector<double> bc;
    std::vector<double> vb;



    while (finn2 >> aa >> bb >> cc >> dd >> ee >> ff)
    {
        bface.push_back(aa);
        bcon1.push_back(bb);
        bcon2.push_back(cc);
        bc.push_back(dd);
        vb.push_back(ee);
        //std::cout << y << " " << z << " 0" << endl;
    }

    nbface = bface.size();
    vector<double> rhspo(npoin, 0);
    vector<double> Z(npoin, 0);

    //cout << nbface;
    //cout << rhspo[0][0];

    for (int i = 0; i < nbface; i++)
    {

        //local nodes data
        loc1_x = x[bcon1[i] - 1];
        loc1_y = y[bcon1[i] - 1];

        loc2_x = x[bcon2[i] - 1];
        loc2_y = y[bcon2[i] - 1];

        len = sqrt(pow((loc2_y - loc1_y), 2) + pow((loc2_x - loc1_x), 2));

        xx = bcon1[i] - 1;
        yy = bcon2[i] - 1;

        m = atan2(-(loc2_x - loc1_x) , (loc2_y - loc1_y));

        //Element Matrix
        if (bc[i] == 4)
        {
            rhspo[xx] = rhspo[xx] + ((cos(m) * vb[i] * len) / 2);
            rhspo[yy] = rhspo[yy] + ((cos(m) * vb[i] * len) / 2);
            // rhspo[xx] = rhspo[xx] + ((((loc2_y - loc1_y) / abs(loc2_y - loc1_y)) * vb[i] * len) / 2);
           // rhspo[yy] = rhspo[yy] + ((((loc2_y - loc1_y) / abs(loc2_y - loc1_y)) * vb[i] * len) / 2);
            //cout << yy << " " << bcon2[i] << "\n";
            //cout << vb[i];
        }
    }


    // Initialize Phi values

    double  Vinf = 1;

    for (int i = 0; i < npoin; i++)
    {
        phi[i] = Vinf * x[i];
        //cout << phi_old[i] << "\n";
    }




    // Jacobian Iteration Method
    double Res;
    std::vector<double> res; //residual vector
    res.resize(npoin, 0); //init residual

/*       for (int i = 0; i < npoin; i++)
       {
           for (int j =0; j < npoin; j++)
           {
               Z[i] = Z[i] + lhspo[i][j] * phi[j];
           }
           cout << Z[i] << "\n";
       }
  */ 


 //   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    for (int k = 0; k < 2000; k++)
    {
        double sum = 0;
        //sum = 0;
        //cout << Z[i] << "\n";

        for (int j = 0; j < npoin; j++)
        {
            //cout << phi_old[j] << "\t" << phi[j] << "\n";\
            




            phi_old[j] = phi[j];
        //  Z[j] = std::inner_product(lhspo[j].begin(), lhspo[j].end(), phi.begin(), 0);
            Z[j] = 0;
            for (int i = 0; i < npoin; i++)
            {
                Z[j] = Z[j] + lhspo[j][i] * phi[i];
            }
            //cout << Z[j] << "\n";
            //cout << loc_var1 << "\n";

            
            phi[j] = ((rhspo[j] - Z[j]) / D_mat[j]) + phi[j];
            res[j] = phi[j] - phi_old[j];
            //diff[j] = pow((phi[j] - phi_old[j]), 2);
            //sum = sum + diff[j];
            // Res = sqrt(sum) / npoin;
            //cout << Res << "\n";
        }
        
        for (int i = 1; i < npoin; i++)
        {
            sum = sum + sqrt(pow(res[i], 2));
        }
        sum = sum / npoin;
        cout << sum << endl;
       

    }
    
    double a1, a2, a3, b1, b2, b3;
    vector<vector<double>> Vx_temp(npoin, vector<double>(nelem, 0));
    vector<vector<double>> Wt_temp(npoin, vector<double>(nelem, 0));
    vector<vector<double>> Vy_temp(npoin, vector<double>(nelem, 0));

    for (int i = 0; i < nelem; i++) 
    {
        ip1 = con1[i];
        ip2 = con2[i];
        ip3 = con3[i];
        //cout << i << " \t";

        //loc1_x = x[con1[i]-1];
        //cout << loc1_x << "\n";
        //local nodes data

        loc1_x = x[con1[i] - 1];
        loc1_y = y[con1[i] - 1];

        loc2_x = x[con2[i] - 1];
        loc2_y = y[con2[i] - 1];

        loc3_x = x[con3[i] - 1];
        loc3_y = y[con3[i] - 1];

        a1 = loc2_y - loc3_y;
        a2 = loc3_y - loc1_y;
        a3 = loc1_y - loc2_y;

        b1 = loc2_x - loc3_x;
        b2 = loc3_x - loc1_x;
        b3 = loc1_x - loc2_x;

        D = loc1_x * (loc2_y - loc3_y) - loc1_y * (loc2_x - loc3_x) + (loc2_x * loc3_y - loc2_y * loc3_x);

        ip1 = ip1 - 1;
        ip2 = ip2 - 1;
        ip3 = ip3 - 1;

        Vx_temp[ip1][i] = (phi[ip1] * (a1 / D)) + (phi[ip2] * (a2 / D)) + (phi[ip3] * (a3 / D));
        Vx_temp[ip2][i] = Vx_temp[ip1][i];
        Vx_temp[ip3][i] = Vx_temp[ip1][i];

        Vy_temp[ip1][i] = (phi[ip1] * (b1 / D)) + (phi[ip2] * (b2 / D)) + (phi[ip3] * (b3 / D));
        Vy_temp[ip2][i] = Vy_temp[ip1][i];
        Vy_temp[ip3][i] = Vy_temp[ip1][i];


        Wt_temp[ip1][i] = D / 2;
        Wt_temp[ip2][i] = D / 2;
        Wt_temp[ip3][i] = D / 2;

    }

    vector<double> Vx(npoin, 0);
    vector<double> Vy(npoin, 0);
    vector<double> Vz(npoin, 0);

    for(int i =0; i<npoin;i++)
    {
        double loc_var1 = 0;
        double loc_var2 = 0;
        double loc_var3 = 0;

        for (int j = 0; j< nelem ; j++)
        {
            loc_var1 = loc_var1 + (Vx_temp[i][j] * Wt_temp[i][j]);
            loc_var2 = loc_var2 + (Vy_temp[i][j] * Wt_temp[i][j]);
            loc_var3 = loc_var3 + Wt_temp[i][j];
        }

        Vx[i] = loc_var1 / loc_var3;
        Vy[i] = loc_var2 / loc_var3;

    }

    vector<double> Ur(npoin, 0);
    vector<double> Ut(npoin, 0);
    vector<double> Vr(npoin, 0);
    vector<double> Vt(npoin, 0);
    vector<double> Vmag(npoin, 0);
    vector<double> Umag(npoin, 0);
    vector<double> Er(npoin, 0);
    double theta, r, sum = 0, Err;

    for (int i = 0; i < npoin; i++)
    {
        r = sqrt(pow(x[i], 2)+pow(y[i], 2));
        theta = atan2((y[i]), x[i]);

        Ur[i] = cos(theta) * (1 - (0.7071 / (pow(r, 2))));
        Ut[i] = -sin(theta) * (1 + (0.7071 / (pow(r, 2))));
        Vr[i] = Vx[i] * cos(theta) + Vy[i] * sin(theta);
        Vt[i] = Vx[i] * sin(theta) - Vy[i] * cos(theta);

        Vmag[i] = sqrt(pow(Vr[i], 2) + pow(Vt[i], 2));
        Umag[i] = sqrt(pow(Ur[i], 2) + pow(Ut[i], 2));
        Er[i] = pow((Vmag[i] - Umag[i]), 2);
        sum = sum + Er[i];
    }
    
    Err = sqrt(sum)/npoin;
    cout << Err;

/*

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

*/

    // 



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
            fw1 << rhspo[i] << "\t";
            fw1 << "\n";
        }
        fw1.close();
    }

    else cout << "Problem with opening file";


    /// Printing out Phi
    ofstream
        fw2("C:\\Users\\omjag\\Downloads\\phi.txt");
    std::ofstream::out;

    if (fw2.is_open())
    {
        for (int i = 0; i < npoin; i++)
        {
            fw2 << phi[i];
            fw2 << "\n";
        }
        fw2.close();
    }

    else cout << "Problem with opening file";

    ofstream
        fw3("C:\\Users\\omjag\\Downloads\\Vx_Vy.txt");
    std::ofstream::out;

    if (fw3.is_open())
    {
        for (int i = 0; i < npoin; i++)
        {
            fw3 << Vx[i] << " " << Vy[i] << " " << Vz[i];
            fw3 << "\n";
        }
        fw3.close();
    }

    else cout << "Problem with opening file";

    ofstream
        fw4("C:\\Users\\omjag\\Downloads\\x_Vx_Vy.dat");
    std::ofstream::out;
    
    if (fw4.is_open())
    {
        for (int i = 0; i < nbface; i++)
        {
            fw4 << x[bcon1[i] - 1] << " " << Vx[bcon1[i] - 1] << " " << Vy[bcon1[i] - 1];
            fw4 << "\n";
        }
        fw4.close();
    }

    else cout << "Problem with opening file";

    ofstream
        fw6("C:\\Users\\omjag\\Downloads\\Cy_x_Vx_Vy.dat");
    std::ofstream::out;

    if (fw6.is_open())
    {
        for (int i = 0; i < nbface; i++)
        {
            if (bc[i] == 2)
            {
                fw6 << x[bcon1[i] - 1] << " " << sqrt(pow(Vx[bcon1[i] - 1],2) + pow(Vy[bcon1[i] - 1],2)) << "\n";
                fw6 << "\n";
            }
            
        }
        fw6.close();
    }

    else cout << "Problem with opening file";


    ofstream
        fw5("C:\\Users\\omjag\\Downloads\\Cy_1_VTK.vtk");
    std::ofstream::out;

    if (fw5.is_open())
    {
        fw5 << "# vtk DataFile Version 2.0" << "\n";
        fw5 << "Unstructured Grid" << "\n";
        fw5 << "ASCII" << "\n";
        fw5 << "DATASET UNSTRUCTURED_GRID" << "\n";
        fw5 << "\n";
        fw5 << "POINTS " << npoin << " float";
        fw5 << "\n";
        
        for (int i = 0; i < npoin; i++)
        {
            fw5 << x[i] << " " << y[i] << " 0";
            fw5 << "\n";
        }

        fw5 << "\n";
        fw5 << "CELLS " << nelem << " " << (nelem * 4) << "\n";
        
        for (int i = 0; i < nelem; i++)
        {
            fw5 << 3 << " " << (con1[i] - 1) << " " << (con2[i] - 1) << " " << (con3[i] - 1) << "\n";
        }

        fw5 << "\n";
        fw5 << "CELL_TYPES " << nelem << "\n";

        for (int i = 0; i < nelem; i++)
        {
            fw5 << "5" << "\n";
        }

        fw5 << "\n";
        fw5 << "POINT_DATA " << npoin << "\n";
        fw5 << "SCALARS Velocity_Potential double 1" << "\n";
        fw5 << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < npoin; i++)
        {
            fw5 << phi[i] << "\n";
        }

        fw5 << "\n";
        fw5 << "FIELD FieldData 1" << "\n";
        fw5 << "Velocity 3 " << npoin << " double" << "\n";

        for (int i = 0; i < npoin; i++)
        {
            fw5 << Vx[i] << " " << Vy[i] << " 0" << "\n";
        }


        fw5.close();
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
