// Hw3.cpp : This file contains the 'main' function. Program execution begins and ends there.
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
    double PL, VL, rho_L, PR, VR, rho_R, P, P2, CR, V2, gamma, alpha, rho_2, P3, V3, rho_3, C3, t;

    t = 0.2;
    gamma = 1.4;

    PL = 1;
    rho_L = 1;
    VL = 0;
    PR = 0.1;
    VR = 0;
    rho_R = 0.125;

    P = (((PR + PL) / 2) / PR);
    cout << "P " << P << "\n";

    double A, B, Num3, Den3, Num4, Den4, f_x, f_x1, aa, bb, pw, CL;
    CL = sqrt(gamma * (PL / rho_L));
    CR = sqrt(gamma * (PR / rho_R));

    double a1, a2;

    for (int i = 0; i < 5 ;i++) 
    {
        a2 = sqrt(((gamma + 1) / 2) * (P - 1)) + 1;
        a1 = 1 - ((gamma - 1) * (VR - VL + (CR * (P - 1)) / (gamma * a2)) / (2 * CL));

        A = pow(a1, -(2 * gamma) / (gamma - 1));
        Num3 = CR * ((gamma + 1) / 2) * (P - 1);
        Den3 = 2 * gamma * sqrt(((gamma + 1) / 2) * (P - 1)) * pow(a2, 2);
        B = P * gamma * ((CR / (gamma * a2)) - (Num3 / Den3)) / (CL * pow(a1, ((2 * gamma / (gamma - 1)) + 1)));

        f_x1 = A + B;

        aa = (gamma - 1) / (2 * CL);
        Num4 = P - 1;
        Den4 = sqrt(((gamma + 1) / 2) * ((P)-1)) + 1;
        bb = VL - VR - ((CR / gamma) * (Num4 / Den4));
        pw = -2 * gamma / (gamma - 1);
        f_x = (P * pow((1 + aa * (bb)), (pw))) - (PL / PR);

        P = P - (f_x / f_x1);

        cout << "P " << P << "\n";
        cout << "f_x " << f_x << "\n";
        cout << "f_x1 " << f_x1 << "\n";
    }
    
    
    


    // Region 2
    // Pressure 2
    P2 = P * PR;
    cout << "P2 " << P2 << "\n";
    // Vel 2
    double Num1, Den1;
    Num1 = P - 1;
    Den1 = sqrt((((gamma + 1) / 2) * (P - 1))) + 1;
    
    V2 = VR + (CR / gamma) * (Num1 / Den1);
    cout << "V2 " << V2 << "\n";

    //rho 2
    alpha = (gamma + 1) / (gamma - 1);
    rho_2 = rho_R * ((1 + (alpha * P))/(alpha + P));
    cout << "rho_2 " << rho_2 << "\n";
    
    // Region 3
    P3 = P2;
    V3 = V2;
    cout << "P3 " << P3 << "\n";
    cout << "V3 " << V3 << "\n";
    rho_3 = pow(((PL / (pow(rho_L,gamma))) * P3), (1 / gamma));
    C3 = sqrt(gamma * (P3 / rho_3));
    cout << "rho_3 " << rho_3 << "\n";
    cout << "C3 " << C3 << "\n";

    // Region 4
    //lambda 0
    double xh, xt,rho_4, P4, V4, C4, xc, xs, xl, xr, sh = 0.5 ;

    //CL = sqrt(gamma * (PL / rho_L));
    xh = (VL - CL) * t;
    xt = (V3 - C3) * t;
    cout << "CL " << CL << "\n";
    
    V4 = (2 / (gamma + 1)) * ((xh / t) + CL + ((gamma - 1) * VL / 2));
    cout << "V4 " << V4 << "\n";

    C4 = V4 - (xh / t);
    cout << "C4 " << C4 << "\n";

    P4 = PL * pow((C4 / CL), ((2 * gamma) / (gamma - 1)));
    cout << "P4 " << P4 << "\n";
    rho_4 = pow(((PL / (pow(rho_L, gamma))) * P4), (1 / gamma));
    cout << "\n" << "rho_4 " << rho_4 << "\n";

    xc = V2 * t;
    xs = (((rho_R * VR) - (rho_2 * V2)) / (rho_R - rho_2))*t;
    xl = 0;
    xr = 1;

    cout << "\n" << "xh " << xh << "\n";
    cout << "\n" << "xt " << xt << "\n";
    cout << "\n" << "xc " << xc << "\n";
    cout << "\n" << "xs " << xs << "\n";

    ofstream
        fw("C:\\Users\\omjag\\Downloads\\H3_Vel.dat");
    std::ofstream::out;

    if (fw.is_open())
    {
        fw << xl  << "\t" << VL << "\n";
        fw << xh + sh << "\t" << V4 << "\n";
        fw << xt + sh << "\t" << V3 << "\n";
        fw << xc + sh << "\t" << V3 << "\n";
        fw << xs + sh << "\t" << V2 << "\n";
        fw << xr  << "\t" << VR << "\n";
        fw.close();
    }

    else cout << "Problem with opening file";

    ofstream
        fw1("C:\\Users\\omjag\\Downloads\\H3_P.dat");
    std::ofstream::out;

    if (fw1.is_open())
    {
        fw1 << xl  << "\t" << PL << "\n";
        fw1 << xh + sh << "\t" << P4 << "\n";
        fw1 << xt + sh << "\t" << P3 << "\n";
        fw1 << xc + sh << "\t" << P3 << "\n";
        fw1 << xs + sh << "\t" << P2 << "\n";
        fw1 << xr  << "\t" << PR << "\n";
        fw1.close();
    }

    else cout << "Problem with opening file";

    ofstream
        fw2("C:\\Users\\omjag\\Downloads\\H3_rho.dat");
    std::ofstream::out;

    if (fw2.is_open())
    {
        fw2 << xl  << "\t" << rho_L << "\n";
        fw2 << xh + sh << "\t" << rho_4 << "\n";
        fw2 << xt + sh << "\t" << rho_3 << "\n";
        fw2 << xc + sh << "\t" << rho_3 << "\n";
        fw2 << xc + sh << "\t" << rho_2 << "\n";
        fw2 << xs + sh << "\t" << rho_2 << "\n";
        fw2 << xs + sh << "\t" << rho_R << "\n";
        fw2 << xr  << "\t" << rho_R << "\n";
        fw2.close();
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
