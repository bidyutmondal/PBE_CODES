#include <bits/stdc++.h>

using namespace std;

// define the class size;
int class_size = 30;

vector<double> x_arr()
{
    double x_p_min = pow(10, -6);
    double x_p_max = pow(10, 3);
    vector<double> x_p(class_size+1);
    x_p[0] = x_p_min;
    x_p[class_size] = x_p_max;
    for(int i=1;i<class_size;i++){
        double alpha = log(x_p_min)+i*((log(x_p_max)-log(x_p_min))/class_size);
        x_p[i] = exp(alpha);
    }
    vector<double> x(class_size);
    for(int i=0;i<class_size;i++){
        x[i] = (x_p[i]+x_p[i+1])/2;
    }
    return x;
}
// f = e^-x
// write generalised form of f and beta

vector<double> f_arr(vector<double> &x)
{
    vector<double> f(class_size);
    for (int i = 0; i < class_size; i++)
    {
        f[i] = exp(-1 * (i + 1) * x[i]);
    }
    return f;
}

vector<vector<int>> beta_arr(vector<double> &x)
{
    int max_vol = 100; // define max_size
    vector<vector<int>> beta(class_size, vector<int>(class_size));
    for (int i = 0; i < class_size; i++)
    {
        for (int j = 0; j < class_size; j++)
        {
            if (x[i] + x[j] <= max_vol)
                beta[i][j] = 1;
            else
                beta[i][j] = 0;
        }
    }
    return beta;
}

double weight_const_b(vector<double> &x, int i, int j, int k)
{
    if (x[i] > x[j] + x[k])
        return 0.0;
    double xlj = 0;
    for (auto e : x)
    {
        if (e > x[j] + x[k])
        {
            xlj = e;
            break;
        }
    }
    return (x[j] + x[k]) / (2 * xlj - (x[j] + x[k]));
}

double weight_const_d(vector<double> &x, int i, int j, int k)
{
    if (x[i] > x[j] + x[k])
        return 0.0;
    double xlj = 0;
    for (auto e : x)
    {
        if (e > x[j] + x[k])
        {
            xlj = e;
            break;
        }
    }
    return (xlj) / (xlj - (x[j] + x[k]));
}

int main()
{
    vector<double> x = x_arr();
    cout << "\nX\n";
    for (auto e : x)
    {
        cout << e << " ";
    }
    cout << endl
         << "F \n";
    vector<double> f = f_arr(x);
    for (auto e : f)
    {
        cout << e << " ";
    }
    cout << endl;
    // write a generalised form of beta
    vector<vector<int>> beta = beta_arr(x); // beta taken such that it follow the NFVS modification
    vector<vector<double>> f_time;
    f_time.push_back(f);

    // f = f + del(t)*pbe
    int t_final = 100;
    for (double n = 1; n < t_final; n += 1) // n is taken as time in equation
    {
        for (int i = 1; i < class_size; i++) // i is cell no
        {
            double birth_term = 0.0;
            for (int j = 1; j < i; j++)
            {
                for (int k = 1; k < i; k++)
                {
                    if (x[i-1]<x[j] + x[k] < x[i])
                        birth_term += 1 * f_time[n-1][j] * f_time[n-1][k] * (((x[j] - x[j - 1]) * (x[k] - x[k - 1])) / (x[i] - x[i - 1]));
                }
            }
            double death_term = 0.0;
            for (int j = 1; j < class_size; j++) // upto capital I
            {
                death_term += beta[i][j] * f_time[n-1][i] * f_time[n-1][j] * (x[j] - x[j - 1]); // class size var can be changed to something
            }
            // f = f + del(t)(1.0/2*birth_term-death_term)
            f[i] = f_time[n-1][i] + 0.01 *(0.5 *0.00002*birth_term - death_term); // 0.5 to ommit symmetry
        }
        f_time.push_back(f);
    }

    cout << "number of particles:\n";
    for (int i = 0; i < f_time[99].size(); i++)
    {
        cout << "classSize --> f : "<<x[i]<<" --> "<<f_time[99][i]<<"\n";
    }
    cout << endl;

    //writing in file
    ofstream numberData;
    numberData.open("number_of_particles_data_sfvm.txt");
    for (int i = 0; i < f_time[99].size()-10; i++)
    {
        numberData<<x[i]<<" "<<f_time[99][i]*x[i]<<"\n";
    }
    numberData.close();
    return 0;
}
