#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <algorithm>
#include <thread>

#define _USE_MATH_DEFINES
#define GNUPLOT "gnuplot -persist"
#define NX 720
 
void add_m_m(double*, const double*, const double*);
void sum_m_m(double*, const double*, const double*);
void mul_m_v(double*, const double*, const double*);
void add_v_v(double*, const double*, const double*);
void mul_s_v(double*, const double, const double*);
void mul_s_m(double*, const double, const double*);
void inv_m(double*, const double*);

const double tl = 2000; // 線路長 track length
const double max_vel = 40; // 最高速度 maximum velocity[km/h]
const double sim_time = 50;   // simulation time
const double t_i = 0;
const double t_f = sim_time;
const double v_i = 0;
const double v_f = 0;
const double x_i = 0;
const double x_f = tl;
const double weight = 353; // [ton]
// train's resistances correspond to train's speed [km/h]
const double a = 14.2974;
const double b = 0.1949;
const double c = 0.8095;
// Track gradient
const double p1 = 0;
const double p2 = 25;
const double p3 = 0;
const double z1 = 200;
const double z2 = 350;
const double z3 = 800;
const double g = 9.8; // gravitational constant [m/s^2]
const double l = 30;  // train's length[m]
std::vector<double> u_notch{-1, -0.9, -0.5, -0.1, 0, 0.1, 0.5, 0.9, 1};
const int notch_count = u_notch.size();
const double pen1 = std::pow(10, 6);
const double pen2 = std::pow(10, 100);

// 全部2次元
void add_m_m(double *m_ans, const double *m1, const double *m2)
{
    m_ans[0] = m1[0] + m2[0];
    m_ans[1] = m1[1] + m2[1];
    m_ans[2] = m1[2] + m2[2];
    m_ans[3] = m1[3] + m2[3];
}

void sub_m_m(double *m_ans, const double *m1, const double *m2)
{
    m_ans[0] = m1[0] - m2[0];
    m_ans[1] = m1[1] - m2[1];
    m_ans[2] = m1[2] - m2[2];
    m_ans[3] = m1[3] - m2[3];
}

void mul_m_v(double *v_ans, const double *m, const double *v)
{
    v_ans[0] = m[0] * v[0] + m[1] * v[1];
    v_ans[1] = m[2] * v[0] + m[3] * v[1];
}

void add_v_v(double *v_ans, const double *v1, const double *v2)
{
    v_ans[0] = v1[0] + v2[0];
    v_ans[1] = v1[1] + v2[1];
}
// スカラー倍
void mul_s_v(double *v_ans, const double s, const double *v)
{
    v_ans[0] = s * v[0];
    v_ans[1] = s * v[1];
}

void mul_s_m(double *m_ans, const double s, const double *m)
{
    m_ans[0] = s * m[0];
    m_ans[1] = s * m[1];
    m_ans[2] = s * m[2];
    m_ans[3] = s * m[3];
}
// 逆行列
void inv_m(double *m_ans, const double *m)
{
    // m_ans == mのときは計算がうまくいかなくなるから別の変数を用いる。
    double a = m[0];
    double b = m[1];
    double c = m[2];
    double d = m[3];
    double det = a * d - b * c;

    if (det == 0) {
        printf("det = 0\n");
        exit(1);
    }

    m_ans[0] = d / det;
    m_ans[1] = -b / det;
    m_ans[2] = -c / det;
    m_ans[3] = a / det;
}

double ms2kmh(double v)
{
    return 3.6 * v;
}

double kmh2ms(double v)
{
    return v / 3.6;
}

void solve_next_state(double *next_xv, const double *curr_xv, const double *A, const double *B, const double dt, const double F0)
{
    double identity[4] = {1, 0, 0, 1};
    double C[4];
    double D[4];
    double E[2];
    double G[2];
    double H[2];

    // 式: next_xv = (I - dt * A / 2).inverse() * (((I + dt * A / 2) * curr_xv) + (dt * B * F0));
    mul_s_m(C, 0.5 * dt, A);  // C = 0.5 * dt * A
    add_m_m(D, identity, C);
    sub_m_m(C, identity, C);    // C = I - C
    inv_m(C, C);

    mul_m_v(E, D, curr_xv);
    mul_s_v(G, F0 * dt, B);
    add_v_v(H, E, G);
    mul_m_v(next_xv, C, H);
}

void traction(double *F, double *dFdv, const double u, const double v)
{
    double f_a = 3.3 * 1000 / 3.6; // [N/t]
    double v_a1 = 40 / 3.6;    // [km/3.5*1000/3.6h]->[m/s]
    double v_a2 = 80 / 3.6;   // [km/3.5*1000/3.6h]->[m/s]

    // regenative
    double f_b = 4 * 1000 / 3.6;    // [N/t]
    double v_b = 80 / 3.6;    // [km/h]->[m/s]

    *F = 0;
    *dFdv = 0;


    if (u > 0) {
        if (v < v_a1) {
            *F = f_a * u;
            *dFdv = 0;
        }
        else if (v < v_a2) {
            *F = f_a * u * v_a1 /v;
            *dFdv = -f_a * u * v_a1 / std::pow(v, 2);
        }
        else {
            *F = f_a * u * v_a1 * v_a2 / std::pow(v, 2);
            *dFdv = -2 * f_a * u * v_a1 * v_a2 / std::pow(v, 3);
        }
    }
    else if (u < 0) {
        if (v < v_b) {
            *F = f_b * u;
            *dFdv = 0;
        }
        else {
            *F = f_b * u * std::pow(v_b / v, 2);
            *dFdv = -2 * f_b * u * std::pow(v_b, 2) / std::pow(v, 3);
        }
    }
    else {
        *F = 0;
        *dFdv = 0;
    }

    //printf("traction %f %f %f %f %f %f %f %f %f\n", u, v, f_a, v_a1, v_a2, f_b, v_b, F, dFdv);
}

void difsolve(double *next_x, double *next_v, double *dW, double curr_x, double curr_v, double dt, const double u)
{
    curr_v = kmh2ms(curr_v);
    double F;
    double dFdv;
    // calculate traction force
    traction(&F, &dFdv, u, curr_v);
    // calculation for train's resistance
    double R;
    double dRdv;
    if (curr_v > 0) {
        R = a + (b * 3.6 * curr_v) + ((c * std::pow(3.6 * curr_v, 2)) / weight);
        dRdv = (2 * c * std::pow(3.6, 2) * curr_v / weight) + 3.6 * b;
    }
    else {
        R = 30; // starting resistance
        dRdv = 0;
    }
    // [N/ton]->[N]
    F *= weight;
    dFdv *= weight;  // [N/(m/s)]
    R *= weight ; // [N/ton]->[N]
    dRdv *= weight;  // [N(m/s)]
    //printf("difsolve%f %f %f %f\n", F, dFdv, R, dRdv);

    // calculating force used to accelerate / deccelerate the train
    // Running resistance for straight line
    // F0 = F - (dFdv*curr_v) - R + (dRdv*curr_v); //[N]
    // beta=0;

    //traction affected by the gradient
    double R12;
    double dR12dx;
    double R23;
    double dR23dx;
    double F0;
    double alpha;
    double beta;
    if ((z1 <= curr_x) && (curr_x < z2)) {
        if (curr_x < z1 + l) { //affected by the former slope
           R12 = weight * g * (p1 * (z1 - (curr_x - l)) + p2 * (curr_x - z1)) / l;//[N]
           dR12dx = weight * g * (-p1 + p2) / l;//[N/m]
        }
        else {
            R12 = weight * g * p2;
            dR12dx = 0;
        }
        F0 = F - (dFdv * curr_v) - R + (dRdv * curr_v) - R12 + (dR12dx * curr_x);
        beta = -dR12dx;
    }
    else if ((z2 <= curr_x) && (curr_x < z3)) {
        if (curr_x < z2+l) { //affected by the former slope
            R23 = weight * g * (p2 * (z2 - (curr_x - l)) + p3 * (curr_x - z2)) / l;//[N]
            dR23dx = weight * g * (-p2 + p3)/l;//[N/m]
        }
        else {
            R23 = weight * g * p3;
            dR23dx = 0;
        }       
        F0 = F - (dFdv * curr_v) - R + (dRdv * curr_v) - R23 + (dR23dx * curr_x);
        beta = -dR23dx;
    }
    else {
        F0 = F - (dFdv * curr_v) - R + (dRdv * curr_v); //[N] traction for straight line
        beta = 0;
    }

    F0 /= (weight * 1000); //[N]->[N/ton]->[N/kg]=[m/s/s]

    alpha = dFdv - dRdv; //[N/(m/s)]
    alpha /= (1000 * weight); // //[N/(m/s)] -> [N/ton/(m/s)] -> [N/kg/(m/s)]=[(m/s/s)/(m/s)]
    beta = beta / (1000 * weight); // [N/(m/s)] -> [N/ton/(m/s)] -> [N/kg/(m/s)]=[(m/s/s)/(m/s)]
/*
    Eigen::Matrix2d A;
    A(0, 0) = 0;
    A(1, 0) = beta;
    A(0, 1) = 1;
    A(1, 1) = alpha;
    Eigen::Vector2d B;
    B(0) = 0;
    B(1) = 1;
    Eigen::Vector2d curr_xv;
    curr_xv(0) = curr_x,
    curr_xv(1) = curr_v;
    Eigen::Vector2d next_xv;
    Eigen::Matrix2d I = Eigen::MatrixXd::Identity(2, 2);
    I(0, 0) = 1;
    I(1, 0) = 0;
    I(0, 1) = 0;
    I(1, 1) = 1;
    */

    double A[4] = {0, 1, beta, alpha};
    double B[2] = {0, 1};
    double curr_xv[2] = {curr_x, curr_v};
    double next_xv[2];

    //solving the next state
    solve_next_state(next_xv, curr_xv, A, B, dt, F0);

    //printf("difsolve%f %f %f %f\n", next_xv(0), next_xv(1), curr_xv(0), curr_xv(1));
    if (next_xv[1] < 0) {
        dt *= (curr_xv[1] / (curr_xv[1] - next_xv[1]));
        solve_next_state(next_xv, curr_xv, A, B, dt, F0);
    }
/*
    //solving the next state
    next_xv = (I - dt * A / 2).inverse() * (((I + dt * A / 2) * curr_xv) + (dt * B * F0));

    //printf("difsolve%f %f %f %f\n", next_xv(0), next_xv(1), curr_xv(0), curr_xv(1));
    if (next_xv(1) < 0) {
        dt *= (curr_xv(1) / (curr_xv(1) - next_xv(1)));
        next_xv = (I - dt * A / 2).inverse() * (((I + dt * A / 2) * curr_xv) + (dt * B * F0));
    }
*/
    double dx = next_xv[0] - curr_x;
    double dv = next_xv[1] - curr_xv[1];
    double vn = next_xv[1]; //[m/s]
    double Fn;
    double dev_null;
    traction(&Fn, &dev_null, u, curr_v); //[N]
    Fn *= weight;

/*
    double a = (dFdv * dv) / dt;
    double b = F;
    double c = dv / dt;
    double d = curr_v;
    */

/*
    Eigen::Matrix2d tmptmp;
    tmptmp << 1, 0, 
              0, 3.6;
    next_xv = tmptmp * next_xv;
    */
    next_xv[1] = ms2kmh(next_xv[1]);
    // efficiency
    double eta_a = 0.9;
    double eta_b = 0.8;
    //

    ////energy calculation
    if (u > 0) {
        *dW  = 0.5 * (F * curr_v + Fn * vn ) * dt / eta_a / 3600 * 0.001 * 1.1;  //[m/s]*[N]*[s]=[J]=1/(3.6*10^6)[kWh]
        //*dW  = 0.5*(F*curr_v + (F+dFdv*dv)*vn ) * dt / eta / 3600 / 1000 *1.1;  
        //*dW  = (F + 0.5*dFdv*dv) * dx / eta / 3600 / 1000;  
        //*dW  =(a*c)/3*dt^3+(a*d+b*c)/2*dt^2+(b*d)*dt / eta / 3600 / 1000;  
    }
    else if (u == 0) {
        *dW  = 0;
    }
    else {
        *dW  =  0.5 * (F * curr_v + Fn * vn ) * dt * eta_b / 3600 * 0.001 * 1.1;
        //*dW  = 0.5*(F*curr_v + (F+dFdv*dv)*vn ) * dt * eta / 3600 / 1000 *1.1 ;   
        //*dW  = (F + 0.5*dFdv*dv) * dx * eta / 3600 / 1000;     
        //*dW  = (a*c)/3*dt^3+(a*d+b*c)/2*dt^2+(b*d)*dt * eta / 3600 / 1000;     
    }

//printf("%f %f %f\n", *dW, next_xv(0), next_xv(1));
    *next_x = next_xv[0];
    *next_v = next_xv[1];
}

double val2latno(double n, const std::vector<double> *latt, int latt_size)
{
    n = std::max(std::min(n, (*latt)[latt_size - 1]), (*latt)[0]);

//std::cout << static_cast<double>(latt.size()) << std::endl;
/*
// section?
    int section_start = 0;
    int section_end = latt_size - 1;

        //printf("%f sect0 %d sect1 %d\n", n, section_start, section_end);
    int i;
    while (section_end - section_start > 1) {
        i = section_start + ((section_end - section_start) >> 1);
        //
        // sect(0)-----i-----sect(1)
        //
        if (latt[i] <= n) {
            section_start = i;
        }
        else {
            section_end = i;
        }
        //printf("%d %d %d %f\n", i, section_start, section_end, latt[i]);
    }
//std::cout << "n "<< n << " i " << i;

    // wakaranai
    int tmp = std::min(i + 3, latt_size);
    for (int j = std::max(0, i - 3); j < tmp; j++) {
        //puts("b");
        if (latt[j] <= n) {
            i = j;
        }
    }

//std::cout << " return " << i << std::endl;
*/
    for (int i = 0; i < latt_size - 1; i++) {
        if ((*latt)[i] <= n && n < (*latt)[i + 1]) {
            return i;
        }
    }
    //std::cout << latt_size << " " << latt[i] << " " << n << " " << i << std::endl;
    return latt_size - 1;
}

double interpolate(const std::vector<std::vector<std::vector<double>>> *label, int time, double next_x, double next_v, const std::vector<double> *x_latt, const std::vector<double> *v_latt, const int x_latt_size, const int v_latt_size)
{
    int k = val2latno(next_x, x_latt, x_latt_size);
    //printf(" ");
    int j = val2latno(next_v, v_latt, v_latt_size);
if (time == label->size() - 4 + 1) {
//printf("k=%d j=%d time=%d next_x=%f, next_v=%f %f %f\n", k, j, time, next_x, next_v, x_latt[k], v_latt[j]);
}
    /*
     |      |          |
    y2---(k,j+1)---(k+1,j+1)
     |      |          |
     |      |          |
     |      |          |
    y1----(k,j)-----(k+1,j)
     |      |          |
         --x1---------x2--
         */
    /*
    puts("");
    for (double u : label[0][0]) {
        printf("%f ", u);
    }
    puts("");
    for (double u : label[0][1]) {
        printf("%f ", u);
    }
    puts("");
*/
//printf("label size %lu %lu\n", label.size(), label[0].size());
    double u_11 = (*label)[time][k][j];   // (k,j)
    double u_12 = (*label)[time][std::min(k + 1, x_latt_size - 1)][j]; // (k+1,j)
    double u_21 = (*label)[time][k][std::min(j + 1, v_latt_size - 1)]; // (k,j+1)
    double u_22 = (*label)[time][std::min(k + 1, x_latt_size - 1)][std::min(j + 1, v_latt_size - 1)];   // (k+1,j+1)
        //puts("d");

    double g1;
    double g2;

    if (k == x_latt->size() - 1) {
        g1 = u_11;  // u_12 = u_11
        g2 = u_21;  // u_22 = u_21
    }
    else {
        g1 = ((u_12 * (next_x - (*x_latt)[k])) + (u_11 * ((*x_latt)[std::min(k + 1, x_latt_size - 1)] - next_x))) / ((*x_latt)[std::min(k + 1, x_latt_size - 1)] - (*x_latt)[k]);
        g2 = ((u_22 * (next_x - (*x_latt)[k])) + (u_21 * ((*x_latt)[std::min(k + 1, x_latt_size - 1)] - next_x))) / ((*x_latt)[std::min(k + 1, x_latt_size - 1)] - (*x_latt)[k]);
    }

    double u;
//if (time == label.size() - 4 + 1) {
////printf("u11=%f u12=%f u21=%f u22=%f g1=%f g2=%f\n", u_11, u_12, u_21, u_22, g1, g2);
//}

    if (j == v_latt->size() - 1) {
        u = g1; // g1 = g2
    }
    else {
        u = ((g2 * (next_v - (*v_latt)[j])) + (g1 * ((*v_latt)[std::min(j + 1, v_latt_size - 1)] - next_v))) / ((*v_latt)[std::min(j + 1, v_latt_size - 1)] - (*v_latt)[j]);
    }
//printf("%f\n", u);
    return u;
}

void calc(int i, int j, int dt, const std::vector<double> *x_latt, const std::vector<double> *v_latt, const std::vector<double> *t_latt, std::vector<std::vector<std::vector<double>>> *e_opt, std::vector<std::vector<std::vector<double>>> *u_opt)
{
    for (int k = 0; k < v_latt->size(); k++) {
        double e_opt0 = pen1;
        double u_opt0 = 0;
        double e_opt1 = 0;
        for (int u = 0; u < notch_count; u++) {
            double next_x;
            double next_v;
            double dW;
            difsolve(&next_x, &next_v, &dW, (*x_latt)[j], (*v_latt)[k], dt, u_notch[u]);
            //if (t_size - 2 - i < 3)
            //printf("notch %f, x %f v %f dW %f\n", u_notch[u], next_x, next_v, dW);
            // speed limitation線路の形状による速度制限
            if (((*x_latt)[j] < 100) && ((*v_latt)[k] >= 60)) {
                e_opt1 = pen1;   
            }
            else if (((*x_latt)[j] >= 100 && (*x_latt)[j] < 450 && (*v_latt)[k] >= 40)) {
                e_opt1 = pen1;   
            }
            else if (((*x_latt)[j] >= 450 && (*x_latt)[j] < 550 && (*v_latt)[k] >= 25)) {
                e_opt1 = pen1;   
            }
            else if (((*x_latt)[j] >= 550 && (*x_latt)[j] < 800 && (*v_latt)[k] >= 40)) {
                e_opt1 = pen1;   
            }
            /*  
            else if ((x_latt[j] > x_latt[x_latt.size() - 1])) {
                e_opt1 = pen2;   
            }*/
            else {
                e_opt1 = dW + interpolate(e_opt, i + 1, next_x, next_v, x_latt, v_latt, x_latt->size(), v_latt->size()); // linear interpolation by interpol program
            }
            /*
            if (i==t_size-4) {
            //printf("notch %f, x %f v %f dW %f\n", u_notch[u], next_x, next_v, dW);
                //std::cout << e_opt1 << std::endl;
                
    for (double u : u_opt[t_size - 4][0]) {
    printf("%f\n ", u);
    }
    return 0;
                //std::cout << u +1 << " " << dW << " " << e_opt0 << " " << u_notch[u] << std::endl;
            }*/
            //e_opt1 = dW + interpol(e_par,X_n,x_latt,v_latt,x_nola,v_nola); %linear interpolation by interpol program
            if (e_opt1 < e_opt0) { // finding the optimum notch by evaluating the energy 最小を記憶
                e_opt0 = e_opt1;
                u_opt0 = u_notch[u];
            }
        }
        // 消費エネルギーが最小になるノッチと、そのときの消費エネルギーが求まった。
        (*e_opt)[i][j][k] = e_opt0;
        (*u_opt)[i][j][k] = u_opt0;
        //printf("%f %f\n", e_opt0, u_opt0);
    }
}

void assign(int time, int n, int x_latt_size, int nproc, int dt, const std::vector<double> *x_latt, const std::vector<double> *v_latt, const std::vector<double> *t_latt, std::vector<std::vector<std::vector<double>>> *e_opt, std::vector<std::vector<std::vector<double>>> *u_opt)
{
    for (int j = n; j < x_latt_size; j += nproc) {
        calc(time, j, dt, x_latt, v_latt, t_latt, e_opt, u_opt);
    }
}

const double penalty_v = 1;
const double penalty_x = 1;

int main()
{
    // 位置、速度、時間の格子 lattice for position, speed and time
    //int tmp = ((tl - 200) / 5 + 1) + 190 + 18 + 20 + 1; // ((tl-10)-(tl-199))/1+1=190
    std::vector<double> x_latt;
    //printf("%d\n", tmp);
    //tmp = 81 + (max_vel - 21) + 1;
    //printf("%d\n", tmp);
    std::vector<double> v_latt;
    std::vector<double> t_latt(sim_time + 1);
    for (double i = 0; i < tl - 200; i+=5) {x_latt.push_back(i);}
    for (double i = tl - 200; i < tl - 10; i+=1) {x_latt.push_back(i);}
    for (double i = tl - 10; i < tl - 1; i+=0.5) {x_latt.push_back(i);}
    for (double i = tl - 1; i <= tl + 1; i+=0.1) {x_latt.push_back(i);}
    x_latt.push_back(tl + 1);
    x_latt.push_back(tl + 10);
    //printf("%d\n", x_latt.size());

    for (double i = 0; i < 20; i+=0.25) {v_latt.push_back(i);}
    for (double i = 20; i <= max_vel; i++) {v_latt.push_back(i);}
    //printf("%d\n", v_latt.size());
    
    for (double i = 0; i <= sim_time; i++) {
        t_latt[i] = i;
    }

    const int x_size = x_latt.size();
    const int v_size = v_latt.size();
    const int t_size = t_latt.size();

    std::vector<std::vector<std::vector<double>>> e_opt(t_latt.size(), std::vector<std::vector<double>>(x_latt.size(), std::vector<double>(v_latt.size(), 0)));
    std::vector<std::vector<std::vector<double>>> u_opt(t_latt.size(), std::vector<std::vector<double>>(x_latt.size(), std::vector<double>(v_latt.size(), 0)));

    printf("t %d x %d v %d\n", e_opt.size(), e_opt[0].size(), e_opt[0][0].size());

    printf("start backward search.\n");

    for (int i = 0; i < x_size; i++) {
        for (int j = 0; j < v_size; j++) {
            e_opt[t_size - 1][i][j] = penalty_x * std::abs(x_f - x_latt[i]) + penalty_v * std::abs(v_f - v_latt[j]);
        }
    }

    int x_latt_size = x_latt.size();
    int v_latt_size = v_latt.size();

    unsigned int nproc = std::thread::hardware_concurrency();   // number of processer(実際はスレッド数)
    std::vector<std::thread*> threads(nproc);


    auto start = std::chrono::system_clock::now();
    for (int i = t_size - 2; i >= 0; i--) {
        printf("\rcalculating %d / %d", t_size - i, t_size);
        fflush(stdout);
        double dt = t_latt[i + 1] - t_latt[i];
                //std::cout << j << std::endl;
        for (int m = 0; m < nproc; m++) {
            threads[m] = new std::thread(assign, i, m, x_latt_size, nproc, dt, &x_latt, &v_latt, &t_latt, &e_opt, &u_opt);
        }
        for (std::thread *t : threads) {
            t->join();
        }
        for (std::thread *t : threads) {
            delete t;
        }

        /*
                if (i == 59) {
                    //std::cout << e_opt[i][j][k] << " " << u_opt[i][j][k] << std::endl;
                    break;
                }
                */
    }
    auto end = std::chrono::system_clock::now();
    auto dur = end - start;
    auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    std::cout << std::endl << msec << " ms" << std::endl;

/*
    for (double u : u_opt[t_size - 4][0]) {
        printf("%f\n ", u);
    }
    return 0;
  
    for (double u : u_opt[0][180]) {
        //printf("%f\n ", u);
    }
    puts("");
    for (double u : u_opt[50][100]) {
        printf("%f\n ", u)
    }
    */

    // forward search

    double energy=0;

    std::vector<double> u_solved(t_latt.size(), 0);
    std::vector<double> t_solved(t_latt.size(), 0);
    std::vector<double> v_solved(t_latt.size(), 0);
    v_solved[0] = v_i;
    std::vector<double> x_solved(t_latt.size(), 0);
    x_solved[0] = x_i;
    std::vector<double> e_solved(t_latt.size(), 0);
    std::vector<double> p_solved(t_latt.size(), 0);

    printf("Start of Forward Search \n");

    for (int i = 0; i < t_latt.size() - 1; i++) {   //start searching from the intial time
        double dt = t_latt[i+1] - t_latt[i];

        u_solved[i] = interpolate(&u_opt, i, x_solved[i], v_solved[i], &x_latt, &v_latt, x_latt_size, v_latt_size);    //linear interpolation by intepol program
        //printf("%d %f %f %f\n ", i, u_solved[i], x_solved[i], v_solved[i]);
      
        double next_x;
        double next_v;
        double dW;
        difsolve(&next_x, &next_v, &dW, x_solved[i], v_solved[i], dt, u_solved[i]); //solving the next state by difsolve program
      
        x_solved[i + 1] = next_x;
        v_solved[i + 1] = next_v;
        
        energy += dW;
        e_solved[i + 1] = energy;
        p_solved[i] = dW;
    }
    // 逆方向探索をすると、ある地点からの終状態へ行くまでの最適なノッチが全地点で計算されて、初期値を決めれば計算されたノッチを用いて軌道を生成できる
    // したがってロバスト性が高い ということか

    u_solved[t_latt.size() - 1] = 0;

    printf("End of Forward Search \n");


    for (double u : u_solved) {
        //std::cout << u << "\n";
    }
    puts("");

    printf("x_error is %f   v_error is %f \n",x_solved[t_size - 1]-x_f,v_solved[t_size - 1]-v_f);
    printf("Energy is %f kWh \n", energy);

    double v_solved_max = std::ceil(*max_element(v_solved.begin(), v_solved.end())/10)*10;

    FILE *gp;
    if ((gp = popen(GNUPLOT, "w")) == NULL) {
        fprintf(stderr, "ERROR: cannot open \"%s\".\n", GNUPLOT);
        exit(1);
    }

    fprintf(gp, "set multiplot\n\
    set lmargin screen 0.1\n\
    set rmargin screen 0.9\n\
    set tmargin screen 0.9\n\
    set bmargin screen 0.7\n");

    fprintf(gp, "set xrange [%f:%f]\n", 0., x_f);
    fprintf(gp, "set yrange [%f:%f]\n", *std::min_element(u_notch.begin(), u_notch.end()), *std::max_element(u_notch.begin(), u_notch.end()));
    fprintf(gp, "plot '-' with lines linetype 1 title \"Notch\"\n");

    for (int i = 0; i < t_size; i++) {
        fprintf(gp, "%f\t%f\n", x_solved[i], u_solved[i]);
    }
    fprintf(gp, "e\n");

    fprintf(gp, "set lmargin screen 0.1\n\
    set rmargin screen 0.9\n\
    set tmargin screen 0.6\n\
    set bmargin screen 0.4\n");

    fprintf(gp, "set xrange [%f:%f]\n", 0., x_f);
    fprintf(gp, "set yrange [%f:%f]\n", *std::min_element(e_solved.begin(), e_solved.end()), *std::max_element(e_solved.begin(), e_solved.end()));
    fprintf(gp, "plot '-' with lines linetype 1 title \"energy\"\n");

    for (int i = 0; i < t_size; i++) {
        fprintf(gp, "%f\t%f\n", x_solved[i], e_solved[i]);
    }
    fprintf(gp, "e\n");

    fprintf(gp, "set lmargin screen 0.1\n\
    set rmargin screen 0.9\n\
    set tmargin screen 0.3\n\
    set bmargin screen 0.1\n\
    set ylabel \"velocity\"\n\
    set ytics nomirror\n\
    set y2label \"time\"\n\
    set y2tics\n");

    fprintf(gp, "set xrange [%f:%f]\n", 0., x_f);
    fprintf(gp, "set yrange [%f:%f]\n", *std::min_element(v_solved.begin(), v_solved.end()), *std::max_element(v_solved.begin(), v_solved.end()));
    fprintf(gp, "set y2range [%f:%f]\n", *std::min_element(t_solved.begin(), t_solved.end()), *std::max_element(t_solved.begin(), t_solved.end()));
    fprintf(gp, "plot '-' with lines linetype 1 title \"energy\"\n");

    for (int i = 0; i < t_size; i++) {
        fprintf(gp, "%f\t%f\n", x_solved[i], v_solved[i]);
    }
    fprintf(gp, "axis x1y1\n");
    fprintf(gp, "replot '-' with lines linetype 1 title \"time\"\n");
    for (int i = 0; i < t_size; i++) {
        fprintf(gp, "%f\t%f\n", x_solved[i], t_solved[i]);
    }
    fprintf(gp, "axis x1y2\n");
    fprintf(gp, "e\n");
/////////縦軸を２つにしてグラフ描画する方法がわからぬ
    if (pclose(gp) == EOF) {
        fprintf(stderr, "Error: cannot close \"%s\".\n", GNUPLOT);
        exit(2);
    }

    return 0;
}
