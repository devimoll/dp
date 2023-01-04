#include <iostream>
#include <cmath>
#include <sched.h>
#include <vector>
#include <chrono>
#include <algorithm>
#include <thread>

#define _USE_MATH_DEFINES
#define GNUPLOT "gnuplot -persist"

// 最初にgnuplotの存在を確認 あとデータはcsvでも出力
// 勾配と速度制限とグラフ描画の設定がマジックすぎる
// difsolve関数が長い

// constants that can be changed as needed
const double tl = 800; // 線路長 track length
const double max_vel = 40; // 最高速度 maximum velocity[km/h]
const double sim_time = 120;   // simulation time
const double weight = 32; // [ton]
// train's resistances correspond to train's speed [km/h]
const double a = 14.2974;
const double b = 0.1949;
const double c = 0.8095;
// Track gradient
const double p1 = 0;
const double p2 = 20;
const double p3 = 0;
const double z1 = 300;
const double z2 = 500;
const double z3 = 800;
const double g = 9.8; // gravitational constant [m/s^2]
const double l = 24;  // train's length[m]
std::vector<double> u_notch{-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
const double pen1 = std::pow(10, 6);
const double pen2 = std::pow(10, 100);
const double penalty_v = 1;
const double penalty_x = 1;
// スレッド数
const int nproc = std::thread::hardware_concurrency();   // number of processer(実際はスレッド数)
//const int nproc = 1;


// do not edit under from this line
void add_m_m(double*, const double*, const double*);
void sum_m_m(double*, const double*, const double*);
void mul_m_v(double*, const double*, const double*);
void add_v_v(double*, const double*, const double*);
void mul_s_v(double*, const double, const double*);
void mul_s_m(double*, const double, const double*);
void inv_m(double*, const double*);
void solve_next_state(double*, const double*, const double*, const double*, const double, const double);
double ms2kmh(double);
double kmh2ms(double);
void traction(double*, double*, const double, const double, const double);
void difsolve(double*, double*, double*, double, double, double, const double, const double);
double value2index(double, const std::vector<double>&, int);
double interpolate(const std::vector<std::vector<double>>&, double, double);
void calc(int, int, int);
void assign_calculation(int, int, int);
int calc_x_size();
int calc_v_size();
int calc_t_size();
void init_x_latt();
void init_v_latt();
void init_t_latt();
void do_backward_search();
void do_forward_search();
void draw_graphs();


const int notch_count = u_notch.size();
const int x_size = calc_x_size();
const int v_size = calc_v_size();
const int t_size = calc_t_size();
const double t_i = 0;
const double t_f = sim_time;
const double v_i = 0;
const double v_f = 0;
const double x_i = 0;
const double x_f = tl;
// 位置、速度、時間の格子 lattice for position, speed and time
std::vector<double> x_latt(x_size);
std::vector<double> v_latt(v_size);
std::vector<double> t_latt(t_size);

std::vector<std::vector<std::vector<double>>> e_opt(t_size, std::vector<std::vector<double>>(x_size, std::vector<double>(v_size, 0)));
std::vector<std::vector<std::vector<double>>> u_opt(t_size, std::vector<std::vector<double>>(x_size, std::vector<double>(v_size, 0)));
std::vector<std::vector<std::vector<double>>> i_opt(t_size, std::vector<std::vector<double>>(x_size, std::vector<double>(v_size, 0)));

std::vector<double> u_solved(t_size, 0);
std::vector<double> v_solved(t_size, 0);
std::vector<double> x_solved(t_size, 0);
std::vector<double> e_solved(t_size, 0);
std::vector<double> p_solved(t_size, 0);
std::vector<double> i_solved(t_size, 0);


int main()
{
    printf("t %d x %d v %d\n", t_size, x_size, v_size);

    init_x_latt();
    init_v_latt();
    init_t_latt();

    // backward search with processing time measurement
    // 処理時間を計測しつつ後ろ向き探索をする
    printf("start of backward search.\n");
    auto start_time = std::chrono::system_clock::now();
    do_backward_search();
    auto end_time = std::chrono::system_clock::now();
    printf("end of backward search.\n");
    auto dur = end_time - start_time;
    auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    std::cout << msec << " ms" << std::endl;

    // forward search
    printf("Start of Forward Search \n");
    do_forward_search();
    printf("End of Forward Search \n");

    // 誤差表示
    double x_error = x_solved.back() - x_f;
    double v_error = v_solved.back() - v_f;
    printf("x_error = %f\n"
        "v_error = %f\n"
        "energy = %f kWh\n",
        x_error, v_error, e_solved.back());

    draw_graphs();

    return 0;
}

// 全部2次元 NULLチェックも要素数チェックもないから注意
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

inline double ms2kmh(double v)
{
    return 3.6 * v;
}

inline double kmh2ms(double v)
{
    return v / 3.6;
}

double fc_voltage(double i)
{
    return -0.00001 * i + 600;
    return -0.4 * i + 600;
}

void traction(double *F, double *dFdv, const double u, const double v, const double current)
{
    double f_a = 3.3 * 1000 / 3.6; // [N/t]
    double v_a1 = kmh2ms(40);    // [km/3.5*1000/3.6h]->[m/s]
    double v_a2 = kmh2ms(80);   // [km/3.5*1000/3.6h]->[m/s]

    // regenative
    double f_b = 4 * 1000 / 3.6;    // [N/t]
    double v_b = kmh2ms(80);    // [km/h]->[m/s]

    *F = 0;
    *dFdv = 0;

    // 堀内さんP.18式(2.22), 野田さんP.21
    double Vn = 600;    // 電動機の定格電圧
    double Vin = fc_voltage(current);    // 入力電圧
    //printf("%f ", Vin);
    v_a1 *= (Vin / Vn);
    v_a2 *= (Vin / Vn);
    v_b *= (Vin / Vn);

    // f_a, f_bとuの積が先行研究でのf_aやf_0か
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

// dt秒後のxとvと、dt秒間のdWを計算して、第1~3引数に値を入れる。
void difsolve(double *next_x, double *next_v, double *dW, double curr_x, double curr_v, double dt, const double u, double current)
{
    curr_v = kmh2ms(curr_v);
    // calculate traction force
    // 引張力
    double F;
    double dFdv;
    traction(&F, &dFdv, u, curr_v, current);
    // calculation for train's resistance
    // 列車抵抗
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
    beta /= (1000 * weight); // [N/(m/s)] -> [N/ton/(m/s)] -> [N/kg/(m/s)]=[(m/s/s)/(m/s)]

    double A[4] = {0, 1, beta, alpha};  /*  0    1
                                          beta alpha*/
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

//if(u==1)
  //      printf("u=%f, %f %f %f %f %f %f", u, curr_xv[0], curr_xv[1], dt, F0, next_xv[0], next_xv[1]);

    double dx = next_xv[0] - curr_xv[0];
    double dv = next_xv[1] - curr_xv[1];
    double vn = next_xv[1]; //[m/s]
    double Fn, dev_null;
    traction(&Fn, &dev_null, u, curr_v, current); // dFdvは不要
    Fn *= weight;

/*
    double a = (dFdv * dv) / dt;
    double b = F;
    double c = dv / dt;
    double d = curr_v;
    */

    next_xv[1] = ms2kmh(next_xv[1]);
    // efficiency
    double eta_a = 0.9;
    //double eta_b = 0.8;
    // 回生なし
    double eta_b = 0;
    //

    ////energy calculation 今と次の瞬間の電力の平均に時間をかけて電力量
    if (u > 0) {
        *dW  = 0.5 * (F * curr_v + Fn * vn ) * dt / eta_a / 3600 / 1000 * 1.1;  //[m/s]*[N]*[s]=[J]=1/(3.6*10^6)[kWh]
        //*dW  = 0.5*(F*curr_v + (F+dFdv*dv)*vn ) * dt / eta / 3600 / 1000 *1.1;  
        //*dW  = (F + 0.5*dFdv*dv) * dx / eta / 3600 / 1000;  
        //*dW  =(a*c)/3*dt^3+(a*d+b*c)/2*dt^2+(b*d)*dt / eta / 3600 / 1000;  
    }
    else if (u == 0) {
        *dW  = 0;
    }
    else {
        *dW  =  0.5 * (F * curr_v + Fn * vn ) * dt * eta_b / 3600 / 1000 * 1.1;

        //*dW  = 0.5*(F*curr_v + (F+dFdv*dv)*vn ) * dt * eta / 3600 / 1000 *1.1 ;   
        //*dW  = (F + 0.5*dFdv*dv) * dx * eta / 3600 / 1000;     
        //*dW  = (a*c)/3*dt^3+(a*d+b*c)/2*dt^2+(b*d)*dt * eta / 3600 / 1000;     
    }

//printf("%f %f %f\n", *dW, next_xv(0), next_xv(1));
    *next_x = next_xv[0];
    *next_v = next_xv[1];
//if(u==1)
  //      printf("%f %f %f\n", *dW, next_xv[0], next_xv[1]);
}

// 二分探索により、その値よりも小さいかつ最も近い格子のindexを求める
double value2index(double n, const std::vector<double> &latt, int latt_size)
{
    if (n >= latt.back()) {
        return latt_size - 1;
    }
    if (n <= latt[0]) {
        return 0;
    }

    int section_start = 0;
    int section_end = latt_size - 1;
    int d = section_end - section_start;
    int i;
    while (d > 1) {
        i = section_start + (d >> 1);
        //
        // sect(0)-----i-----sect(1)
        //
        if (latt[i] < n) {
            section_start = i;
        }
        else if (n < latt[i]) {
            section_end = i;
        }
        else {
            return i;
        }
        //printf("%d %d %d %f\n", i, section_start, section_end, latt[i]);
        d = section_end - section_start;
    }

    return latt[i] <= n ? i : i - 1;
}

double interpolate(const std::vector<std::vector<double>> &label, double next_x, double next_v)
{
    int k = value2index(next_x, x_latt, x_size);
    int j = value2index(next_v, v_latt, v_size);
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
    double u_11 = label[k][j];   // (k,j)
    double u_12 = label[std::min(k + 1, x_size - 1)][j]; // (k+1,j)
    double u_21 = label[k][std::min(j + 1, v_size - 1)]; // (k,j+1)
    double u_22 = label[std::min(k + 1, x_size - 1)][std::min(j + 1, v_size - 1)];   // (k+1,j+1)

    double g1;
    double g2;

    if (k == x_size - 1) {
        g1 = u_11;  // u_12 = u_11
        g2 = u_21;  // u_22 = u_21
    }
    else {
        g1 = ((u_12 * (next_x - x_latt[k])) + (u_11 * (x_latt[std::min(k + 1, x_size - 1)] - next_x))) / (x_latt[std::min(k + 1, x_size - 1)] - x_latt[k]);
        g2 = ((u_22 * (next_x - x_latt[k])) + (u_21 * (x_latt[std::min(k + 1, x_size - 1)] - next_x))) / (x_latt[std::min(k + 1, x_size - 1)] - x_latt[k]);
    }

    double u;

    if (j == v_size - 1) {
        u = g1; // g1 = g2
    }
    else {
        u = ((g2 * (next_v - v_latt[j])) + (g1 * (v_latt[std::min(j + 1, v_size - 1)] - next_v))) / (v_latt[std::min(j + 1, v_size - 1)] - v_latt[j]);
    }

    return u;
}

void calc(int step, int j, int dt)
{
    for (int k = 0; k < v_size; ++k) {
        double e_opt0 = pen1;
        double u_opt0 = 0;
        double e_opt1 = 0;
        double i_opt0 = 0;
        double i_opt1 = 0;
        for (int u = 0; u < notch_count; ++u) {
            double next_x;
            double next_v;
            double dW;
            difsolve(&next_x, &next_v, &dW, x_latt[j], v_latt[k], dt, u_notch[u], i_opt[step + 1][j][k]);
            double hokidenryoku = 0.01;
            i_opt1 = (dW + hokidenryoku) / (dt / 3600.) * 1000 / fc_voltage(i_opt[step + 1][j][k]); // [A]
            // 電流が範囲外
            if (i_opt1 < 0 || 500 < i_opt1) {
                e_opt1 = pen2;
            }
            //else if (i_opt1 < 20) {
            //    e_opt1 = pen1;
            //}
            //else if (i_opt1 > 200) {
            //    e_opt1 = pen1;
            //}
            // 電圧が低すぎる
            else if (fc_voltage(i_opt1) < 400) {
                e_opt1 = pen2;
            }
            //printf("\n%f, %f\n", dW, i_opt1);
            // speed limitation線路の形状による速度制限
            else if (x_latt[j] < 100 && v_latt[k] >= 60) {
                e_opt1 = pen1;   
            }
            else if (x_latt[j] >= 100 && x_latt[j] < 450 && v_latt[k] >= 90) {
                e_opt1 = pen1;   
            }
            else if (x_latt[j] >= 450 && x_latt[j] < 550 && v_latt[k] >= 95) {
                e_opt1 = pen1;   
            }
            else if (x_latt[j] >= 550 && x_latt[j] < 800 && v_latt[k] >= 90) {
                e_opt1 = pen1;   
            }
            /*  
            else if ((x_latt[j] > x_latt[x_latt.size() - 1])) {
                e_opt1 = pen2;   
            }*/
            else {
                e_opt1 = dW + interpolate(e_opt[step + 1], next_x, next_v); // linear interpolation by interpol program
            }
            if (e_opt1 < e_opt0) { // finding the optimum notch by evaluating the energy 最小を記憶
                e_opt0 = e_opt1;
                u_opt0 = u_notch[u];
                i_opt0 = i_opt1;
            }
        }
        // 消費エネルギーが最小になるノッチと、そのときの消費エネルギーが求まった。
        e_opt[step][j][k] = e_opt0;
        u_opt[step][j][k] = u_opt0;
        i_opt[step][j][k] = i_opt0;
        //printf("\n%f\n", i_opt[step][j][k]);
    }
}

// 各スレッドに計算を割り当てる
void assign_calculation(int step, int n, int dt)
{
    for (int j = n; j < x_size; j += nproc) {
        calc(step, j, dt);
    }
}

int calc_x_size()
{
    /*
    0 <= x < tl - 200: 5ずつ(tl - 200) / 5 個
    tl - 200 <= x < tl - 10: 1ずつ 190個
    tl - 10 <= x < tl - 1: 0.5ずつ 18個
    tl - 1 <= x <= tl + 1: 0.1ずつ 21個
    x = tl + 10: 1個
    */
    if (tl < 200) {
        printf("track length must be greater than 200.\n");
        exit(1);
    }

    return (tl - 200) * 0.2 + 230;
}

int calc_v_size()
{
    if (max_vel < 20) {
        printf("max velocity must be greater than 20.\n");
        exit(1);
    }
    return max_vel + 61;  // (v - 20) + 1 + (20 * 4)
}

int calc_t_size()
{
    return sim_time + 1;
}

void init_x_latt()
{
    double i = 0;
    for (double &x : x_latt) {
        x = i;
        if (i < tl - 200) {
            i += 5;
        }
        else if (i < tl - 10) {
            ++i;
        }
        else if (i < tl - 1) {
            i += 0.5;
        }
        else if (i < tl + 1) {
            i += 0.1;
        }
        else {
            x_latt.back() = tl + 10;
            break;
        }
    }
}

void init_v_latt()
{
    double i = 0;
    for (double &v : v_latt) {
        v = i;
        
        if (i < 20) {
            i += 0.25;
        }
        else {
            ++i;
        }
    }
}

void init_t_latt()
{
    double i = 0;
    for (double &t : t_latt) {
        t = i;
        ++i;
    }
}

void do_backward_search()
{
    std::vector<std::thread*> threads(nproc);

    // 終状態
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < v_size; ++j) {
            e_opt.back()[i][j] = penalty_x * std::abs(x_f - x_latt[i]) + penalty_v * std::abs(v_f - v_latt[j]);
            i_opt.back()[i][j] = 0;
        }
    }
    // 最後から2番目より最初まで探索
    for (int i = t_size - 2; i >= 0; --i) {
        printf("\rcalculating %d / %d", t_size - i, t_size);
        fflush(stdout);

        double dt = t_latt[i + 1] - t_latt[i];

        // multi threading
        for (int n = 0; n < nproc; ++n) {
            threads[n] = new std::thread(assign_calculation, i, n, dt);
        }
        for (std::thread *thread : threads) {
            thread->join();
        }
        for (std::thread *thread : threads) {
            delete thread;
        }
    }
    printf("\n");
}

void do_forward_search()
{
    v_solved[0] = v_i;
    x_solved[0] = x_i;
    i_solved[0] = 0;

    //start searching from the intial time
    for (int t = 0; t < t_size - 1; ++t) {
        double dt = t_latt[t + 1] - t_latt[t];
        u_solved[t] = interpolate(u_opt[t], x_solved[t], v_solved[t]);
        difsolve(&x_solved[t + 1], &v_solved[t + 1], &p_solved[t], x_solved[t], v_solved[t], dt, u_solved[t], i_solved[t]);
        e_solved[t + 1] = e_solved[t] + p_solved[t];
        double hokidenryoku = 0.01;
        i_solved[t + 1] = (p_solved[t] + hokidenryoku) / (dt / 3600) * 1000 / fc_voltage(i_solved[t]); // [A]
    }
    // 逆方向探索をすると、ある地点からの終状態へ行くまでの最適なノッチが全地点で計算されて、初期値を決めれば計算されたノッチを用いて軌道を生成できる
    // したがってロバスト性が高い ということか

    u_solved.back() = 0;
}

void draw_graphs()
{
    double x_error = x_solved.back() - x_f;
    double v_error = v_solved.back() - v_f;

    // 軸の最大値 10刻みにする
    int v_axis_max = std::ceil(*max_element(v_solved.begin(), v_solved.end()) * 0.1) * 10;
    // 1刻みにする
    int e_axis_max = std::ceil(*max_element(e_solved.begin(), e_solved.end()));

    FILE *gp;
    if ((gp = popen(GNUPLOT, "w")) == NULL) {
        fprintf(stderr, "ERROR: cannot open \"%s\".\n", GNUPLOT);
        exit(1);
    }
    fprintf(gp, "set terminal wxt size 600,800\n");

    // 誤差やペナルティなど
    fprintf(gp, "set label 1 at screen 0.1, 0.95 \"energy= %f kWh, p_x= %f, p_v= %f, e_x= %f, e_v= %f\"\n",
    e_solved.back(), penalty_x, penalty_v, x_error, v_error);

    // 運転曲線
    fprintf(gp, "set multiplot\n\
    set lmargin screen 0.1\n\
    set rmargin screen 0.88\n\
    set tmargin screen 0.923\n\
    set bmargin screen 0.77\n\
    set xlabel \"Distance [m]\"\n\
    set xtics nomirror\n\
    set ylabel \"Speed [km/h]\"\n\
    set ytics nomirror\n\
    set y2tics nomirror\n\
    set y2label \"Time [s]\"\n\
    ");
    fprintf(gp, "set xrange [%d:%f]\n", 0, x_f);
    fprintf(gp, "set yrange [%d:%d]\n", 0, v_axis_max);
    fprintf(gp, "set y2range [%d:%f]\n", 0, sim_time);
    fprintf(gp, "plot '-' with lines linetype 1 title \"Speed\" axis x1y1, ");
    fprintf(gp, "'-' with lines linetype 2 title \"Time\" axis x1y2\n");
    for (int i = 0; i < t_size; ++i) {
        fprintf(gp, "%f\t%f\n", x_solved[i], v_solved[i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < t_size; ++i) {
        fprintf(gp, "%f\t%f\n", x_solved[i], t_latt[i]);
    }
    fprintf(gp, "e\n");

    // エネルギ
    fprintf(gp, "set lmargin screen 0.1\n\
    set lmargin screen 0.1\n\
    set rmargin screen 0.88\n\
    set tmargin screen 0.69\n\
    set bmargin screen 0.54\n\
    set ylabel \"Energy [kWh]\"\n\
    set noy2label\n\
    ");
    fprintf(gp, "set xrange [%d:%f]\n", 0, x_f);
    fprintf(gp, "set yrange [%d:%d]\n", 0, e_axis_max);
    fprintf(gp, "plot '-' with lines linetype 1 notitle\n");

    for (int i = 0; i < t_size; ++i) {
        fprintf(gp, "%f\t%f\n", x_solved[i], e_solved[i]);
    }
    fprintf(gp, "e\n");

    // 電流、電圧
    std::vector<double> voltage_solved(t_size);
    for (int i = 0; i < t_size; ++i) {
        voltage_solved[i] = fc_voltage(i_solved[i]);
    }
    fprintf(gp, "set lmargin screen 0.1\n\
    set lmargin screen 0.1\n\
    set rmargin screen 0.88\n\
    set tmargin screen 0.46\n\
    set bmargin screen 0.31\n\
    set ylabel \"Current [A]\"\n\
    set y2label \"Voltage [V]\"\n\
    ");
    fprintf(gp, "set xrange [%d:%f]\n", 0, x_f);
    fprintf(gp, "set yrange [%f:%f]\n", *std::min_element(i_solved.begin(), i_solved.end()), *std::max_element(i_solved.begin(), i_solved.end()));
    fprintf(gp, "set y2range [%f:%f]\n", *std::min_element(voltage_solved.begin(), voltage_solved.end()), *std::max_element(voltage_solved.begin(), voltage_solved.end()));
    fprintf(gp, "plot '-' with lines linetype 1 title \"Current\" axis x1y1, ");
    fprintf(gp, "'-' with lines linetype 2 title \"Voltage\" axis x1y2\n");

    for (int i = 0; i < t_size; ++i) {
        fprintf(gp, "%f\t%f\n", x_solved[i], i_solved[i]);
    }
    fprintf(gp, "e\n");
    for (int i = 0; i < t_size; i++) {
        fprintf(gp, "%f\t%f\n", x_solved[i], voltage_solved[i]);
    }
    fprintf(gp, "e\n");

    // ノッチ
    fprintf(gp, "set lmargin screen 0.1\n\
    set rmargin screen 0.9\n\
    set lmargin screen 0.1\n\
    set rmargin screen 0.88\n\
    set tmargin screen 0.23\n\
    set bmargin screen 0.077\n\
    set ylabel \"Notch\"\n\
    set xzeroaxis\n\
    set noy2tics\n\
    set noy2label\n\
    ");
    fprintf(gp, "set xrange [%d:%f]\n", 0, x_f);
    fprintf(gp, "set yrange [%f:%f]\n", *std::min_element(u_notch.begin(), u_notch.end()), *std::max_element(u_notch.begin(), u_notch.end()));
    fprintf(gp, "plot '-' with lines linetype 1 notitle\n");

    for (int i = 0; i < t_size; i++) {
        fprintf(gp, "%f\t%f\n", x_solved[i], u_solved[i]);
    }
    fprintf(gp, "e\n");

    if (pclose(gp) == EOF) {
        fprintf(stderr, "Error: cannot close \"%s\".\n", GNUPLOT);
        exit(2);
    }
}

