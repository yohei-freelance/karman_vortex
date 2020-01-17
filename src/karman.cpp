#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

const double Re = 500.0;  // Reynolds Number
const double cfl = 0.15;  // CFL Number

/*SOR Pamameters*/
const double omegap = 0.500;//1.00;
const int maxitp = 100;
const double errorp = 0.00001;

/* No. of Time Steps*/
const int nlast = 500;

/* set x-grid parameters*/
const int mx = 401;   // x軸格子点数(1~401),x=30.0がmx=401に対応
const int i_1 = 169;
const int i_2 = 173;   // x軸格子点数基準での,角柱の左端
const int i_3 = 167;  // x軸格子点数基準での,角柱の右端
const int i_4 = 169;
const int i_5 = 165;
const int i_6 = 167;
const int i_7 = 163;
const int i_8 = 165;
const int i_9 = 161;
const int i_10 = 163;
const int i_11 = 173;
const int i_12 = 175;
const int i_13 = 175;
const int i_14 = 177;
const int i_15 = 177;
const int i_16 = 179;
const int i_17 = 179;
const int i_18 = 181;

/* set y-grid parameters*/
const int my = 201;   // y軸格子点数(1~201)
const int j_1 = 91;
const int j_2 = 111;   // y軸格子点数基準での,角柱の下端
const int j_3 = 92;  // y軸格子点数基準での,角柱の上端
const int j_4 = 110;
const int j_5 = 93;
const int j_6 = 109;
const int j_7 = 95;
const int j_8 = 107;   // y軸格子点数基準での,角柱の下端
const int j_9 = 99;  // y軸格子点数基準での,角柱の上端
const int j_10 = 103;
const int j_11 = 92;
const int j_12 = 110;
const int j_13 = 93;
const int j_14 = 109;   // y軸格子点数基準での,角柱の下端
const int j_15 = 95;  // y軸格子点数基準での,角柱の上端
const int j_16 = 107;
const int j_17 = 99;
const int j_18 = 103;

/* set delta x,y,t*/
const double dx = 1.0 / (i_2 - i_1);
const double dy = 1.0 / (j_2 - j_1);
const double dt = cfl * fmin(dx, dy);

/*配列の定義*/
double x[mx + 1];
double y[my + 1];
float u[mx + 2][my + 2], v[mx + 2][my + 2], p[mx + 2][my + 2];

void make_xygrid(double x[mx + 1], double y[my + 1])
{
	double icent = (i_1 + i_2) / 2;
	double jcent = (j_1 + j_2) / 2;
	for (int i = 1; i <= mx; i++)
	{
		x[i] = dx * double(i - icent);
	}
	for (int j = 1; j <= my; j++)
	{
		y[j] = dy * double(j - jcent);
	}
    double icenti = (i_3 + i_4) / 2;
	double jcenti = (j_3 + j_4) / 2;
	for (int i = 1; i <= mx; i++)
	{
		x[i] = dx * double(i - icenti);
	}
	for (int j = 1; j <= my; j++)
	{
		y[j] = dy * double(j - jcenti);
	}
    double icenta = (i_5 + i_6) / 2;
	double jcenta = (j_5 + j_6) / 2;
	for (int i = 1; i <= mx; i++)
	{
		x[i] = dx * double(i - icenta);
	}
	for (int j = 1; j <= my; j++)
	{
		y[j] = dy * double(j - jcenta);
	}
	double icenta1 = (i_7 + i_8) / 2;
	double jcenta1 = (j_7 + j_8) / 2;
	for (int i = 1; i <= mx; i++)
	{
		x[i] = dx * double(i - icenta1);
	}
	for (int j = 1; j <= my; j++)
	{
		y[j] = dy * double(j - jcenta1);
	}
	double icenta2 = (i_9 + i_10) / 2;
	double jcenta2 = (j_9 + j_10) / 2;
	for (int i = 1; i <= mx; i++)
	{
		x[i] = dx * double(i - icenta2);
	}
	for (int j = 1; j <= my; j++)
	{
		y[j] = dy * double(j - jcenta2);
	}
	double icenta3 = (i_11 + i_12) / 2;
	double jcenta3 = (j_11 + j_12) / 2;
	for (int i = 1; i <= mx; i++)
	{
		x[i] = dx * double(i - icenta3);
	}
	for (int j = 1; j <= my; j++)
	{
		y[j] = dy * double(j - jcenta3);
	}
	double icenta4 = (i_13 + i_14) / 2;
	double jcenta4 = (j_13 + j_14) / 2;
	for (int i = 1; i <= mx; i++)
	{
		x[i] = dx * double(i - icenta4);
	}
	for (int j = 1; j <= my; j++)
	{
		y[j] = dy * double(j - jcenta4);
	}
	double icenta5 = (i_15 + i_16) / 2;
	double jcenta5 = (j_15 + j_16) / 2;
	for (int i = 1; i <= mx; i++)
	{
		x[i] = dx * double(i - icenta5);
	}
	for (int j = 1; j <= my; j++)
	{
		y[j] = dy * double(j - jcenta5);
	}
	double icenta6 = (i_17 + i_18) / 2;
	double jcenta6 = (j_17 + j_18) / 2;
	for (int i = 1; i <= mx; i++)
	{
		x[i] = dx * double(i - icenta6);
	}
	for (int j = 1; j <= my; j++)
	{
		y[j] = dy * double(j - jcenta6);
	}
}

void init_condition(float u[mx + 2][my + 2], float v[mx + 2][my + 2],
                    float p[mx + 2][my + 2])
{
	for (int i = 1; i <= mx; i++)
	{
		for (int j = 1; j <= my; j++)
		{
			u[i][j] = 1.0;
			v[i][j] = 0.05;
			p[i][j] = 0.0;
		}
	}
}

void bcfor_p(float p[mx + 2][my + 2])
{
	for (int j = 1; j <= my; j++)
	{
		p[1][j] = 0.0;
		p[mx][j] = 0.0;
	}
	for (int i = 1; i <= mx; i++)
	{
		p[i][1] = 0.0;
		p[i][my] = 0.0;
	}
	p[i_1][j_1] = p[i_1 - 1][j_1 - 1];
	p[i_1][j_2] = p[i_1 - 1][j_2 + 1];
	p[i_2][j_1] = p[i_2 + 1][j_1 - 1];
	p[i_2][j_2] = p[i_2 + 1][j_2 + 1];
    p[i_3][j_3] = p[i_3 - 1][j_3 - 1];
	p[i_3][j_4] = p[i_3 - 1][j_4 + 1];
	p[i_4][j_3] = p[i_4 + 1][j_3 - 1];
	p[i_4][j_4] = p[i_4 + 1][j_4 + 1];
    p[i_5][j_5] = p[i_5 - 1][j_5 - 1];
	p[i_5][j_6] = p[i_5 - 1][j_6 + 1];
	p[i_6][j_5] = p[i_6 + 1][j_5 - 1];
	p[i_6][j_6] = p[i_6 + 1][j_6 + 1];
	p[i_7][j_7] = p[i_7 - 1][j_7 - 1];
	p[i_7][j_8] = p[i_7 - 1][j_8 + 1];
	p[i_8][j_7] = p[i_8 + 1][j_7 - 1];
	p[i_8][j_8] = p[i_8 + 1][j_8 + 1];
	p[i_9][j_9] = p[i_9 - 1][j_9 - 1];
	p[i_9][j_10] = p[i_9 - 1][j_10 + 1];
	p[i_10][j_9] = p[i_10 + 1][j_9 - 1];
	p[i_10][j_10] = p[i_10 + 1][j_10 + 1];
	p[i_11][j_11] = p[i_11 - 1][j_11 - 1];
	p[i_11][j_12] = p[i_11 - 1][j_12 + 1];
	p[i_12][j_11] = p[i_12 + 1][j_11 - 1];
	p[i_12][j_12] = p[i_12 + 1][j_12 + 1];
	p[i_13][j_13] = p[i_13 - 1][j_13 - 1];
	p[i_13][j_14] = p[i_13 - 1][j_14 + 1];
	p[i_14][j_13] = p[i_14 + 1][j_13 - 1];
	p[i_14][j_14] = p[i_14 + 1][j_14 + 1];
	p[i_15][j_15] = p[i_15 - 1][j_15 - 1];
	p[i_15][j_16] = p[i_15 - 1][j_16 + 1];
	p[i_16][j_15] = p[i_16 + 1][j_15 - 1];
	p[i_16][j_16] = p[i_16 + 1][j_16 + 1];
	p[i_17][j_17] = p[i_17 - 1][j_17 - 1];
	p[i_17][j_18] = p[i_17 - 1][j_18 + 1];
	p[i_18][j_17] = p[i_18 + 1][j_17 - 1];
	p[i_18][j_18] = p[i_18 + 1][j_18 + 1];
	for (int j = j_1 + 1; j <= j_2 - 1; j++)
	{
		p[i_1][j] = p[i_1 - 1][j];
		p[i_2][j] = p[i_2 + 1][j];
	}
	for (int i = i_1 + 1; i <= i_2 - 1; i++)
	{
		p[i][j_1] = p[i][j_1 - 1];
		p[i][j_2] = p[i][j_2 + 1];
	}
    for (int j = j_3 + 1; j <= j_4 - 1; j++)
	{
		p[i_3][j] = p[i_3 - 1][j];
		p[i_4][j] = p[i_4 + 1][j];
	}
	for (int i = i_3 + 1; i <= i_4 - 1; i++)
	{
		p[i][j_3] = p[i][j_3 - 1];
		p[i][j_4] = p[i][j_4 + 1];
	}
    for (int j = j_5 + 1; j <= j_6 - 1; j++)
	{
		p[i_5][j] = p[i_5 - 1][j];
		p[i_6][j] = p[i_6 + 1][j];
	}
	for (int i = i_5 + 1; i <= i_6 - 1; i++)
	{
		p[i][j_5] = p[i][j_5 - 1];
		p[i][j_6] = p[i][j_6 + 1];
	}
	for (int j = j_7 + 1; j <= j_8 - 1; j++)
	{
		p[i_7][j] = p[i_7 - 1][j];
		p[i_8][j] = p[i_8 + 1][j];
	}
	for (int i = i_7 + 1; i <= i_8 - 1; i++)
	{
		p[i][j_7] = p[i][j_7 - 1];
		p[i][j_8] = p[i][j_8 + 1];
	}
	for (int j = j_9 + 1; j <= j_10 - 1; j++)
	{
		p[i_9][j] = p[i_9 - 1][j];
		p[i_10][j] = p[i_10 + 1][j];
	}
	for (int i = i_9 + 1; i <= i_10 - 1; i++)
	{
		p[i][j_9] = p[i][j_9 - 1];
		p[i][j_10] = p[i][j_10 + 1];
	}
	for (int j = j_11 + 1; j <= j_12 - 1; j++)
	{
		p[i_11][j] = p[i_11 - 1][j];
		p[i_12][j] = p[i_12 + 1][j];
	}
	for (int i = i_11 + 1; i <= i_12 - 1; i++)
	{
		p[i][j_11] = p[i][j_11 - 1];
		p[i][j_12] = p[i][j_12 + 1];
	}
	for (int j = j_13 + 1; j <= j_14 - 1; j++)
	{
		p[i_13][j] = p[i_13 - 1][j];
		p[i_14][j] = p[i_14 + 1][j];
	}
	for (int i = i_13 + 1; i <= i_14 - 1; i++)
	{
		p[i][j_13] = p[i][j_13 - 1];
		p[i][j_14] = p[i][j_14 + 1];
	}
	for (int j = j_15 + 1; j <= j_16 - 1; j++)
	{
		p[i_15][j] = p[i_15 - 1][j];
		p[i_16][j] = p[i_16 + 1][j];
	}
	for (int i = i_15 + 1; i <= i_16 - 1; i++)
	{
		p[i][j_15] = p[i][j_15 - 1];
		p[i][j_16] = p[i][j_16 + 1];
	}
	for (int j = j_17 + 1; j <= j_18 - 1; j++)
	{
		p[i_17][j] = p[i_17 - 1][j];
		p[i_18][j] = p[i_18 + 1][j];
	}
	for (int i = i_17 + 1; i <= i_18 - 1; i++)
	{
		p[i][j_17] = p[i][j_17 - 1];
		p[i][j_18] = p[i][j_18 + 1];
	}
}

void bcfor_v(float u[mx + 2][my + 2], float v[mx + 2][my + 2])
{
	for (int j = 1; j <= my; j++)
	{
		u[1][j] = 1.0;
		v[1][j] = 0.0;
		u[0][j] = 1.0;
		v[0][j] = 0.0;
		u[mx][j] = 2.0 * u[mx - 1][j] - u[mx - 2][j];
		v[mx][j] = 2.0 * v[mx - 1][j] - v[mx - 2][j];
		u[mx + 1][j] = 2.0 * u[mx][j] - u[mx - 1][j];
		v[mx + 1][j] = 2.0 * v[mx][j] - v[mx - 1][j];
	}
	for (int i = 1; i <= mx; i++)
	{
		u[i][1] = 2.0 * u[i][2] - u[i][3];
		v[i][1] = 2.0 * v[i][2] - v[i][3];
		u[i][0] = 2.0 * u[i][1] - u[i][2];
		v[i][0] = 2.0 * v[i][1] - v[i][2];
		u[i][my] = 2.0 * u[i][my - 1] - u[i][my - 2];
		v[i][my] = 2.0 * v[i][my - 1] - v[i][my - 2];
		u[i][my + 1] = 2.0 * u[i][my] - u[i][my - 1];
		v[i][my + 1] = 2.0 * v[i][my] - v[i][my - 1];
	}
	for (int i = i_1; i <= i_2; i++)
	{
		for (int j = j_1; j <= j_2; j++)
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
	}
    for (int i = i_3; i <= i_4; i++)
	{
		for (int j = j_3; j <= j_4; j++)
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
	}
    for (int i = i_5; i <= i_6; i++)
	{
		for (int j = j_5; j <= j_6; j++)
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
	}
	for (int i = i_7; i <= i_8; i++)
	{
		for (int j = j_7; j <= j_8; j++)
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
	}
	for (int i = i_9; i <= i_10; i++)
	{
		for (int j = j_9; j <= j_10; j++)
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
	}
	for (int i = i_11; i <= i_12; i++)
	{
		for (int j = j_11; j <= j_12; j++)
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
	}
	for (int i = i_13; i <= i_14; i++)
	{
		for (int j = j_13; j <= j_14; j++)
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
	}
	for (int i = i_15; i <= i_16; i++)
	{
		for (int j = j_15; j <= j_16; j++)
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
	}
	for (int i = i_17; i <= i_18; i++)
	{
		for (int j = j_17; j <= j_18; j++)
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
	}
}

void poisson_eq(float u[mx + 2][my + 2], float v[mx + 2][my + 2],
                float p[mx + 2][my + 2], double dx, double dy, double dt)
{
	double rhs[mx + 2][my + 2];
	for (int i = 2; i <= mx - 1; i++)
	{
		for (int j = 2; j <= my - 1; j++)
		{
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6) || (i_7 < i && i < i_8 && j_7 < j && j < j_8) || (i_9 < i && i < i_10 && j_9 < j && j < j_10) || (i_11 < i && i < i_12 && j_11 < j && j < j_12) || (i_13 < i && i < i_14 && j_13 < j && j < j_14) || (i_15 < i && i < i_16 && j_15 < j && j < j_16) || (i_17 < i && i < i_18 && j_17 < j && j < j_18)) 
			{
				continue;
			}
			float ux = (u[i + 1][j] - u[i - 1][j]) / (2.0 * dx);
			float uy = (u[i][j + 1] - u[i][j - 1]) / (2.0 * dy);
			float vx = (v[i + 1][j] - v[i - 1][j]) / (2.0 * dx);
			float vy = (v[i][j + 1] - v[i][j - 1]) / (2.0 * dy);
			rhs[i][j] = (ux + vy) / dt - (ux * ux + 2.0 * uy * vx + vy * vy);
		}
	}

	for (int itr = 1; itr <= maxitp; itr++)
	{
		double res = 0.0;
		for (int i = 2; i <= mx - 1; i++)
		{
			for (int j = 2; j <= my - 1; j++)
			{
				if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6) || (i_7 < i && i < i_8 && j_7 < j && j < j_8) || (i_9 < i && i < i_10 && j_9 < j && j < j_10) || (i_11 < i && i < i_12 && j_11 < j && j < j_12) || (i_13 < i && i < i_14 && j_13 < j && j < j_14) || (i_15 < i && i < i_16 && j_15 < j && j < j_16) || (i_17 < i && i < i_18 && j_17 < j && j < j_18))
				{
					continue;
				}
				float dp = (p[i + 1][j] + p[i - 1][j]) / (dx * dx) +
				            (p[i][j + 1] + p[i][j - 1]) / (dy * dy) - rhs[i][j];
				dp = dp / (2.0 / (dx * dx) + 2.0 / (dy * dy)) - p[i][j];
				res += dp * dp;
				p[i][j] = p[i][j] + omegap * dp;
			}
		}
		bcfor_p(p);
		res = sqrt(res / double(mx * my));
		if (res < errorp) break;
	}
}

void velocity_eq(float u[mx + 2][my + 2], float v[mx + 2][my + 2],
                 float p[mx + 2][my + 2], double dx, double dy, double dt)
{
	float urhs[mx + 2][my + 2], vrhs[mx + 2][my + 2];
	for (int i = 2; i <= mx - 1; i++)
	{
		for (int j = 2; j <= my - 1; j++)
		{
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6) || (i_7 < i && i < i_8 && j_7 < j && j < j_8) || (i_9 < i && i < i_10 && j_9 < j && j < j_10) || (i_11 < i && i < i_12 && j_11 < j && j < j_12) || (i_13 < i && i < i_14 && j_13 < j && j < j_14) || (i_15 < i && i < i_16 && j_15 < j && j < j_16) || (i_17 < i && i < i_18 && j_17 < j && j < j_18))
			{
				continue;
			}
			urhs[i][j] = -(p[i + 1][j] - p[i - 1][j]) / (2.0 * dx);
			vrhs[i][j] = -(p[i][j + 1] - p[i][j - 1]) / (2.0 * dy);
		}
	}

	for (int i = 2; i <= mx - 1; i++)
	{
		for (int j = 2; j <= my - 1; j++)
		{
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6) || (i_7 < i && i < i_8 && j_7 < j && j < j_8) || (i_9 < i && i < i_10 && j_9 < j && j < j_10) || (i_11 < i && i < i_12 && j_11 < j && j < j_12) || (i_13 < i && i < i_14 && j_13 < j && j < j_14) || (i_15 < i && i < i_16 && j_15 < j && j < j_16) || (i_17 < i && i < i_18 && j_17 < j && j < j_18))
			{
				continue;
			}
			urhs[i][j] =
				urhs[i][j] +
				(u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / (Re * dx * dx) +
				(u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / (Re * dy * dy);
			vrhs[i][j] =
				vrhs[i][j] +
				(v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) / (Re * dx * dx) +
				(v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) / (Re * dy * dy);
		}
	}

	for (int j = j_1 + 1; j <= j_2 - 1; j++)
	{
		u[i_1 + 1][j] = 2.0 * u[i_1][j] - u[i_1 - 1][j];
		u[i_2 - 1][j] = 2.0 * u[i_2][j] - u[i_2 + 1][j];
		v[i_1 + 1][j] = 2.0 * v[i_1][j] - v[i_1 - 1][j];
		v[i_2 - 1][j] = 2.0 * v[i_2][j] - v[i_2 + 1][j];
	}
    for (int j = j_3 + 1; j <= j_4 - 1; j++)
	{
		u[i_3 + 1][j] = 2.0 * u[i_3][j] - u[i_3 - 1][j];
		u[i_4 - 1][j] = 2.0 * u[i_4][j] - u[i_4 + 1][j];
		v[i_3 + 1][j] = 2.0 * v[i_3][j] - v[i_3 - 1][j];
		v[i_4 - 1][j] = 2.0 * v[i_4][j] - v[i_4 + 1][j];
	}
    for (int j = j_5 + 1; j <= j_6 - 1; j++)
	{
		u[i_5 + 1][j] = 2.0 * u[i_5][j] - u[i_5 - 1][j];
		u[i_6 - 1][j] = 2.0 * u[i_6][j] - u[i_6 + 1][j];
		v[i_5 + 1][j] = 2.0 * v[i_5][j] - v[i_5 - 1][j];
		v[i_6 - 1][j] = 2.0 * v[i_6][j] - v[i_6 + 1][j];
	}
	for (int j = j_7 + 1; j <= j_8 - 1; j++)
	{
		u[i_7 + 1][j] = 2.0 * u[i_7][j] - u[i_7 - 1][j];
		u[i_8 - 1][j] = 2.0 * u[i_8][j] - u[i_8 + 1][j];
		v[i_7 + 1][j] = 2.0 * v[i_7][j] - v[i_7 - 1][j];
		v[i_8 - 1][j] = 2.0 * v[i_8][j] - v[i_8 + 1][j];
	}
	for (int j = j_9 + 1; j <= j_10 - 1; j++)
	{
		u[i_9 + 1][j] = 2.0 * u[i_9][j] - u[i_9 - 1][j];
		u[i_10 - 1][j] = 2.0 * u[i_10][j] - u[i_10 + 1][j];
		v[i_9 + 1][j] = 2.0 * v[i_9][j] - v[i_9 - 1][j];
		v[i_10 - 1][j] = 2.0 * v[i_10][j] - v[i_10 + 1][j];
	}
	for (int j = j_11 + 1; j <= j_12 - 1; j++)
	{
		u[i_11 + 1][j] = 2.0 * u[i_11][j] - u[i_11 - 1][j];
		u[i_12 - 1][j] = 2.0 * u[i_12][j] - u[i_12 + 1][j];
		v[i_11 + 1][j] = 2.0 * v[i_11][j] - v[i_11 - 1][j];
		v[i_12 - 1][j] = 2.0 * v[i_12][j] - v[i_12 + 1][j];
	}
	for (int j = j_13 + 1; j <= j_14 - 1; j++)
	{
		u[i_13 + 1][j] = 2.0 * u[i_13][j] - u[i_13 - 1][j];
		u[i_14 - 1][j] = 2.0 * u[i_14][j] - u[i_14 + 1][j];
		v[i_13 + 1][j] = 2.0 * v[i_13][j] - v[i_13 - 1][j];
		v[i_14 - 1][j] = 2.0 * v[i_14][j] - v[i_14 + 1][j];
	}
	for (int j = j_15 + 1; j <= j_16 - 1; j++)
	{
		u[i_15 + 1][j] = 2.0 * u[i_15][j] - u[i_15 - 1][j];
		u[i_16 - 1][j] = 2.0 * u[i_16][j] - u[i_16 + 1][j];
		v[i_15 + 1][j] = 2.0 * v[i_15][j] - v[i_15 - 1][j];
		v[i_16 - 1][j] = 2.0 * v[i_16][j] - v[i_16 + 1][j];
	}
	for (int j = j_17 + 1; j <= j_18 - 1; j++)
	{
		u[i_17 + 1][j] = 2.0 * u[i_17][j] - u[i_17 - 1][j];
		u[i_18 - 1][j] = 2.0 * u[i_18][j] - u[i_18 + 1][j];
		v[i_17 + 1][j] = 2.0 * v[i_17][j] - v[i_17 - 1][j];
		v[i_18 - 1][j] = 2.0 * v[i_18][j] - v[i_18 + 1][j];
	}

	for (int i = 2; i <= mx - 1; i++)
	{
		for (int j = 2; j <= my - 1; j++)
		{
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6) || (i_7 < i && i < i_8 && j_7 < j && j < j_8) || (i_9 < i && i < i_10 && j_9 < j && j < j_10) || (i_11 < i && i < i_12 && j_11 < j && j < j_12) || (i_13 < i && i < i_14 && j_13 < j && j < j_14) || (i_15 < i && i < i_16 && j_15 < j && j < j_16) || (i_17 < i && i < i_18 && j_17 < j && j < j_18))
			{
				continue;
			}
			urhs[i][j] = urhs[i][j] -
			             u[i][j] *
			             (-u[i + 2][j] + 8.0 * (u[i + 1][j] - u[i - 1][j]) +
			              u[i - 2][j]) /
			             (12.0 * dx) -
			             abs(u[i][j]) *
			             (u[i + 2][j] - 4.0 * u[i + 1][j] + 6.0 * u[i][j] -
			              4.0 * u[i - 1][j] + u[i - 2][j]) /
			             (4.0 * dx);
			vrhs[i][j] = vrhs[i][j] -
			             u[i][j] *
			             (-v[i + 2][j] + 8.0 * (v[i + 1][j] - v[i - 1][j]) +
			              v[i - 2][j]) /
			             (12.0 * dx) -
			             abs(u[i][j]) *
			             (v[i + 2][j] - 4.0 * v[i + 1][j] + 6.0 * v[i][j] -
			              4.0 * v[i - 1][j] + v[i - 2][j]) /
			             (4.0 * dx);
		}
	}

	for (int i = i_1 + 1; i <= i_2 - 1; i++)
	{
		u[i][j_1 + 1] = 2.0 * u[i][j_1] - u[i][j_1 - 1];
		u[i][j_2 - 1] = 2.0 * u[i][j_2] - u[i][j_2 + 1];
		v[i][j_1 + 1] = 2.0 * v[i][j_1] - v[i][j_1 - 1];
		v[i][j_2 - 1] = 2.0 * v[i][j_2] - v[i][j_2 + 1];
	}
    for (int i = i_3 + 1; i <= i_4 - 1; i++)
	{
		u[i][j_3 + 1] = 2.0 * u[i][j_3] - u[i][j_3 - 1];
		u[i][j_4 - 1] = 2.0 * u[i][j_4] - u[i][j_4 + 1];
		v[i][j_3 + 1] = 2.0 * v[i][j_3] - v[i][j_3 - 1];
		v[i][j_4 - 1] = 2.0 * v[i][j_4] - v[i][j_4 + 1];
	}
    for (int i = i_5 + 1; i <= i_6 - 1; i++)
	{
		u[i][j_5 + 1] = 2.0 * u[i][j_5] - u[i][j_5 - 1];
		u[i][j_6 - 1] = 2.0 * u[i][j_6] - u[i][j_6 + 1];
		v[i][j_5 + 1] = 2.0 * v[i][j_5] - v[i][j_5 - 1];
		v[i][j_6 - 1] = 2.0 * v[i][j_6] - v[i][j_6 + 1];
	}
	for (int i = i_7 + 1; i <= i_8 - 1; i++)
	{
		u[i][j_7 + 1] = 2.0 * u[i][j_7] - u[i][j_7 - 1];
		u[i][j_8 - 1] = 2.0 * u[i][j_8] - u[i][j_8 + 1];
		v[i][j_7 + 1] = 2.0 * v[i][j_7] - v[i][j_7 - 1];
		v[i][j_8 - 1] = 2.0 * v[i][j_8] - v[i][j_8 + 1];
	}
	for (int i = i_9 + 1; i <= i_10 - 1; i++)
	{
		u[i][j_9 + 1] = 2.0 * u[i][j_9] - u[i][j_9 - 1];
		u[i][j_10 - 1] = 2.0 * u[i][j_10] - u[i][j_10 + 1];
		v[i][j_9 + 1] = 2.0 * v[i][j_9] - v[i][j_9 - 1];
		v[i][j_10 - 1] = 2.0 * v[i][j_10] - v[i][j_10 + 1];
	}
	for (int i = i_11 + 1; i <= i_12 - 1; i++)
	{
		u[i][j_11 + 1] = 2.0 * u[i][j_11] - u[i][j_11 - 1];
		u[i][j_12 - 1] = 2.0 * u[i][j_12] - u[i][j_12 + 1];
		v[i][j_11 + 1] = 2.0 * v[i][j_11] - v[i][j_11 - 1];
		v[i][j_12 - 1] = 2.0 * v[i][j_12] - v[i][j_12 + 1];
	}
	for (int i = i_13 + 1; i <= i_14 - 1; i++)
	{
		u[i][j_13+ 1] = 2.0 * u[i][j_13] - u[i][j_13 - 1];
		u[i][j_14 - 1] = 2.0 * u[i][j_14] - u[i][j_14 + 1];
		v[i][j_13 + 1] = 2.0 * v[i][j_13] - v[i][j_13 - 1];
		v[i][j_14 - 1] = 2.0 * v[i][j_14] - v[i][j_14 + 1];
	}
	for (int i = i_15 + 1; i <= i_16 - 1; i++)
	{
		u[i][j_15 + 1] = 2.0 * u[i][j_15] - u[i][j_15 - 1];
		u[i][j_16 - 1] = 2.0 * u[i][j_16] - u[i][j_16 + 1];
		v[i][j_15 + 1] = 2.0 * v[i][j_15] - v[i][j_15 - 1];
		v[i][j_16 - 1] = 2.0 * v[i][j_16] - v[i][j_16 + 1];
	}
	for (int i = i_17 + 1; i <= i_18 - 1; i++)
	{
		u[i][j_17 + 1] = 2.0 * u[i][j_17] - u[i][j_17 - 1];
		u[i][j_18 - 1] = 2.0 * u[i][j_18] - u[i][j_18 + 1];
		v[i][j_17 + 1] = 2.0 * v[i][j_17] - v[i][j_17 - 1];
		v[i][j_18 - 1] = 2.0 * v[i][j_18] - v[i][j_18+ 1];
	}
	for (int i = 2; i <= mx - 1; i++)
	{
		for (int j = 2; j <= my - 1; j++)
		{
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6) || (i_7 < i && i < i_8 && j_7 < j && j < j_8) || (i_9 < i && i < i_10 && j_9 < j && j < j_10) || (i_11 < i && i < i_12 && j_11 < j && j < j_12) || (i_13 < i && i < i_14 && j_13 < j && j < j_14) || (i_15 < i && i < i_16 && j_15 < j && j < j_16) || (i_17 < i && i < i_18 && j_17 < j && j < j_18))
			{
				continue;
			}

			urhs[i][j] = urhs[i][j] -
			             v[i][j] *
			             (-u[i][j + 2] + 8.0 * (u[i][j + 1] - u[i][j - 1]) +
			              u[i][j - 2]) /
			             (12.0 * dy) -
			             abs(v[i][j]) *
			             (u[i][j + 2] - 4.0 * u[i][j + 1] + 6.0 * u[i][j] -
			              4.0 * u[i][j - 1] + u[i][j - 2]) /
			             (4.0 * dy);
			vrhs[i][j] = vrhs[i][j] -
			             v[i][j] *
			             (-v[i][j + 2] + 8.0 * (v[i][j + 1] - v[i][j - 1]) +
			              v[i][j - 2]) /
			             (12.0 * dy) -
			             abs(v[i][j]) *
			             (v[i][j + 2] - 4.0 * v[i][j + 1] + 6.0 * v[i][j] -
			              4.0 * v[i][j - 1] + v[i][j - 2]) /
			             (4.0 * dy);
		}
	}

	for (int i = 2; i <= mx - 1; i++)
	{
		for (int j = 2; j <= my - 1; j++)
		{
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6) || (i_7 < i && i < i_8 && j_7 < j && j < j_8) || (i_9 < i && i < i_10 && j_9 < j && j < j_10) || (i_11 < i && i < i_12 && j_11 < j && j < j_12) || (i_13 < i && i < i_14 && j_13 < j && j < j_14) || (i_15 < i && i < i_16 && j_15 < j && j < j_16) || (i_17 < i && i < i_18 && j_17 < j && j < j_18))
			{
				continue;
			}
			u[i][j] = u[i][j] + dt * urhs[i][j];
			v[i][j] = v[i][j] + dt * vrhs[i][j];
		}
	}
}

int file_write_data(double x[mx + 1], double y[my + 1],
                    float p[mx + 2][my + 2])
{
	FILE *fp = fopen("data9.csv", "w");
	for (int i = 1; i <= mx; i++)
	{
		for (int j = 1; j <= my; j++)
		{
			fprintf(fp, "%f %f %f\n", x[i], y[j], p[i][j]);
		}
		fprintf(fp, "\n");
	}
	return 0;
}

int main()
{
	make_xygrid(x, y);
	init_condition(u, v, p);
	bcfor_p(p);
	bcfor_v(u, v);

	std::chrono::system_clock::time_point start_time,
	                                      end_time; //時間計測用変数を確保
	start_time = std::chrono::system_clock::now(); //計測スタート:
	for (int n = 1; n <= nlast; n++)
	{
		poisson_eq(u, v, p, dx, dy, dt);
		bcfor_p(p);
		velocity_eq(u, v, p, dx, dy, dt);
		bcfor_v(u, v);
		end_time = std::chrono::system_clock::now(); //
		double elapsed_time =
			std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
			                                                      start_time)
			.count() / //ここまでにかかった時間を計算
			1000.0; //[s]
		std::cout << "step:" << n << "/" << nlast
		          << " elapsed_time:" << elapsed_time << "[s]" << std::endl;
	}

	file_write_data(x, y, p);
	return 0;
}