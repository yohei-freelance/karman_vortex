#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

const double Re = 70.0;  // Reynolds Number
const double cfl = 0.2;  // CFL Number

/*SOR Pamameters*/
const double omegap = 1.00;
const int maxitp = 100;
const double errorp = 0.00001;

/* No. of Time Steps*/
const int nlast = 10000;

/* set x-grid parameters*/
const int mx = 401;   // x軸格子点数(1~401),x=30.0がmx=401に対応
const int i_2 = 106;   // x軸格子点数基準での,角柱の左端
const int i_3 = 96;  // x軸格子点数基準での,角柱の右端
const int i_1 = 96;
const int i_4 = 106;
const int i_5 = 116;
const int i_6 = 126;

/* set y-grid parameters*/
const int my = 201;   // y軸格子点数(1~201)
const int j_2 = 96;   // y軸格子点数基準での,角柱の下端
const int j_3 = 106;  // y軸格子点数基準での,角柱の上端
const int j_1 = 86;
const int j_4 = 116;
const int j_5 = 96;
const int j_6 = 106;

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
}

void init_condition(float u[mx + 2][my + 2], float v[mx + 2][my + 2],
                    float p[mx + 2][my + 2])
{
	for (int i = 1; i <= mx; i++)
	{
		for (int j = 1; j <= my; j++)
		{
			u[i][j] = 0.7;
			v[i][j] = 0.0;
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
	for (int i = i_1; i <= i_2; i++)
	{
		for (int j = j_1; j <= j_2; j++)
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
}

void poisson_eq(float u[mx + 2][my + 2], float v[mx + 2][my + 2],
                float p[mx + 2][my + 2], double dx, double dy, double dt)
{
	double rhs[mx + 2][my + 2];
	for (int i = 2; i <= mx - 1; i++)
	{
		for (int j = 2; j <= my - 1; j++)
		{
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6)) 
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
				if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6))
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
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6))
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
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6))
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

	for (int i = 2; i <= mx - 1; i++)
	{
		for (int j = 2; j <= my - 1; j++)
		{
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6))
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
	for (int i = 2; i <= mx - 1; i++)
	{
		for (int j = 2; j <= my - 1; j++)
		{
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6))
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
			if ((i_1 < i && i < i_2 && j_1 < j && j < j_2) || (i_3 < i && i < i_4 && j_3 < j && j < j_4) || (i_5 < i && i < i_6 && j_5 < j && j < j_6))
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
	FILE *fp = fopen("data5.csv", "w");
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
