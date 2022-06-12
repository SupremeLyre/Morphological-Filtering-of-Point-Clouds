#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

const int N = 2E6 + 5;
const double inf = 0x3f3f3f3f;

//四至
double l = inf;
double b = inf;
double r = 0;
double up = 0;

//窗口尺寸
int m, n;
int cellsize = 1;

//最大窗口
int winsize = 2;

struct Point
{
    double x, y, z;
};
vector<Point> p;
void read()
{
    ifstream file("test.txt", ios_base::in);
    while (!file.eof())
    {
        double x, y, z;
        file >> x >> y >> z;
        p.push_back({x, y, z});
    }
    file.clear();
    file.close();
}

void get_bounder(double &left, double &right, double &up, double &bottom)
{

    for (auto i : p)
    {
        left = min(left, i.x);
        bottom = min(bottom, i.y);
        right = max(right, i.x);
        up = max(up, i.y);
    }
}
void set_window(int &m, int &n, int &cellsize)
{
    m = (int)((up - b) / cellsize) + 1;
    n = (int)((r - l) / cellsize) + 1;
}

vector<Point> push_window(double &left, double &right, double &up, double &bottom)
{
    vector<Point> res(m * n);
    for (int i = 0; i < m * n; i++)
    {
        res[i].x = 0;
        res[i].y = 0;
        res[i].z = 0;
    }
    for (auto i : p)
    {
        int col = (int)((i.x - left) / cellsize);
        int row = (int)((i.y - bottom) / cellsize);
        if (res[row * m + col].z != 0)
        {
            if (res[row * m + col].z > i.z)
            {
                res[row * m + col].x = i.x;
                res[row * m + col].y = i.y;
                res[row * m + col].z = i.z;
            }
        }
        else if (res[row * m + col].z == 0)
        {
            res[row * m + col].x = i.x;
            res[row * m + col].y = i.y;
            res[row * m + col].z = i.z;
        }
    }
    return res;
}

vector<double> Erosion(vector<Point> &PointGrid) //腐蚀
{
    vector<double> res(m * n);
    for (int i = 0; i < n; i++) // row
    {
        for (int j = 0; j < m; j++) // col
        {
            int midrow = i;
            int midcol = j;
            double Minz = PointGrid[i * m + j].z;
            for (int x = -winsize >> 1; x <= winsize >> 1; x++)
            {
                for (int y = -winsize >> 1; y <= winsize >> 1; y++)
                {
                    if (midrow + x < 0 || midrow + x >= n || midcol + y < 0 || midcol + y >= m)
                        continue;
                    if (PointGrid[(midrow + x) * m + (midcol + y)].z < Minz && PointGrid[(midrow + x) * m + (midcol + y)].z != 0)
                    {
                        Minz = PointGrid[(midrow + x) * m + (midcol + y)].z;
                    }
                }
            }
            res[i * m + j] = Minz;
        }
    }
    return res;
}

vector<double> Dilation(vector<double> &MinPointGrid) //膨胀
{
    vector<double> res(m * n);
    for (int i = 0; i < n; i++) // row
    {
        for (int j = 0; j < m; j++) // col
        {
            int midrow = i;
            int midcol = j;
            double Maxz = MinPointGrid[i * m + j];
            for (int x = -winsize >> 1; x <= winsize >> 1; x++)
            {
                for (int y = -winsize >> 1; y <= winsize >> 1; y++)
                {
                    if (midrow + x < 0 || midrow + x >= n || midcol + y < 0 || midcol + y >= m)
                        continue;
                    if (MinPointGrid[(midrow + x) * m + (midcol + y)] > Maxz && MinPointGrid[(midrow + x) * m + (midcol + y)] != 0)
                    {
                        Maxz = MinPointGrid[(midrow + x) * m + (midcol + y)];
                    }
                }
            }
            res[i * m + j] = Maxz;
        }
    }
    return res;
}

int main()
{
    read();
    printf("size:%ld\n", p.size());
    get_bounder(l, r, up, b);
    printf("left:%lf right:%lf up:%lf bottom:%lf\n", l, r, up, b);
    set_window(m, n, cellsize);
    printf("m = %d n = %d\n", m, n);
    vector<Point> PointGrid = push_window(l, r, up, b);
    printf("size of PointGrid:%ld\n", PointGrid.size());
    vector<double> MinPointGrid = Erosion(PointGrid);
    printf("MinPointGrid get.\n");
    vector<double> MaxPointGrid = Dilation(MinPointGrid);
    printf("MaxPointGrid get.\n");
    ofstream out;
    out.open("result.txt", ios_base::out);
    int sum = 0;
    for (auto i : p)
    {
        int col = (int)((i.x - l) / cellsize);
        int row = (int)((i.y - b) / cellsize);
        if (fabs(MaxPointGrid[row * m + col] - i.z) < 0.3)
        {
            out << fixed << setprecision(3) << i.x << " " << i.y << " " << i.z << endl;
            sum++;
        }
    }
    out.close();
    if (!sum)
    {
        printf("no point.\n");
    }
    else
    {
        printf("%d points in total.\n", sum);
    }
    return 0;
}